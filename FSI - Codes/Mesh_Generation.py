## Importing required modules
import tkinter as tk
from tkinter import *
import math

## Read the input files for inputs
def Read_input(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Initialize variables to store parsed data
    global NDIM         # Dimension of the problem
    global DOMAIN       # Dimension of the domain
    global cell_dim     # Dimension of the cell
    global n_cells      # Number of cells per dimension in a block
    global least_dim    # Least dimension of the cell

    # Calculated Fields
    global block_dim        # Dimension of the block
    global max_lvl          # Max level of the tree
    global coarse_grid_dim  # Number of blocks in coarse grid
    global total_cells      # Total Number of cells (including ghost cells)

    # Iterate through lines to find relevant information
    for line in lines:
        line = line.strip()

        ## General Information
        if line.startswith("Number of Dimensions NDIM"):
            NDIM = int(line.split('=')[1].strip())
        elif line.startswith("Dimension of Domain"):
            domain_str = line.split('=')[1].strip()[1:-1]
            DOMAIN = list(float(num) for num in domain_str.split(','))

        ## Block Information
        elif line.startswith("Dimension of the cell"):
            cell_dim_str = line.split('=')[1].strip()[1:-1]
            cell_dim = list(float(num) for num in cell_dim_str.split(','))
        elif line.startswith("Number of Cells per Dimension of the Block"):
            n_cell_str = line.split('=')[1].strip()[1:-1]
            n_cells = list(int(num) for num in n_cell_str.split(','))

        ## Tree Infortmation
        elif line.startswith("Least Dimension of cells needed"):
            least_dim = float(line.split('=')[1].strip())

    block_dim = list(n_cells[i]*cell_dim[i] for i in range(NDIM))
    max_lvl = math.floor(math.log2(min(cell_dim) / least_dim))
    coarse_grid_dim = list(int(DOMAIN[i]/block_dim[i]) for i in range(NDIM))
    total_cells = (n_cells[0] + 2)*(n_cells[1] + 2)

    return


## Class "TBlock" contains the information of the blocks used in the tree
class TBlock:
    def __init__(self, idx, lvl, position):
        self.idx = idx                  # Index of the Block
        self.lvl = lvl                  # Level at which the block is at
        self.is_refined = False         # Refinement flag
        self.position = position        # Position of center of block
        self.ptr_parent = []            # Index of the parent of the block
        self.ptr_children = []          # List of indeces to children blocks
        self.ptr_neighbor = []          # List of indeces to neighboring blocks
        # Cells within each block inlcuding ghost cells and neighbour info
        self.cells    = [None] * total_cells
        self.cell_nbs = [None] * (2*(n_cells[0] + n_cells[1] + 2))


## Class "TLevel" contains the information of the indeces of the blocks at each level of the tree
class TLevel:
    def __init__(self):
        self.idx_block = []  # Index of the box at that level


## Class "TTree" contains the information of the Tree for BB-AMR
class TTree:

    ## Initialisation Function
    def __init__(self, n_cells, max_lvl, min_lvl):
        self.n_cells = n_cells          # Number of cells per dimension in block
        self.max_lvl = max_lvl          # Maximum refinement level of the tree based on min dim of cell
        self.min_lvl = min_lvl          # Minimum refinement level of the tree 
        self.levels = []                # Dictionary to store levels
        self.blocks = []                # List to store all blocks in the tree
        self.block_next_idx = 0         # ID of the next block to be alloted
        self.max_refine_lvl = max_lvl   # Max level to which the mesh will be refined based on user entry

        ## Creating "TLevel" objects for each level to store IDs of the blocks within that level
        for lvl in range(max_lvl+1):
            self.levels.append(TLevel())

    ## Function to Create a new Block Within the Tree
    def Create_block(self, lvl, position, parent_idx, continuous_cell_numbering):
        idx = self.block_next_idx
        self.levels[lvl].idx_block.append(self.block_next_idx)
        self.blocks.append(TBlock(idx = idx,
                                  lvl = lvl,
                                  position = position))
        self.blocks[idx].ptr_parent = parent_idx
        if continuous_cell_numbering:
            self.blocks[idx].cells = [i for i in range(idx*total_cells, (idx+1)*total_cells)]
        else:
            self.blocks[idx].cells = [i for i in range(total_cells)]
        self.block_next_idx += 1
  
    ## Function to get the neighbor Idx at Coarse Grid
    def Get_coarse_grid_nb(self):
        
        ## Block Neighbors ##

        temp_idx = 0
        for i in range(coarse_grid_dim[0]):
            for j in range(coarse_grid_dim[1]):
                self.blocks[temp_idx].ptr_neighbor = [ temp_idx - coarse_grid_dim[1] - 1,   # left bottom
                                                       temp_idx - coarse_grid_dim[1],       # left middle
                                                       temp_idx - coarse_grid_dim[1] + 1,   # left top
                                                       temp_idx +  1,                       # top middle
                                                       temp_idx + coarse_grid_dim[1] + 1,   # right top
                                                       temp_idx + coarse_grid_dim[1],       # right middle
                                                       temp_idx + coarse_grid_dim[1] - 1,   # right bottom
                                                       temp_idx -  1 ]                      # bottom middle
                
                # Dealing with blocks at boundaries (Blocks which are outside the boundary is given (-1)
                if i == 0:                           # Left boundary
                    self.blocks[temp_idx].ptr_neighbor[:3] = [-1, -1, -1]
                if i == coarse_grid_dim[0] - 1:     # Right boundary
                    self.blocks[temp_idx].ptr_neighbor[4:7] = [-1, -1, -1]
                if j == 0:                          # Bottom boundary
                    self.blocks[temp_idx].ptr_neighbor[6:8] = [-1, -1]
                    self.blocks[temp_idx].ptr_neighbor[0] = -1
                if j == coarse_grid_dim[1] - 1:     # Top boundary
                    self.blocks[temp_idx].ptr_neighbor[2:5] = [-1, -1, -1]

                temp_idx += 1


        ## Cell Neighbors ##
        for idx in range(coarse_grid_dim[0]*coarse_grid_dim[1]):
            # Left Bottom
            nb_idx = self.blocks[idx].ptr_neighbor[0]
            if nb_idx == -1:
                self.blocks[idx].cell_nbs[0] = -1
            else:
                self.blocks[idx].cell_nbs[0] = self.blocks[nb_idx].cells[(n_cells[0]+1)*(n_cells[1]+2)-2]

            ## Left Middle
            nb_idx = self.blocks[idx].ptr_neighbor[1]
            if nb_idx == -1:
                for i in range(1, n_cells[1]+1): self.blocks[idx].cell_nbs[i] = -1
            else: 
                for i in range(1, n_cells[1]+1):
                    self.blocks[idx].cell_nbs[i] = self.blocks[nb_idx].cells[(n_cells[0])*(n_cells[1]+2)+i]

            ## Left Top
            nb_idx = self.blocks[idx].ptr_neighbor[2]
            if nb_idx == -1:
                self.blocks[idx].cell_nbs[n_cells[1]+1] = -1
            else:    
                self.blocks[idx].cell_nbs[n_cells[1]+1] = self.blocks[nb_idx].cells[(n_cells[0])*(n_cells[1]+2)+1]

            ## Top Middle
            nb_idx = self.blocks[idx].ptr_neighbor[3]
            if nb_idx == -1:
                for i in range(1, n_cells[0]+1): self.blocks[idx].cell_nbs[(n_cells[1]+1)+i] = -1
            else: 
                for i in range(1, n_cells[0]+1):
                    self.blocks[idx].cell_nbs[(n_cells[1]+1)+i] = self.blocks[nb_idx].cells[(n_cells[1]+2)*i+1]

            ## Right Top
            nb_idx = self.blocks[idx].ptr_neighbor[4]
            if nb_idx == -1:
                self.blocks[idx].cell_nbs[n_cells[1]+n_cells[0]+2] = -1
            else:
                self.blocks[idx].cell_nbs[n_cells[1]+n_cells[0]+2] = self.blocks[nb_idx].cells[(n_cells[1]+3)]

            ## Right Middle
            nb_idx = self.blocks[idx].ptr_neighbor[5]
            if nb_idx == -1:
                for i in range(1, n_cells[1]+1): self.blocks[idx].cell_nbs[(n_cells[1]+n_cells[0]+2)+i] = -1
            else: 
                for i in range(1, n_cells[1]+1):
                    self.blocks[idx].cell_nbs[(n_cells[1]+n_cells[0]+2)+i] = self.blocks[nb_idx].cells[((n_cells[1]+2)*2-1)-i]

            ## Right Bottom
            nb_idx = self.blocks[idx].ptr_neighbor[6]
            if nb_idx == -1:
                self.blocks[idx].cell_nbs[n_cells[1]*2+n_cells[0]+3] = -1
            else:
                self.blocks[idx].cell_nbs[n_cells[1]*2+n_cells[0]+3] = self.blocks[nb_idx].cells[((n_cells[1]+2)*2-2)]

            ## Bottom Middle
            nb_idx = self.blocks[idx].ptr_neighbor[7]
            if nb_idx == -1:
                for i in range(1, n_cells[0]+1): self.blocks[idx].cell_nbs[(n_cells[1]*2+n_cells[0]+3)+i] = -1
            else: 
                for i in range(1, n_cells[0]+1):
                    self.blocks[idx].cell_nbs[(n_cells[1]*2+n_cells[0]+3)+i] = self.blocks[nb_idx].cells[((n_cells[1]+2)*(n_cells[0]+2)-(n_cells[1]+2)*i)-2]

    ## Function to Create Coarse Grid
    def Create_coarse_grid(self, continuous_cell_numbering):
        if NDIM == 2:
            ## Creating the blocks at coarse Grid
            for i in range(coarse_grid_dim[0]):
                for j in range(coarse_grid_dim[1]):
                    position = list([(0.5+float(i))*float(block_dim[0]), (0.5+float(j))*float(block_dim[1])])
                    self.Create_block(lvl=0, 
                                      position=position, 
                                      parent_idx=None,
                                      continuous_cell_numbering=continuous_cell_numbering)
        
            ## Finding the pointers at coarse Grid
            self.Get_coarse_grid_nb()

        save_idx.append(self.block_next_idx)

    ## Function to get neighbor Idx of child blocks
    def Get_child_block_nb(self, parent_idx):

        #  for ptr_children[0, 1, 2, 3]
        nb_block_0_ = [None] * 8     
        nb_block_1_ = [None] * 8     
        nb_block_2_ = [None] * 8     
        nb_block_3_ = [None] * 8 

        ## Left bottom neighbor    
        nb_idx = self.blocks[parent_idx].ptr_neighbor[0]
        if nb_idx == -1:
            nb_block_0_[0] = -1
        else:
            if self.blocks[nb_idx].is_refined:
                nb_block_0_[0] = self.blocks[nb_idx].ptr_children[3]
            else:
                nb_block_0_[0] = nb_idx
        
        ## Left middle neighbor
        nb_idx = self.blocks[parent_idx].ptr_neighbor[1]
        if nb_idx == -1:
            nb_block_0_[1:3] = [-1] * 2
            nb_block_1_[0:2] = [-1] * 2
        elif self.blocks[nb_idx].is_refined:
            nb_block_0_[1] = self.blocks[nb_idx].ptr_children[2]
            nb_block_0_[2] = self.blocks[nb_idx].ptr_children[3]
            nb_block_1_[0] = self.blocks[nb_idx].ptr_children[2]
            nb_block_1_[1] = self.blocks[nb_idx].ptr_children[3]
        else:
            nb_block_0_[1:3] = [nb_idx] * 2
            nb_block_1_[0:2] = [nb_idx] * 2

        ## Left top neighbor
        nb_idx = self.blocks[parent_idx].ptr_neighbor[2]
        if nb_idx == -1:
            nb_block_1_[2] = -1
        else:
            if self.blocks[nb_idx].is_refined:
                nb_block_1_[2] = self.blocks[nb_idx].ptr_children[2]
            else:
                nb_block_1_[2] = nb_idx

        ## Top Middle neighbor
        nb_idx = self.blocks[parent_idx].ptr_neighbor[3]
        if nb_idx == -1:
            nb_block_1_[3:5] = [-1] * 2
            nb_block_3_[2:4] = [-1] * 2
        elif self.blocks[nb_idx].is_refined:
            nb_block_1_[3] = self.blocks[nb_idx].ptr_children[0]
            nb_block_1_[4] = self.blocks[nb_idx].ptr_children[2]
            nb_block_3_[2] = self.blocks[nb_idx].ptr_children[0]
            nb_block_3_[3] = self.blocks[nb_idx].ptr_children[2]
        else:
            nb_block_1_[3:5] = [nb_idx] * 2
            nb_block_3_[2:4] = [nb_idx] * 2
        
        ## Right top neighbor
        nb_idx = self.blocks[parent_idx].ptr_neighbor[4]
        if nb_idx == -1:
            nb_block_3_[4] = -1
        else:
            if self.blocks[nb_idx].is_refined:
                nb_block_3_[4] = self.blocks[nb_idx].ptr_children[0]
            else:
                nb_block_3_[4] = nb_idx

        ## Right Middle neighbor
        nb_idx = self.blocks[parent_idx].ptr_neighbor[5]
        if nb_idx == -1:
            nb_block_3_[5:7] = [-1] * 2
            nb_block_2_[4:6] = [-1] * 2
        elif self.blocks[nb_idx].is_refined:
            nb_block_3_[5] = self.blocks[nb_idx].ptr_children[1]
            nb_block_3_[6] = self.blocks[nb_idx].ptr_children[0]
            nb_block_2_[4] = self.blocks[nb_idx].ptr_children[1]
            nb_block_2_[5] = self.blocks[nb_idx].ptr_children[0]
        else:
            nb_block_3_[5:7] = [nb_idx] * 2
            nb_block_2_[4:6] = [nb_idx] * 2

        ## Right bottom neighbor
        nb_idx = self.blocks[parent_idx].ptr_neighbor[6]
        if nb_idx == -1:
            nb_block_2_[6] = -1
        else:
            if self.blocks[nb_idx].is_refined:
                nb_block_2_[6] = self.blocks[nb_idx].ptr_children[1]
            else:
                nb_block_2_[6] = nb_idx


        ## Bottom Middle neighbor
        nb_idx = self.blocks[parent_idx].ptr_neighbor[7]
        if nb_idx == -1:
            nb_block_2_[7] = -1
            nb_block_2_[0] = -1
            nb_block_0_[6:8] = [-1] * 2
        elif self.blocks[nb_idx].is_refined:
            nb_block_2_[7] = self.blocks[nb_idx].ptr_children[3]
            nb_block_2_[0] = self.blocks[nb_idx].ptr_children[1]
            nb_block_0_[6] = self.blocks[nb_idx].ptr_children[3]
            nb_block_0_[7] = self.blocks[nb_idx].ptr_children[1]
        else:
            nb_block_2_[7] = nb_idx
            nb_block_2_[0] = nb_idx
            nb_block_0_[6:8] = [nb_idx] * 2

        ## Internal neighbors
        nb_block_0_[3] = self.blocks[parent_idx].ptr_children[1]
        nb_block_0_[4] = self.blocks[parent_idx].ptr_children[3]
        nb_block_0_[5] = self.blocks[parent_idx].ptr_children[2]

        nb_block_1_[5] = self.blocks[parent_idx].ptr_children[3]
        nb_block_1_[6] = self.blocks[parent_idx].ptr_children[2]
        nb_block_1_[7] = self.blocks[parent_idx].ptr_children[0]

        nb_block_2_[1] = self.blocks[parent_idx].ptr_children[0]
        nb_block_2_[2] = self.blocks[parent_idx].ptr_children[1]
        nb_block_2_[3] = self.blocks[parent_idx].ptr_children[3]

        nb_block_3_[7] = self.blocks[parent_idx].ptr_children[2]
        nb_block_3_[0] = self.blocks[parent_idx].ptr_children[0]
        nb_block_3_[1] = self.blocks[parent_idx].ptr_children[1]

        ## Setting the neighbours of the child blocks
        child_idx = self.blocks[parent_idx].ptr_children
        self.blocks[child_idx[0]].ptr_neighbor = nb_block_0_
        self.blocks[child_idx[1]].ptr_neighbor = nb_block_1_
        self.blocks[child_idx[2]].ptr_neighbor = nb_block_2_
        self.blocks[child_idx[3]].ptr_neighbor = nb_block_3_

    ## Function to create Child Blocks after refinement
    def Create_child_blocks(self, idx, continuous_cell_numbering):

        lvl = self.blocks[idx].lvl         # Level of the parent
        if lvl >= self.max_refine_lvl or self.blocks[idx].is_refined:
            return
        
        self.blocks[idx].is_refined = True

        if NDIM == 2:

            ## Create the blocks
            for i in range(2):
                for j in range(2):
                    center = self.blocks[idx].position
                    height = block_dim[1]/4/(2**lvl)
                    width  = block_dim[0]/4/(2**lvl)
                    position = [center[0] - width + i*width*2, center[1] - height + j*height*2]
                    self.Create_block(lvl+1, position, idx, continuous_cell_numbering)
                    self.blocks[idx].ptr_children = [i for i in range(self.block_next_idx-(2**NDIM), self.block_next_idx)]

            ## Get the nieghbors for the child blocks and its cells
            self.Get_child_block_nb(idx)
            for idx_child in self.blocks[idx].ptr_children:
                self.Get_cell_nb(idx_child)

            ## Update neighbor's "ptr_neighbor" list
            for nb_idx in self.blocks[idx].ptr_neighbor:
                if self.blocks[nb_idx].is_refined:
                    self.Get_child_block_nb(nb_idx)
                    for idx_child in self.blocks[nb_idx].ptr_children:
                        self.Get_cell_nb(idx_child)

            ## Refine Neighbor block to maintain 2:1 scale
            for nb_idx in self.blocks[idx].ptr_neighbor:
                if nb_idx < 0:
                    continue
                if self.blocks[nb_idx].is_refined:
                    continue
                elif (self.blocks[idx].lvl - self.blocks[nb_idx].lvl) >= 1:
                    print(f'\tThe Block [{nb_idx}] is refined now due to [{idx}]')
                    self.Create_child_blocks(nb_idx, continuous_cell_numbering)
                self.Get_child_block_nb(idx)
                for idx_child in self.blocks[idx].ptr_children:
                    self.Get_cell_nb(idx_child)

    ## Function to get cell neighbours of every block
    def Get_cell_nb(self, idx):

        ## Cell Neighbors default (neighbors assumed at same level)
        
        # Left Bottom
        nb_idx = self.blocks[idx].ptr_neighbor[0]
        if nb_idx == -1:
            self.blocks[idx].cell_nbs[0] = -1
        else:
            self.blocks[idx].cell_nbs[0] = self.blocks[nb_idx].cells[(n_cells[0]+1)*(n_cells[1]+2)-2]

        ## Left Middle
        nb_idx = self.blocks[idx].ptr_neighbor[1]
        if nb_idx == -1:
            for i in range(1, n_cells[1]+1): self.blocks[idx].cell_nbs[i] = -1
        else: 
            for i in range(1, n_cells[1]+1):
                self.blocks[idx].cell_nbs[i] = self.blocks[nb_idx].cells[(n_cells[0])*(n_cells[1]+2)+i]

        ## Left Top
        nb_idx = self.blocks[idx].ptr_neighbor[2]
        if nb_idx == -1:
            self.blocks[idx].cell_nbs[n_cells[1]+1] = -1
        else:    
            self.blocks[idx].cell_nbs[n_cells[1]+1] = self.blocks[nb_idx].cells[(n_cells[0])*(n_cells[1]+2)+1]

        ## Top Middle
        nb_idx = self.blocks[idx].ptr_neighbor[3]
        if nb_idx == -1:
            for i in range(1, n_cells[0]+1): self.blocks[idx].cell_nbs[(n_cells[1]+1)+i] = -1
        else: 
            for i in range(1, n_cells[0]+1):
                self.blocks[idx].cell_nbs[(n_cells[1]+1)+i] = self.blocks[nb_idx].cells[(n_cells[1]+2)*i+1]

        ## Right Top
        nb_idx = self.blocks[idx].ptr_neighbor[4]
        if nb_idx == -1:
            self.blocks[idx].cell_nbs[n_cells[1]+n_cells[0]+2] = -1
        else:
            self.blocks[idx].cell_nbs[n_cells[1]+n_cells[0]+2] = self.blocks[nb_idx].cells[(n_cells[1]+3)]

        ## Right Middle
        nb_idx = self.blocks[idx].ptr_neighbor[5]
        if nb_idx == -1:
            for i in range(1, n_cells[1]+1): self.blocks[idx].cell_nbs[(n_cells[1]+n_cells[0]+2)+i] = -1
        else: 
            for i in range(1, n_cells[1]+1):
                self.blocks[idx].cell_nbs[(n_cells[1]+n_cells[0]+2)+i] = self.blocks[nb_idx].cells[((n_cells[1]+2)*2-1)-i]

        ## Right Bottom
        nb_idx = self.blocks[idx].ptr_neighbor[6]
        if nb_idx == -1:
            self.blocks[idx].cell_nbs[n_cells[1]*2+n_cells[0]+3] = -1
        else:
            self.blocks[idx].cell_nbs[n_cells[1]*2+n_cells[0]+3] = self.blocks[nb_idx].cells[((n_cells[1]+2)*2-2)]

        ## Bottom Middle
        nb_idx = self.blocks[idx].ptr_neighbor[7]
        if nb_idx == -1:
            for i in range(1, n_cells[0]+1): self.blocks[idx].cell_nbs[(n_cells[1]*2+n_cells[0]+3)+i] = -1
        else: 
            for i in range(1, n_cells[0]+1):
                self.blocks[idx].cell_nbs[(n_cells[1]*2+n_cells[0]+3)+i] = self.blocks[nb_idx].cells[((n_cells[1]+2)*(n_cells[0]+2)-(n_cells[1]+2)*i)-2]



        ## Cell Neighbors (when neighbor is at higher level)

        ## Case 1
        nb_1 = self.blocks[idx].ptr_neighbor[0]
        nb_2 = self.blocks[idx].ptr_neighbor[1]
        if (nb_1 == nb_2) and (nb_1 != -1):
            self.blocks[idx].cell_nbs[0] = self.blocks[nb_1].cells[(n_cells[1]+2)*(n_cells[0])+int(0.5*(n_cells[1]))]
            for i in range(1, n_cells[1]+1):
                self.blocks[idx].cell_nbs[i] = self.blocks[nb_1].cells[(n_cells[1]+2)*(n_cells[0])+int(0.5*(n_cells[1])+1)+int((i-0.5)/2)]

        ## Case 2
        nb_1 = self.blocks[idx].ptr_neighbor[1]
        nb_2 = self.blocks[idx].ptr_neighbor[2]
        if (nb_1 == nb_2) and (nb_1 != -1):
            for i in range(1, n_cells[1]+1):
                self.blocks[idx].cell_nbs[i] = self.blocks[nb_1].cells[(n_cells[1]+2)*(n_cells[0]) + 1 + int((i-0.5)/2)]
            self.blocks[idx].cell_nbs[n_cells[1]+1] = self.blocks[nb_1].cells[(n_cells[1]+2)*(n_cells[0])+int(0.5*(n_cells[1])+1)]

        ## Case 3
        nb_1 = self.blocks[idx].ptr_neighbor[2]
        nb_2 = self.blocks[idx].ptr_neighbor[3]
        if (nb_1 == nb_2) and (nb_1 != -1):
            self.blocks[idx].cell_nbs[n_cells[1]+1] = self.blocks[nb_1].cells[(n_cells[1]+2)*int(n_cells[0]/2)+1]
            for i in range(1, n_cells[0]+1):
                self.blocks[idx].cell_nbs[n_cells[1]+1+i] = self.blocks[nb_1].cells[(n_cells[1]+2)*int(n_cells[0]/2)+1+int(i/2+0.5)*(n_cells[1]+2)]

        ## Case 4
        nb_1 = self.blocks[idx].ptr_neighbor[3]
        nb_2 = self.blocks[idx].ptr_neighbor[4]
        if (nb_1 == nb_2) and (nb_1 != -1):
            for i in range(1, n_cells[0]+1):
                self.blocks[idx].cell_nbs[n_cells[1]+1+i] = self.blocks[nb_1].cells[int(i/2+0.5)*(n_cells[1]+2)+1]
            self.blocks[idx].cell_nbs[n_cells[1]+n_cells[0]+2] = self.blocks[nb_1].cells[(n_cells[1]+2)*int(n_cells[0]/2+1)+1]

        ## Case 5
        nb_1 = self.blocks[idx].ptr_neighbor[4]
        nb_2 = self.blocks[idx].ptr_neighbor[5]
        if (nb_1 == nb_2) and (nb_1 != -1):
            self.blocks[idx].cell_nbs[n_cells[1]+n_cells[0]+2] = self.blocks[nb_1].cells[(n_cells[1]+2)+int(n_cells[1]/2)+1]
            for i in range(1, n_cells[1]+1):
                self.blocks[idx].cell_nbs[n_cells[1]+n_cells[0]+2+i] = self.blocks[nb_1].cells[int((n_cells[1]+2)*1.5)-int(i/2+0.5)]

        ## Case 6
        nb_1 = self.blocks[idx].ptr_neighbor[5]
        nb_2 = self.blocks[idx].ptr_neighbor[6]
        if (nb_1 == nb_2) and (nb_1 != -1):
            for i in range(1, n_cells[1]+1):
                self.blocks[idx].cell_nbs[n_cells[1]+n_cells[0]+2+i] = self.blocks[nb_1].cells[int((n_cells[1]+2)*2)-1-int(i/2+0.5)]
            self.blocks[idx].cell_nbs[2*n_cells[1]+n_cells[0]+3] = self.blocks[nb_1].cells[int((n_cells[1]+2)*1.5)-1]

        ## Case 7
        nb_1 = self.blocks[idx].ptr_neighbor[6]
        nb_2 = self.blocks[idx].ptr_neighbor[7]
        if (nb_1 == nb_2) and (nb_1 != -1):
            self.blocks[idx].cell_nbs[2*n_cells[1]+n_cells[0]+3] = self.blocks[nb_1].cells[(n_cells[1]+2)*int(n_cells[0]/2+2)-2]
            for i in range(1, n_cells[0]+1):
                self.blocks[idx].cell_nbs[2*n_cells[1]+n_cells[0]+3+i] = self.blocks[nb_1].cells[(n_cells[1]+2)*int(n_cells[0]/2+2)-(n_cells[1]+2)*int(i/2+0.5)-2]

        ## Case 8
        nb_1 = self.blocks[idx].ptr_neighbor[7]
        nb_2 = self.blocks[idx].ptr_neighbor[0]
        if (nb_1 == nb_2) and (nb_1 != -1):
            for i in range(1, n_cells[0]+1):
                self.blocks[idx].cell_nbs[2*n_cells[1]+n_cells[0]+3+i] = self.blocks[nb_1].cells[(n_cells[1]+2)*int(n_cells[0]+2)-(n_cells[1]+2)*int(i/2+0.5)-2]
            self.blocks[idx].cell_nbs[0] = self.blocks[nb_1].cells[int((n_cells[1]+2)*(n_cells[0]/2+1))-2]

        return


## Creating the GUI for Meshing
def Create_Mesh_GUI(tree, display_id, display_cell, continuous_cell_numbering):
    
    ## Create root window for the GUI
    def Create_root_window():
        root = tk.Tk()
        root.title("IIT-G FSI Solver  [ BLOCK BASED MESHING ]")
        root.iconbitmap('./IITG_logo.ico')
        root.configure(bg='grey')
        return root

    ## Function to get the window dimension
    def Get_window_dimensions(window_scale):
        # Get screen width and height
        screen_width = root.winfo_screenwidth()
        screen_height = root.winfo_screenheight()
        
        # Set the window size to be half of the screen size
        window_width = int(screen_width * window_scale)
        window_height = int(screen_height * window_scale)
        
        # Calculate the position to center the window
        position_right = int(screen_width/2 - window_width/2)
        position_down = int(screen_height/2 - window_height/2)
        
        # Set the geometry of the window
        return window_width, window_height, position_right, position_down

    ## Function to get canvas dimension
    def Get_canvas_dimension(canvas_scale):
        mesh_width = int(window_width * canvas_scale)
        mesh_height = int(window_height * canvas_scale)
        
        # Calculate the position to center the canvas
        canvas_x = (window_width - float(mesh_width)) / 2
        canvas_y = (window_height - float(mesh_height)) / 2

        return mesh_width, mesh_height, canvas_x, canvas_y

    ## Function to get offset and position scale
    def Get_pos_scale_offset(domain_scale):
        domain_width, domain_height = mesh_width*domain_scale, mesh_height*domain_scale
        pos_scale = min(domain_width/DOMAIN[0], domain_height/DOMAIN[1])
        off_x = (mesh_width - DOMAIN[0]*pos_scale)/2
        off_y = (mesh_height - DOMAIN[1]*pos_scale)/2
        return pos_scale, [off_x, off_y]

    ## Function to output x0, y0, x1, y1 for "create_rectangle" of tk.Canvas
    def Get_block_rectangle_info(offset, position, width, height, pos_scale, dim_scale, mesh_height):
            block_center    = list(pos_scale*position[i] + offset[i] for i in range(NDIM))
            block_width     = dim_scale*width
            block_height    = dim_scale*height
            x0 = block_center[0] - block_width/2
            y0 = mesh_height - (block_center[1] - block_height/2)
            x1 = block_center[0] + block_width/2
            y1 = mesh_height - (block_center[1] + block_height/2)
            return x0, y0, x1, y1

    ## Function to display the blocks as a mesh in the Canvas
    def Display_mesh(display_id, display_cell):

        mesh.create_rectangle(0, 0, mesh_width, mesh_height, fill="grey", outline="grey")

        for lvl in range(max_lvl+1):
            dim_scale = pos_scale/2**lvl
            if math.log(max(coarse_grid_dim)) == 0:
                line_width = 1
            else: line_width = 4/math.log(max(coarse_grid_dim))/(1.5)**lvl
            
            for idx in tree.levels[lvl].idx_block:
                position = tree.blocks[idx].position
                width, height = block_dim[0], block_dim[1]

                x0, y0, x1, y1 = Get_block_rectangle_info(offset, position, width, height, pos_scale, dim_scale, mesh_height)
                mesh.create_rectangle(x0, y0, x1, y1, 
                                    width=line_width)
                if display_id:
                    mesh.create_text((x0+x1)/2, (y0+y1)/2, 
                                     text=idx, 
                                     font=("Purisa", int(line_width*8)))

                if display_cell:
                    for i in range(n_cells[0]):
                        for j in range(n_cells[1]):

                            xc0 = x0 + i*(x1-x0)/n_cells[0] +2
                            yc0 = y0 + j*(y1-y0)/n_cells[1] -2
                            xc1 = x0 + (i+1)*(x1-x0)/n_cells[0] -2
                            yc1 = y0 + (j+1)*(y1-y0)/n_cells[1] +2

                            mesh.create_rectangle(xc0, yc0, xc1, yc1, 
                                                width=line_width/4, outline="#5A5A5A")
                            if display_id:
                                mesh.create_text((xc0+xc1)/2, (yc0+yc1)/2, 
                                                 text=tree.blocks[idx].cells[n_cells[1]*i+j + (n_cells[1]+3+2*i)], 
                                                 font=("Purisa", int(line_width*6)), 
                                                 fill="white")

    ## Function to draw new blocks
    def Draw_blocks(last_id):

        for id in range(last_id, tree.block_next_idx):
            lvl = tree.blocks[id].lvl
            dim_scale = pos_scale/2**lvl
            if math.log(max(coarse_grid_dim)) == 0:
                line_width = 1
            else: line_width = 4/math.log(max(coarse_grid_dim))/(1.5)**lvl
            position = tree.blocks[id].position
            width, height = block_dim[0], block_dim[1]

            x0, y0, x1, y1 = Get_block_rectangle_info(offset, position, width, height, pos_scale, dim_scale, mesh_height)
            mesh.create_rectangle(x0, y0, x1, y1, 
                                width=line_width)
            
            if display_id:
                mesh.create_text((x0+x1)/2, (y0+y1)/2, text=id, font=("Purisa", int(line_width*8)))
                
            if display_cell:
                for i in range(n_cells[0]):
                    for j in range(n_cells[1]):

                        xc0 = x0 + i*(x1-x0)/n_cells[0] +2/lvl
                        yc0 = y0 + j*(y1-y0)/n_cells[1] -2/lvl
                        xc1 = x0 + (i+1)*(x1-x0)/n_cells[0] -2/lvl
                        yc1 = y0 + (j+1)*(y1-y0)/n_cells[1] +2/lvl

                        mesh.create_rectangle(xc0, yc0, xc1, yc1, 
                                            width=line_width/4, outline="#5A5A5A")
                        if display_id:
                            mesh.create_text((xc0+xc1)/2, (yc0+yc1)/2, 
                                             text=tree.blocks[id].cells[n_cells[1]*i+j + (n_cells[1]+3+2*i)], 
                                             font=("Purisa", int(line_width*6)),
                                             fill="white")

    ## Function to print values just to check the code
    def Print_tree(tree):

        print(f'\n ----- Indices at levels ----- \n')
        for lvl in range(max_lvl+1):
            if len(tree.levels[lvl].idx_block)>0:
                print(f'Level {lvl}: {tree.levels[lvl].idx_block}')

        print(f'\n ----- Child Pointers ----- \n')
        for lvl in range(max_lvl+1):
            for id in tree.levels[lvl].idx_block:
                if tree.blocks[id].is_refined:
                    print(f'ID: {id},    children: {tree.blocks[id].ptr_children}')

        print(f'\n ----- Parent Pointers ----- \n') 
        for lvl in range(max_lvl+1):
            for id in tree.levels[lvl].idx_block:
                if tree.blocks[id].ptr_parent != None:
                    print(f'ID: {id},    Parent: {tree.blocks[id].ptr_parent}')

        print(f'\n ----- Neighbour Pointers ----- \n') 
        for lvl in range(max_lvl+1):
            for id in tree.levels[lvl].idx_block:
                print(f'LVL: {lvl},    ID: {id},    Neighbour: {tree.blocks[id].ptr_neighbor}')

    ## Function to get block ID based on (x,y) of click/drag
    def Get_id_from_xy(x, y):
        x = (x - offset[0])/pos_scale
        y = (mesh_height - offset[1] - y )/pos_scale
        is_outside = 0

        # Checking if the click is outside the grid
        if not (x>0 and y>0 and x<DOMAIN[0] and y<DOMAIN[1]):
            is_outside = 1
        
        # Find the block in coarse grid where the mouse is
        i = int(x/block_dim[0])
        j = int(y/block_dim[1])
        refine_id = i*coarse_grid_dim[1] + j

        return refine_id, is_outside, x, y

    ## Function to refine the mesh on clicks
    def Refine_mesh_click():

        def on_click(event):

            # Get the position in domain coordinates
            x, y = event.x, event.y
            refine_id, is_outside, x, y = Get_id_from_xy(x, y)
            if is_outside:
                return

            # If the grid is refined, find the grid which is clicked within
            # to be refined.
            while (tree.blocks[refine_id].is_refined):
                lvl = tree.blocks[refine_id].lvl + 1
                for check_id in tree.blocks[refine_id].ptr_children:
                    check_block_pos = tree.blocks[check_id].position
                    check_block_dim = [i/2**lvl for i in block_dim]

                    # conditions to check if the click is inside block
                    cond_left   = (x >= (check_block_pos[0] - check_block_dim[0]/2))
                    cond_right  = (x <= (check_block_pos[0] + check_block_dim[0]/2))
                    cond_bottom = (y >= (check_block_pos[1] - check_block_dim[1]/2))
                    cond_top    = (y <= (check_block_pos[1] + check_block_dim[1]/2))

                    if cond_left and cond_right and cond_top and cond_bottom:
                        refine_id = check_id


            # Refine the block to be refined
            if tree.blocks[refine_id].lvl < tree.max_refine_lvl:
                last_id = tree.block_next_idx
                print(f'\n\tThe Block [{refine_id}] is refined now due to [click]')
                tree.Create_child_blocks(refine_id, continuous_cell_numbering)
                b_save()
                Draw_blocks(last_id)

        
        mesh.bind("<Button-1>", on_click)

    ## Function to refine the mesh on right-click-drag
    def Refine_mesh_right_click_drag():

        current_id = [-1]
        
        def on_drag(event):

            nonlocal current_id

            # Get the position in domain coordinates
            x, y = event.x, event.y
            refine_id, is_outside, x, y = Get_id_from_xy(x, y)
            if is_outside:
                return

            # If the grid is refined, find the grid which is clicked within
            # to be refined.
            while (tree.blocks[refine_id].is_refined):
                lvl = tree.blocks[refine_id].lvl + 1
                for check_id in tree.blocks[refine_id].ptr_children:
                    check_block_pos = tree.blocks[check_id].position
                    check_block_dim = [i/2**lvl for i in block_dim]

                    # conditions to check if the click is inside block
                    cond_left   = (x >= (check_block_pos[0] - check_block_dim[0]/2))
                    cond_right  = (x <= (check_block_pos[0] + check_block_dim[0]/2))
                    cond_bottom = (y >= (check_block_pos[1] - check_block_dim[1]/2))
                    cond_top    = (y <= (check_block_pos[1] + check_block_dim[1]/2))

                    if cond_left and cond_right and cond_top and cond_bottom:
                        refine_id = check_id


            # Checking when mouse is dragged to new block
            if current_id[0] != refine_id:

                # Refine the block to be refined
                if ((current_id[0] > 0) and 
                    not(tree.blocks[current_id[0]].is_refined) and 
                    (tree.blocks[current_id[0]].lvl<tree.max_refine_lvl)):
                    
                    last_id = tree.block_next_idx
                    print(f'\n\tThe Block [{current_id[0]}] is refined now due to [R-click-drag]')
                    tree.Create_child_blocks(current_id[0], continuous_cell_numbering)
                    Draw_blocks(last_id)
                current_id[0] = refine_id

            else:
                return
        
        def on_release(event):

            b_save()
            nonlocal current_id
            current_id = [-1]

        mesh.bind("<B3-Motion>", on_drag)
        mesh.bind("<ButtonRelease-3>", on_release)
   
    ## Function to refine a selected region left-click-drag       ## Need to work on this
    def Refine_mesh_left_click_drag():
        start_x, start_y = None, None
        drag_box = None

        def on_press(event):
            nonlocal start_x, start_y
            start_x, start_y = event.x, event.y

        def on_drag(event):
            nonlocal start_x, start_y, drag_box
            if start_x is None or start_y is None:
                return
            if drag_box is not None:
                mesh.delete(drag_box)
            drag_box = mesh.create_rectangle(start_x, start_y, event.x, event.y, 
                                             outline='black', fill='lightblue', stipple='gray25')

        def on_release(event):
            nonlocal start_x, start_y, drag_box
            if start_x is None or start_y is None:
                return
            
            if drag_box is not None:
                mesh.delete(drag_box)
            
            end_x, end_y = event.x, event.y
            
            x_min, x_max = min(start_x, end_x), max(start_x, end_x)
            y_min, y_max = min(start_y, end_y), max(start_y, end_y)

            x_min, x_max = (x_min - offset[0])/pos_scale, (x_max - offset[0])/pos_scale
            y_max, y_min = (mesh_height - offset[1] - y_min)/pos_scale, (mesh_height - offset[1] - y_max)/pos_scale
            # print(f'xmin, xmax: {x_min}, {x_max}\nymin, ymax: {y_min}, {y_max}')

            for lvl in range(tree.max_refine_lvl):
                check_block_dim = [i/2**lvl for i in block_dim]
                for id in tree.levels[lvl].idx_block:
                    check_block_pos = tree.blocks[id].position
                    cond_1 = (check_block_pos[0] + check_block_dim[0]/2 <= x_min)
                    cond_2 = (check_block_pos[0] - check_block_dim[0]/2 >= x_max)
                    cond_3 = (check_block_pos[1] + check_block_dim[1]/2 <= y_min)
                    cond_4 = (check_block_pos[1] - check_block_dim[1]/2 >= y_max)

                    if not(cond_1 or cond_2 or cond_3 or cond_4):
                        if not tree.blocks[id].is_refined:
                            last_id = tree.block_next_idx
                            print(f'\n\tThe Block [{id}] is refined now due to [L-click-drag]')
                            tree.Create_child_blocks(id, continuous_cell_numbering)
                            Draw_blocks(last_id)
            
            # Reset start position
            start_x, start_y = None, None
            b_save()

        # Bind mouse events
        mesh.bind("<Button-1>", on_press)
        mesh.bind("<B1-Motion>", on_drag)
        mesh.bind("<ButtonRelease-1>", on_release)

    ## Function to call refine on click using button
    def b_click_refine():
        print("\n-----  Refine on Click  -----")
        Refine_mesh_click()
        ## Refine the Mesh on L-Click-Drag as well
        Refine_mesh_right_click_drag()

    ## Function to call refine on drag using button
    def b_drag_refine():
        print("\n-----  Refine on Drag  -----")
        Refine_mesh_left_click_drag()
        ## Refine the Mesh on L-Click-Drag as well
        Refine_mesh_right_click_drag()

    ## Function to get max refine level information from GUI
    def b_max_lvl_value():
        value = entry_max_lvl.get()
        if int(value) > tree.max_lvl:
            value = tree.max_lvl
        tree.max_refine_lvl = int(value)
        print(f'\n-----  Max_refine_to: {value}  -----')
        entry_max_lvl.delete(0, tk.END)
        entry_max_lvl.insert(0, str(value))

    ## Function to display ID for the blocks and cells
    def b_display_block_id():
        global is_on_block_id
        nonlocal display_id
        if is_on_block_id:
            print(f'\n---- Displaying without ID ----')
            display_id = False
            is_on_block_id = False
            Display_mesh(display_id, display_cell)
            
        else:
            print(f'\n---- Displaying with ID ----')
            display_id = True
            is_on_block_id = True
            Display_mesh(display_id, display_cell)

    ## Function to display ID for the blocks and cells
    def b_display_cells():
        global is_on_cell_display
        nonlocal display_cell
        if is_on_cell_display:
            print(f'\n---- Displaying without Cells ----')
            display_cell = False
            is_on_cell_display = False
            Display_mesh(display_id, display_cell)
            
        else:
            print(f'\n---- Displaying with Cells ----')
            display_cell = True
            is_on_cell_display = True
            Display_mesh(display_id, display_cell)

    ## Function to call display tree details on click
    def b_print_details():
        print(f'\n---- Tree Details ----')
        Print_tree(tree)

    ## Function to save a state of mesh on click
    def b_save():
        print(f'\n---- Saving Mesh ----')
        if len(save_idx) > 0:
            if tree.block_next_idx > save_idx[-1]:
                save_idx.append(tree.block_next_idx)
        else: 
            save_idx.append(tree.block_next_idx)
        print(f'Saved IDs: {save_idx}')

    ## Function to save a state of mesh on click
    def b_undo():
        if len(save_idx) > 0:
            print(f'\n---- Reloading Mesh at Last Save ----')

            last_save_id = save_idx[-1]

            did_undo = False

            # Removing the blocks
            for id_rem in range(last_save_id, tree.block_next_idx):
                did_undo = True
                del(tree.blocks[last_save_id])

            # Removing ids in the level information
            for lvl in range(tree.max_lvl+1):
                ids = tree.levels[lvl].idx_block.copy()
                for id in ids:
                    if id >= last_save_id:
                        tree.levels[lvl].idx_block.remove(id)
            
            # Removing child pointers and refinement flag based on condition
            for lvl in range(max_lvl+1):
                for idx in tree.levels[lvl].idx_block:
                    if len(tree.blocks[idx].ptr_children)>0:
                        if tree.blocks[idx].ptr_children[0]>=last_save_id:
                            tree.blocks[idx].ptr_children = []
                            tree.blocks[idx].is_refined = False

            for lvl in range(1, max_lvl+1):
                for idx in tree.levels[lvl].idx_block:
                    tree.Get_child_block_nb(tree.blocks[idx].ptr_parent)
                    tree.Get_cell_nb(idx)

            
            tree.block_next_idx = save_idx[-1]
            if not did_undo:
                save_idx.pop()
                b_undo()
            mesh.delete("all")
            Display_mesh(display_id, display_cell)
            print(f'Saved IDs: {save_idx}')
        
        if len(save_idx) == 0:
            save_idx.append(tree.block_next_idx)

    ## Function to save a state of mesh on click
    def b_verify():
        left_x   = 0.01*window_width
        top_y    = 0.1*mesh_height
        v_space  = 25
        print(f'\n---- Verifying ----')
        idx = int(entry_verify.get())
        label_block = tk.Label(root, text="Block Neighbor: " + str([i for i in tree.blocks[idx].ptr_neighbor]) + " "*100, bg="grey")
        label_cell = tk.Label(root, text="Cell Neighbor:    "+ str(tree.blocks[idx].cell_nbs) + " "*100, bg="grey")



        label_block.place (x=left_x + 200, y=top_y+ 20*v_space)
        label_cell.place (x=left_x + 200, y=top_y+ 21*v_space)

    ## Function Carrying all the buttons
    def Buttons():

        ## Buttons for selecting method of refinement
        label_refine = tk.Label(root, text=" ------ Refinement ------ ", font=("Arial", 10, "bold"), bg="grey")
        button_max_lvl = tk.Button(root, text="Max Level to refine", command=b_max_lvl_value, width=20)
        button_click = tk.Button(root, text="Refine on Click", command=b_click_refine, width=20)
        button_drag = tk.Button(root, text="Refine on Drag", command=b_drag_refine, width=20)
        
        label_display = tk.Label(root, text=" ------ Display ------ ", font=("Arial", 10, "bold"), bg="grey")
        button_display_block_id = tk.Button(root, text="Display Block ID", command=b_display_block_id, width=20)
        button_display_cells = tk.Button(root, text="Display Cells", command=b_display_cells, width=20)
        
        label_print = tk.Label(root, text=" ------ Print ------ ", font=("Arial", 10, "bold"), bg="grey")
        button_print = tk.Button(root, text="Print Details", command=b_print_details, width=20)

        label_save_undo = tk.Label(root, text=" ------ Save and Undo ------ ", font=("Arial", 10, "bold"), bg="grey")
        button_undo = tk.Button(root, text="Undo", command=b_undo, width=20)

        label_verify = tk.Label(root, text=" ------ Verify Block ------ ", font=("Arial", 10, "bold"), bg="grey")
        button_verify = tk.Button(root, text="Verify", command=b_verify, width=20)


        # Placing the buttons on the window
        left_x   = 0.01*window_width
        right_x  = 0.9*window_width - 50
        top_y    = 0.1*mesh_height
        v_space  = 25

        label_refine.place  (x=left_x, y=top_y)
        entry_max_lvl.place (x=left_x, y=top_y+ 1*v_space)
        button_max_lvl.place(x=left_x, y=top_y+ 2*v_space)
        button_click.place  (x=left_x, y=top_y+ 4*v_space)
        button_drag.place   (x=left_x, y=top_y+ 6*v_space)

        label_display.place          (x=left_x, y=top_y+ 9*v_space)
        button_display_block_id.place(x=left_x, y=top_y+ 10*v_space)
        button_display_cells.place   (x=left_x, y=top_y+ 12*v_space)

        label_print.place (x=left_x, y=top_y+ 15*v_space)
        button_print.place(x=left_x, y=top_y+ 16*v_space)

        label_verify.place (x=left_x, y=top_y+ 19*v_space)
        entry_verify.place (x=left_x, y=top_y+ 20*v_space)
        button_verify.place(x=left_x, y=top_y+ 21*v_space)

        label_save_undo.place(x=right_x, y=top_y)
        button_undo.place    (x=right_x, y=top_y+ 1*v_space)


    # Create the Root Window and Canvas Mesh of TKINTER and set the title and icon
    root = Create_root_window()

    # Set the geometry of the window
    window_width, window_height, position_right, position_down = Get_window_dimensions(0.8)
    root.geometry(f"{window_width}x{window_height}+{position_right}+{position_down}")

    # Create the canvas for the mesh
    mesh_width, mesh_height, canvas_x, canvas_y = Get_canvas_dimension(0.8)
    mesh = tk.Canvas(root, width=mesh_width, height=mesh_height, highlightthickness=0)
    mesh.place(x=canvas_x, y=canvas_y+100, anchor=CENTER)
    mesh.pack()

    # Create the Blocks
    pos_scale, offset = Get_pos_scale_offset(0.8)

    ## Show all the blocks on the Canvas
    Display_mesh(display_id, display_cell)

    ## Add buttons to the Window with its functions
    entry_max_lvl = tk.Entry(root, width=20)
    entry_max_lvl.insert(0, str(tree.max_lvl))
    entry_verify = tk.Entry(root, width=20)
    entry_verify.insert(0, str(4))
    Buttons()
    
    ## Run the GUI
    root.mainloop()


## Main Program
def main():

    ## Reading the input file
    file_path = 'Input.txt'
    Read_input(file_path)

    ## Creating a Tree object as an instance of "TTree" class
    tree = TTree(n_cells = n_cells,
                 max_lvl = max_lvl,
                 min_lvl = 0)
    

    ## Switch for the buttons of block_id and cell display
    global is_on_block_id
    global is_on_cell_display
    is_on_block_id     = False
    is_on_cell_display = False

    continuous_cell_numbering = True

    ## Saving State of Mesh
    global save_idx
    save_idx = []

    ## Creating the coarse grid
    tree.Create_coarse_grid(continuous_cell_numbering=continuous_cell_numbering)


    ## GUI of mesh formed
    Create_Mesh_GUI(tree=tree, 
                    display_id=False, 
                    display_cell=False,
                    continuous_cell_numbering=continuous_cell_numbering)


## Calling the main function
if __name__ == '__main__':
    main()
    print("Program ran successfully")