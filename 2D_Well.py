from http.client import SWITCHING_PROTOCOLS
from tokenize import Double
import numpy as np
import matplotlib.pyplot as plt
from tkinter import *
from matplotlib import cm
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
from matplotlib.ticker import MultipleLocator
import math
from scipy import integrate

# Make an example of perturbation theory?
# Make two potential regions, one for the main potential and one for the perurbation
# Then show the amount for the first and second order corrections (Will probably need to use degenerate states?)

# Make it start at energy level n =1
# Make a button to toggle the auto update (either auto update when potential is changed or have a button to generate solutions)
# Normalize solutions
# Add an option to make the draw potential value times some factor of 10 i.e. 10^-6

class Cell():

    # Colormap to pick the colors in the potential grid
    cmap = cm.cool
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=100.0)

    FILLED_COLOR_BG = "green"
    EMPTY_COLOR_BG = "white"
    FILLED_COLOR_BORDER = "green"
    EMPTY_COLOR_BORDER = "black"

    def __init__(self, master, x, y, size, value = 0):
        # Initilizae each cell
        self.master = master
        self.abs = x
        self.ord = y
        self.size= size
        self.value = value
        self.fill= False

    def _switch(self, value):
        # Swtich the cell between being filled or not and change its value
        if self.fill:
            self.value = 0
        else:
            self.value = value
        self.fill= not self.fill

    def draw(self):
        # Draw the cells on the canvas and choose the color from the colormap based on its value
        if self.master != None :
            fill = matplotlib.colors.rgb2hex(Cell.cmap(Cell.norm(self.value)))  #Cell.FILLED_COLOR_BG
            outline = Cell.FILLED_COLOR_BORDER

            if not self.fill:
                fill = Cell.EMPTY_COLOR_BG
                outline = Cell.EMPTY_COLOR_BORDER

            xmin = self.abs * self.size
            xmax = xmin + self.size
            ymin = self.ord * self.size
            ymax = ymin + self.size

            self.master.create_rectangle(xmin, ymin, xmax, ymax, fill = fill, outline = outline)

# Class for the potential grid
class CellGrid(Canvas):
    def __init__(self,master, parent, rowNumber, columnNumber, cellSize, *args, **kwargs):
        # Initiilize the potential grid
        Canvas.__init__(self, master, width = cellSize * columnNumber , height = cellSize * rowNumber, *args, **kwargs)

        self.master = master
        self.parent = parent # Reference of the main application to get the solve function
        self.sizex = columnNumber
        self.sizey = rowNumber
        self.cellSize = cellSize


        self.cells = [] # List of all the cells in the grid
        for row in range(rowNumber):

            line = []
            for column in range(columnNumber):
                line.append(Cell(self, column, row, cellSize))

            self.cells.append(line)

        #memorize the cells that have been modified to avoid many switching of state during mouse motion.
        self.switched = []

        #bind click action
        self.bind("<Button-1>", self.handleMouseClick)  
        #bind moving while clicking
        self.bind("<B1-Motion>", self.handleMouseMotion)
        #bind release button action - clear the memory of modified cells.
        self.bind("<ButtonRelease-1>", self.update_potential)

        self.draw()

    # Draw function that just loops through and draws each cell
    def draw(self):
        for row in self.cells:
            for cell in row:
                cell.draw()

    # Function to handle the coordinates of a mouse event
    def _eventCoords(self, event):
        row = int(event.y / self.cellSize)
        column = int(event.x / self.cellSize)
        return row, column

    # Function to handle a mouse click event
    def handleMouseClick(self, event):
        row, column = self._eventCoords(event)
        cell = self.cells[row][column]
        cell._switch(self.parent.slider.get())
        cell.draw()
        #add the cell to the list of cell switched during the click
        self.switched.append(cell)

    # Function to handle if themouse moves while clicked
    def handleMouseMotion(self, event):
        row, column = self._eventCoords(event)
        cell = self.cells[row][column]

        if cell not in self.switched and cell.fill != self.switched[0].fill:
            cell._switch(self.parent.slider.get())
            cell.draw()
            self.switched.append(cell)
    
    # Function to resolve for a new drawn potential
    def update_potential(self, event):
        self.switched.clear()
        self.parent.solve()
        self.parent.plot()

    # Function to return all the values of the potential in the grid
    def get_values(self):
        values = []
        for line in self.cells:
            temp = []
            for cell in line:
                temp.append(cell.value)
            values.append(temp)
        values.reverse()
        return values
    
    # Function to clear all the values in the grid
    def clear_values(self):
        for line in self.cells:
            for cell in line:
                cell.value = 0
                cell.fill = False
        self.draw()



# Main function for he whole application
class Application():
    def __init__(self, master, sizex, sizey, cellSize, method, dimension):
        # Initialize the function
        self.master = master
        self.resolution = sizex
        self.one_D_potential = np.zeros(self.resolution)
        self.one_D_x = np.linspace(0, sizex, self.resolution)
        self.sizex = sizex
        self.sizey = sizey
        self.cellSize = cellSize
        self.method = method
        self.dimension = dimension
        self.probability = True
        self.cmap = 'magma'
        self.energy_level = 0
        self.degenerate_level = 0
        self.total_degeneracy = 1
        
        self.two_D = True

        #############################
        #####                   #####
        #####    2D Widgets     #####
        #####                   #####
        #############################

        # Make a slider bar for the different possible values to draw
        self.slider = Scale(self.master, from_=0, to=100, orient='horizontal')
        self.slider.grid(column = 0, row = 0, columnspan= 4)

        # Make a text box to input values instead
        self.draw_value = StringVar()
        self.draw_value.trace_add('write', self.input_command)
        self.potential_input = Entry(self.master, textvariable = self.draw_value)
        self.potential_input.grid(column = 1, row = 1, padx = 10, columnspan= 2)
        self.potential_input_label = Label(self.master, text = "Draw Potential")
        self.potential_input_label.grid(column = 0, row = 1, columnspan= 2)

        # Draw the grid
        self.grid = CellGrid(self.master, self, self.sizey, self.sizex, self.cellSize) # Row Num, Column Num, Cell Size
        self.grid.grid(column = 0, row = 3, columnspan=4, padx = 10)

        # Make a textbox to change sizey
        self.sizey_string = StringVar()
        self.sizey_string.trace_add('write', self.change_sizey)
        self.sizey_input = Entry(self.master, textvariable = self.sizey_string)
        self.sizey_label = Label(self.master, text = "Size y")
        self.sizey_input.grid(column = 1, row = 6, columnspan= 2, sticky=W)
        self.sizey_label.grid(column = 0, row = 6, columnspan= 1, sticky=E)


        # Make a text box to change cell size
        self.cellSize_string = StringVar()
        self.cellSize_string.trace_add('write', self.change_cellsize)
        self.cellSize_input = Entry(self.master, textvariable = self.cellSize_string)
        self.cellSize_input_label = Label(self.master, text = "Cell Size")
        self.cellSize_input.grid(column = 1, row = 7, columnspan= 2, sticky=W)
        self.cellSize_input_label.grid(column = 0, row = 7, columnspan= 1, sticky=E)


        # Button to toggle between solving methods
        self.swap_method = Button(self.master, command = self.toggle_method, text = "Swap Stacking Method")
        self.swap_method.grid(column = 0, row = 8, columnspan= 4)


        #############################
        #####                   #####
        #####    1D Widgets     #####
        #####                   #####
        #############################

        # Make a text box to choose number of x points
        self.resolution_string = StringVar()
        self.resolution_string.trace_add('write', self.change_resolution)
        self.resolution_input = Entry(self.master, textvariable = self.resolution_string)
        self.resolution_label = Label(self.master, text = "Resolution")

        # Make a text box to input a function to apply as a potential
        self.potential_func_input = Entry(self.master, width = 20)
        self.potential_func_label = Label(self.master, text = "Function")
        
        # Make a text box to choose left x value
        self.leftx_input = Entry(self.master, width = 5)
        self.leftx_label = Label(self.master, text = "Left")

        # Make a text box to choose right x value
        self.rightx_input = Entry(self.master, width = 5)
        self.rightx_label = Label(self.master, text = "Right")
        
        # Make a text box to choose potential to add
        self.potential_1D_input_float = DoubleVar()
        self.potential_1D_input = Spinbox(self.master, from_ = 0, to = 10, textvariable=self.potential_1D_input_float)
        self.potential_1D_input_label = Label(self.master, text = "Potential")
        
        # Button to apply the new potential changes
        self.apply_potential = Button(self.master, command = self.apply_potential_1D, text = "Apply")

        ###############################
        #####                     #####
        #####     Both Widgets    #####
        #####                     #####
        ###############################
        

        # Button to clear the potential
        self.clear_potent = Button(self.master, command = self.clear_potential, text = "Clear")
        self.clear_potent.grid(column = 0, row = 2, columnspan= 4)

        # Make a text box to change sizex
        self.sizex_string = StringVar()
        self.sizex_string.trace_add('write', self.change_sizex)
        self.sizex_input = Entry(self.master, textvariable = self.sizex_string)
        self.sizex_input.grid(column = 1, row = 5, columnspan= 2, sticky=W)
        self.sizex_label = Label(self.master, text = "Size x")
        self.sizex_label.grid(column = 0, row = 5, columnspan= 1, sticky=E)
        

        # Button to toggle between wavefunctio and probability density
        self.swap_graph = Button(self.master, command = self.toggle_graph, text = "Swap Graph")
        self.swap_graph.grid(column = 0, row = 9, columnspan= 4)

        
        # Make a text box to choose the energy level
        self.energy_level_string = StringVar()
        self.energy_level_string.trace_add('write', self.change_energy_level)
        self.energy_level_input = Entry(self.master, textvariable = self.energy_level_string, width = 3)
        self.energy_level_input.grid(column = 1, row = 10, sticky=E)
        self.energy_level_input_label = Label(self.master, text = "Energy Level:")
        self.energy_level_input_label.grid(column = 0, row = 10, sticky=W)

        self.prev_energy_button = Button(self.master, command = self.prev_energy, text = "< Prev")
        self.prev_energy_button.grid(column = 2, row = 10)
        self.next_energy_button = Button(self.master, command = self.next_energy, text = "Next >")
        self.next_energy_button.grid(column = 3, row = 10)

        # Make a next and prev button to go between degenerate states
        self.degenerate_state_label = Label(self.master, text = "Degeneracy: {}".format(self.total_degeneracy))
        self.degenerate_state_label.grid(column = 0, row = 11)
        self.prev_degen_button = Button(self.master, command = self.prev_degenerate, text = "< Prev")
        self.prev_degen_button.grid(column = 1, row = 11)
        self.current_degenerate_state_label = Label(self.master, text = str(self.degenerate_level + 1))
        self.current_degenerate_state_label.grid(column = 2, row = 11)
        self.next_degen_button = Button(self.master, command = self.next_degenerate, text = "Next >")
        self.next_degen_button.grid(column = 3, row = 11)

        # Button to toggle between solving methods
        self.swap_dimensions = Button(self.master, command = self.toggle_dimension, text = "2D 2 Much")
        self.swap_dimensions.grid(column = 0, row = 12, columnspan= 4)



        ############## MAKING THE GRAPH WIDGET ################
        # This needs to happen here so I can then clear the plot and reuse it later when the potential is changed
        # the figure that will contain the plot
        self.fig = Figure(figsize=(8, 6),
                    dpi = 100)
    
        # adding the subplot
        self.plot1 = self.fig.add_subplot(111)


        # Solve the initalize empy potential so we can graph it
        self.solve()
        grid = self.vector_to_grid(self.evecs[self.energy_level][self.degenerate_level])

        for i in range(len(grid)):
            for j in range(len(grid[i])):
                grid[i][j] *= grid[i][j]

        img = self.plot1.pcolormesh(grid, cmap ='magma', shading = 'gouraud')
        self.plot1.title.set_text("Energy: {}".format( self.evals[self.energy_level]))
        self.clrbar = self.fig.colorbar(img)
        self.plot12 = self.plot1.twinx()
        self.plot12.get_yaxis().set_visible(False)

        # creating the Tkinter canvas
        # containing the Matplotlib figure
        self.plot_canvas = FigureCanvasTkAgg(self.fig,
                                master = self.master)  
        self.plot_canvas.draw()
    
        # placing the canvas on the Tkinter window
        self.plot_canvas.get_tk_widget().grid(column = 5, row = 0, rowspan=20, padx = 10)
    
        # creating the Matplotlib toolbar
        toolbar = NavigationToolbar2Tk(self.plot_canvas,
                                    self.master, pack_toolbar=False)
        toolbar.update()
        toolbar.grid(column = 5, row =21)



    def is_float(self, element):
        try:
            float(element)
            return True
        except ValueError:
            return False


    # Function to toggle the dimension
    def toggle_dimension(self):
        self.two_D = not self.two_D
        if self.two_D:
            
            # Forget all the 1D widgets we don't need
            self.sizex_input.grid_forget()
            self.sizex_label.grid_forget()
            self.potential_func_input.grid_forget()
            self.potential_func_label.grid_forget()
            self.leftx_input.grid_forget()
            self.leftx_label.grid_forget()
            self.rightx_input.grid_forget()
            self.rightx_label.grid_forget()
            self.apply_potential.grid_forget()
            self.potential_1D_input.grid_forget()
            self.potential_1D_input_label.grid_forget()
            self.resolution_input.grid_forget()
            self.resolution_label.grid_forget()


            # Show all of the 2D widget in their correct spot
            # 100% this will come back and bite me in the ass later when I also have to change these positions
            self.slider.grid(column = 0, row = 0, columnspan= 4)
            self.potential_input.grid(column = 1, row = 1, columnspan= 2)
            self.potential_input_label.grid(column = 0, row = 1, columnspan= 2)
            self.grid.grid(column = 0, row = 3, columnspan=4, padx = 10)
            self.sizey_input.grid(column = 1, row = 6, columnspan= 2, sticky=W)
            self.sizey_label.grid(column = 0, row = 6, columnspan= 1, sticky=E)
            self.cellSize_input.grid(column = 1, row = 7, columnspan= 2, sticky=W)
            self.cellSize_input_label.grid(column = 0, row = 7, columnspan= 1, sticky=E)
            self.swap_method.grid(column = 0, row = 8, columnspan= 4)
            self.sizex_input.grid(column = 1, row = 5, columnspan= 2, sticky=W)
            self.sizex_label.grid(column = 0, row = 5, columnspan= 1, sticky=E)

            self.solve()
            self.energy_level = 0
            self.total_degeneracy = len(self.evecs[self.energy_level])
            self.degenerate_state_label.config(text = "Degeneracy: {}".format(self.total_degeneracy))
            self.current_degenerate_state_label.config(text = "1")
            self.degenerate_level = 0

            self.plot(switch = True)

        else:
            # Forget the stuff we don't need for 1D
            self.slider.grid_forget()
            self.potential_input.grid_forget()
            self.potential_input_label.grid_forget()
            self.grid.grid_forget()
            self.sizey_input.grid_forget()
            self.sizey_label.grid_forget()
            self.cellSize_input.grid_forget()
            self.cellSize_input_label.grid_forget()
            self.swap_method.grid_forget()
            
            # Put on the stuff we need for 1 D
            self.resolution_input.grid(column = 1, row = 3, columnspan=2, sticky=W)
            self.resolution_label.grid(column = 0, row = 3, columnspan=1, sticky=E)
            self.sizex_input.grid(column = 1, row = 4, columnspan= 2, sticky=W)
            self.sizex_label.grid(column = 0, row = 4, columnspan= 1, sticky=E)
            self.potential_func_input.grid(column = 1, row =5, columnspan= 3, sticky=W)
            self.potential_func_label.grid(column = 0, row = 5, columnspan= 1, sticky=E)
            self.leftx_input.grid(column = 1, row =6, columnspan= 1, sticky = W)
            self.leftx_label.grid(column = 0, row = 6, columnspan= 1, sticky= E)
            self.rightx_input.grid(column = 2, row =6,  columnspan= 1, sticky=E)
            self.rightx_label.grid(column = 3, row = 6, columnspan= 1, sticky=W)
            self.potential_1D_input.grid(column = 1, row = 7, columnspan=2, sticky=W)
            self.potential_1D_input_label.grid(column=0, row=7, columnspan=1, sticky=E)
            self.apply_potential.grid(column = 0, row = 8, columnspan = 4)

            self.solve()
            self.energy_level = 0
            self.total_degeneracy = len(self.evecs[self.energy_level])
            self.degenerate_state_label.config(text = "Degeneracy: {}".format(self.total_degeneracy))
            self.current_degenerate_state_label.config(text = "1")
            self.degenerate_level = 0

            self.clrbar.remove()
            
            self.plot(switch = True)


    # Function to change the resolution of the 1D graph
    def change_resolution(self, var, index, mode):
        if self.resolution_string.get().isdigit():
            self.resolution = max(2, int(self.resolution_string.get()))
            self.one_D_x = np.linspace(0, sizex, self.resolution)
            self.clear_potential()


    def apply_potential_1D(self):

        if "x" in self.potential_func_input.get():
            for i in range(len(self.one_D_potential)):
                x = self.one_D_x[i]
                self.one_D_potential[i] += eval(self.potential_func_input.get(), {"x": x,'sqrt': math.sqrt, 'pow': math.pow, "sin" : math.sin})

        if self.leftx_input.get().isdigit() and self.rightx_input.get().isdigit() and self.is_float(self.potential_1D_input.get()):
            il = int((int(self.leftx_input.get())/self.sizex)*self.resolution)
            ir = int((int(self.rightx_input.get())/self.sizex)*self.resolution)
            self.one_D_potential[il:ir] += float(self.potential_1D_input.get())

        self.solve()
        self.plot()
        return

    # Function to clear the potential
    def clear_potential(self):
        if self.two_D:
            self.grid.clear_values()
            self.solve()
            self.plot()
        else:
            self.one_D_potential = np.zeros(self.resolution)
            self.solve()
            self.plot()

    # Function to go to the previous degenerate state
    def prev_degenerate(self):
        if self.degenerate_level > 0:
            self.degenerate_level -= 1
            self.current_degenerate_state_label.config(text = str(self.degenerate_level + 1))
            self.plot()
        return

    # Function to go to the next degenerate state
    def next_degenerate(self):
        if self.degenerate_level < self.total_degeneracy - 1:
            self.degenerate_level += 1
            self.current_degenerate_state_label.config(text = str(self.degenerate_level + 1))
            self.plot()
        return

    # Function to go to the previous energy level
    def prev_energy(self):
        if self.energy_level > 0:
            self.energy_level -= 1
            self.energy_level_input.select_clear()
            self.total_degeneracy = len(self.evecs[self.energy_level])
            self.degenerate_state_label.config(text = "Degeneracy: {}".format(self.total_degeneracy))
            self.degenerate_level = 0
            self.plot()
        return

    # Function to go to the next energy level
    def next_energy(self):
        if self.energy_level < len(self.evecs) - 1:
            self.energy_level += 1
            self.energy_level_input.select_clear()
            self.total_degeneracy = len(self.evecs[self.energy_level])
            self.degenerate_state_label.config(text = "Degeneracy: {}".format(self.total_degeneracy))
            self.degenerate_level = 0
            self.plot()
        return

    # Function to toggle the stacking method used to solve
    def toggle_method(self):
        if self.method == 1:
            self.method = 2
        else:
            self.method = 1
        self.solve()
        self.plot()

    # Function to toggle the graph from and to probability density and wavefunction
    def toggle_graph(self):
        self.probability = not self.probability
        if self.probability:
            self.cmap = 'magma'
        else:
            self.cmap = "PiYG"
        self.solve()
        self.plot()

    # Function to change the size of the graph in the x direction
    def change_sizex(self, var, index, mode):
        if self.sizex_string.get().isdigit():
            if self.two_D:
                self.sizex = int(self.sizex_string.get())
                self.grid.destroy()
                self.grid = CellGrid(self.master, self, self.sizey, self.sizex, self.cellSize) # Row Num, Column Num, Cell Size
                self.grid.grid(column = 0, row = 2, columnspan=2)
            else:
                self.sizex = int(self.sizex_string.get())
                self.one_D_x = np.linspace(0, self.sizex, self.resolution)
                self.clear_potential()
                self.solve()
                self.plot()
        return

    # Function to change the size of the graph in the y direction
    def change_sizey(self, var, index, mode):
        if self.sizey_string.get().isdigit():
            self.sizey = int(self.sizey_string.get())
            self.grid.destroy()
            self.grid = CellGrid(self.master, self, self.sizey, self.sizex, self.cellSize) # Row Num, Column Num, Cell Size
            self.grid.grid(column = 0, row = 2, columnspan=2)
        return


    def change_cellsize(self, var, index, mode):
        if self.cellSize_string.get().isdigit():
            self.cellSize = int(self.cellSize_string.get())
            self.grid.destroy()
            self.grid = CellGrid(self.master, self, self.sizey, self.sizex, self.cellSize) # Row Num, Column Num, Cell Size
            self.grid.grid(column = 0, row = 2, columnspan=2)
        return

    # Function to change the energy level you are plotting
    def change_energy_level(self, var, index, mode):
        if self.energy_level_string.get().isdigit():
            self.energy_level = int(self.energy_level_string.get())
            self.total_degeneracy = len(self.evecs[self.energy_level])
            self.degenerate_state_label.config(text = "Degeneracy: {}".format(self.total_degeneracy))
            self.degenerate_level = 0
            self.plot()

    # Command for what to do when someone types in the draw potential input box
    def input_command(self, var, index, mode):
        if self.draw_value.get().isdigit():
            self.slider.set(int(self.draw_value.get()))
        return

    # Function to get the index in a list of elements for some x,y position in a 2D space
    def vector_coords(self, x, y):
        # The method number is the same as in the make_laplcian function and is described more there
        if self.method == 1:
            return x + y*self.sizex
        elif self.method == 2:
            return y + x*self.sizey
        return 0


    # Function to make the Laplcian part of the matrix for the hamiltonian
    #############
    # Still need to be done: make the dx and dy division on the laplcian
    # Also comment the specific loops better
    ################
    def make_laplacian(self):

        # Method will determine the way we arrange the column vector for the points in the 2D case (This does not change the 1D case)
        # Method = 1 Will have all the rows stacked vertically in the column vector
        # Method = 2 Will have all the colmns stacked vertically in the column vector

        laplacian = []

        # Make the Laplcaian in 2 dimensions
        if self.two_D:
            laplacian = np.eye(self.sizex*self.sizey) # Make the Laplcian a diagonal square matrix with the size the area of the well
            laplacian *= -4
            # Make the Laplacian for the method
            for x in range(self.sizex):
                for y in range(self.sizey):
                    if x != self.sizex - 1:
                        laplacian[self.vector_coords(x, y)][self.vector_coords(x+1, y)] += 1
                    if x != 0:
                        laplacian[self.vector_coords(x, y)][self.vector_coords(x-1, y)] += 1
                    if y != self.sizey - 1:
                        laplacian[self.vector_coords(x, y)][self.vector_coords(x, y+1)] += 1
                    if y != 0:
                        laplacian[self.vector_coords(x, y)][self.vector_coords(x, y-1)] += 1

        # Make the Laplican in 1 dimension
        else:
            dx = self.sizex/self.resolution
            laplacian = np.eye(self.resolution) # Make the Laplcian a diagonal square matrix of the size of the number of steps
            laplacian *= -2
            # Loop through every point and set the laplcian
            for x in range(self.resolution):
                if x != 0:
                    laplacian[x][x-1] += 1
                if x != self.resolution - 1:
                    laplacian[x][x+1] += 1

            laplacian = laplacian/dx

        return laplacian

    # Function to make the potential matrix
    # This could be more efficient by not looping through the whole thing to make the potential but instead
    # have the function that detects when a square changes call something to update the potential for that grid point
    # This would need me to always store the potential though instead of just calculating it when I need it
    def make_potential(self):
        potential = []
        # Make the potential matrix for 2 Dimensions
        if self.two_D:
            potential = np.zeros((self.sizex*self.sizey, self.sizex*self.sizey)) # Make a sizex*sizey square matrix
            values = self.grid.get_values() # Get the values in the potential grid
            # Loop through all the values on the grid and append them to their spot in the potential matrix
            for x in range(self.sizex):
                for y in range(self.sizey):
                    potential[self.vector_coords(x, y)][self.vector_coords(x,y)] += values[y][x]

        # Make the potential matrix for 1 Dimension
        else:  
            potential = np.zeros((self.resolution, self.resolution))
            for i in range(self.resolution):
                potential[i][i] = self.one_D_potential[i]

        return potential

    # Function to solve for the potential
    def solve(self, full = False):

        # Generate the Laplacian and potential and make the Hamiltonian
        laplacian = self.make_laplacian()
        potential = self.make_potential()
        hamiltonian = laplacian + potential

        # Use Library to solve for the hamiltonian
        eval, vec = np.linalg.eigh(hamiltonian)
        np.append(eval, 0)
        np.append(vec, [])
        self.evals = []
        self.evecs = []
        temp = []
        # Loop throug the solutions and add them to list of lists based on degeneracy
        for i in range(len(eval) - 1):
            temp.append(vec[:,i])
            if np.abs(eval[i] - eval[i+1]) > 0.0001:
                self.evecs.append(temp)
                self.evals.append(eval[i])
                temp = []

    
    # Funtion to take a vector answer and convert it into a 2D array to be plotted
    def vector_to_grid(self, vec):
        grid = []
        
        # Just loop through the grid and append the first sizex or sizey values depending on the method
        if self.method == 1:
            temp = []
            for i in vec:
                temp.append(i)
                if len(temp) == self.sizex:
                    grid.append(temp)
                    temp = []

        elif self.method == 2:
            temp = []
            for i in vec:
                temp.append(i)
                if len(temp) == self.sizey:
                    grid.append(temp)
                    temp = []
            grid = np.transpose(grid) # Need to transpose the y method one to get it to plot correctly

        return grid

    # Function to plot the grid solution to the potential
    def plot(self, switch = False):

        if self.two_D:
            plot_vec = [] # Make a vector to copy the values over to (basically a deep copy oof)

            # Either get the probability density or the wavefunction directly
            if self.probability:
                for i in range(len(self.evecs[self.energy_level][self.degenerate_level])):
                    plot_vec.append(self.evecs[self.energy_level][self.degenerate_level][i] * self.evecs[self.energy_level][self.degenerate_level][i])
            else :
                for i in range(len(self.evecs[self.energy_level][self.degenerate_level])):
                    plot_vec.append(self.evecs[self.energy_level][self.degenerate_level][i])

            # Clear the current plot and put on the new one with proper axis and colorbar
            self.plot1.cla()
            if not switch:
                self.clrbar.remove()
            if switch:
                self.plot12.get_yaxis().set_visible(False)
                self.plot12.cla()
            self.plot1.set_xlim([0, max(self.sizex-1, self.sizey-1)])
            self.plot1.set_ylim([0, max(self.sizex-1, self.sizey-1)])
            self.plot1.set_title("Energy: {}".format(self.evals[self.energy_level]))
            img = self.plot1.pcolormesh(self.vector_to_grid(plot_vec), cmap = self.cmap, shading = 'gouraud')
            self.clrbar = self.fig.colorbar(img, ax = self.plot1)
            self.plot_canvas.draw()
        else:
            plot_vec = []

            norm_const = 1/np.sqrt(integrate.simps(self.evecs[self.energy_level][self.degenerate_level]**2, self.one_D_x))

            # Either get the probability density or the wavefunction directly
            if self.probability:
                for i in range(len(self.evecs[self.energy_level][self.degenerate_level])):
                    plot_vec.append((norm_const**2)*self.evecs[self.energy_level][self.degenerate_level][i] * self.evecs[self.energy_level][self.degenerate_level][i])
            else :
                for i in range(len(self.evecs[self.energy_level][self.degenerate_level])):
                    plot_vec.append((norm_const)*self.evecs[self.energy_level][self.degenerate_level][i])


            # Plot1 is plotting the potential
            # Plot12 is plotting the wavefunction/probability
            if switch:
                self.plot12.get_yaxis().set_visible(True)
            if not switch:
                self.plot12.cla()
            self.plot1.cla()
            if self.probability:
                self.plot1.set_title("Probability Density")
                self.plot12.set_ylabel("Probability")
            else:
                self.plot1.set_title("Wavefunction")
                self.plot12.set_ylabel("Wavefunction")
            self.plot1.plot(self.one_D_x, self.one_D_potential)
            self.plot1.set_ylabel("Potential")
            self.plot1.set_ylim([0, None])
            self.plot12.plot(self.one_D_x, plot_vec)
            self.plot_canvas.draw()
        return

# if __name__ == "__main__": lmao
if __name__ == "__main__" :

    sizex = 10
    sizey = 10
    dimension = 2
    cellSize = 15

    app = Tk()

    everything = Application(app, sizex, sizey, cellSize, 2, 2)

    app.mainloop()
