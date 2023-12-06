import numpy as np
from tkinter import *
from PIL import Image, ImageFilter, ImageTk as itk

import time
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from tqdm import tqdm 
#Creates a sizexsize system with an iceblock covering "percent_ice" of the system
# This one creates the initial temperatures of ice and water
def initialize_system(size, percent_ice):
    
    initial_water_temp = 2
    initial_ice_temp = -31
    
    system = np.ones((size, size)) * (273 + initial_water_temp)
    ice = np.zeros((size, size))
    
    # Calculate the range for the square area
    start_row = round(size / 2 - (size / (partial_icebit)))
    end_row = round(size / 2 + (size / (partial_icebit)))
    start_col = round(size / 2 - (size /(partial_icebit)))
    end_col = round(size / 2 + (size / (partial_icebit)))
    
    # Set the lower temperature in the square area
    system[start_row:end_row, start_col:end_col] = 273 + initial_ice_temp # HÄR HAR VI ÄNDRAT
    
    # Set 1s in the square area
    ice[start_row:end_row, start_col:end_col] = 1

    return system, ice

# needs to be scaled so that the daily energytransfer matches the conduction. 
def temperature_changes(maxG, minG, t):
    return (maxG-np.abs(((maxG-minG))*np.cos( np.pi * t*2/Max_time)))*(60*60*24) # HÄR HAR VI ÄNDRAT

# For a given time step each cell is periodically interacts with it's neighbors.
# The interaction: all neighbor blocks are summed into 1 block and the energy-function
# over the heatexchange is used to find their common temperature after the exchange.
def update_temperature(T, ice, solarenergy, time):
    directions = np.array([[-1, 0], [0, -1], [0, 1], [1, 0]])
    row_indices, col_indices = len(T), len(T[0])
    
    T_new = np.zeros_like(T)
    ice_new = np.zeros_like(ice)
    solarenergy = temperature_changes(maxG, minG, time)

    for row in range(row_indices):
        for col in range(col_indices):

            dr = (row + directions[:, 0]) % N
            dc = (col + directions[:, 1]) % N

            CPm=np.zeros(4)
            rad=np.zeros(4)
            for i in range(4):
                # Calculate the coordinates of the neighbor
                neighbor_row = dr[i]
                neighbor_col = dc[i]
                
                # Check the temperature at the neighbor
                if T[neighbor_row, neighbor_col] <= 273:
                    CPm[i] = Cp_ice * m_ice
                    rad[i] = solarenergy*absorption_ice #absorption of ice
                else:
                    CPm[i] = Cp_water * m_water
                    rad[i] = solarenergy*absorption_water #absorption of water

            # Calculate the average value for neighbors
            avg_neighbour_rad = np.sum(rad) / 4
            avg_neighbour_mCp = np.sum(CPm) / 4
            avg_neighbour_temp = np.sum(T[dr, dc]) / 4

            # Chooses the correct radb, m, Cp and Q_phase
            #----
            # For now avg T, Cp, m
            #----
            
            if avg_neighbour_temp > 273 and T[row, col] <= 273: #ice and water interaction -->  melting
                Q_phase = m_ice * Lf
                radb = solarenergy*absorption_ice #absorption of ice
                mb, Cpb = m_ice, Cp_ice
            elif avg_neighbour_temp <= 273 and T[row, col] > 273: #ice and water interaction -->  freezing
                radb = solarenergy*absorption_water #absorption of water
                Q_phase = -m_water * Lf
                mb, Cpb = m_water, Cp_water
            elif avg_neighbour_temp and T[row, col] <= 273: #both are ice --> no melting
                radb = solarenergy*absorption_ice #absorption of ice
                Q_phase = 0
                mb, Cpb = m_ice, Cp_ice
            elif avg_neighbour_temp and T[row, col] > 273: #both are water --> no freezing
                radb = solarenergy*absorption_water #absorption of water
                Q_phase = 0
                mb, Cpb = m_water, Cp_water
                          
            cell=Cpb * mb * T[row, col]+ radb- Q_phase
            avg_neighbour=avg_neighbour_temp * avg_neighbour_mCp+avg_neighbour_rad

            Tf = (cell + avg_neighbour) / (Cpb * mb + avg_neighbour_mCp)
            
            radb, avg_neighbour_rad= 0,0 # DENNA ÄR NERFLYTTAT  

            T_new[row, col] = Tf
            ice_new[row, col] = 0 if Tf > 273 else 1

    return T_new, ice_new, solarenergy

def animate_simulation():
    global T, ice, amp, f, solarenergy
    
    ice_coverage = []

    for time in tqdm(range(Max_time)):
        T, ice, solarenergy = update_temperature(T, ice, solarenergy, time)
        
        ice_coverage.append(np.sum(ice)/(N*N))

        # Set colors for ice and water
        ice_color = 'red'  # Set the color for ice (white)
        water_color = 'blue'  # Set the color for water (blue)

        # Create a custom colormap
        cmap = LinearSegmentedColormap.from_list('ice_water', ['blue', 'green', 'yellow', 'red' ], N=256)
        if time == 0 or time ==100 or time == 200 or time == 360:
            # Plot the temperature grid
            plt.imshow(T, cmap=cmap, interpolation='nearest')
            plt.axis('off')
            plt.title(f'Temperature at t={time}')
            plt.colorbar()
            plt.savefig(f'TempConfiguration{time}')
            plt.show()

        ice_im[ice ==1] = 255
        ice_im[ice ==0] = [0,0,255]
        ice_image = Image.fromarray(np.uint8(ice_im), 'RGB').resize((n, n), Image.NEAREST)
        image_sharp = ice_image.filter(ImageFilter.SHARPEN)
        img_sharp2 = itk.PhotoImage(image_sharp)
        img = itk.PhotoImage(ice_image)
        canvas.create_image(0, 0, anchor=NW, image=img)
        window.title(f"Ice block with seasonal solar energy = {solarenergy} at time step {time}")
        window.update()
        #time.sleep(0.5)
        
    x_val = list(range(Max_time))
    plt.plot(x_val, ice_coverage)
    plt.ylabel('Ice coverage [%]')
    plt.xlabel('Time step [days]')
    plt.title('Ice coverage over a year')
    plt.savefig('IceCoverage')
    plt.show()

# Simulation parameters
N = 100
Max_time = 361
n = 500
S = 1364 # DET SVÄNGDE GANSKA SNYGGT OM AMP = 104
maxG = 204+100 
minG = 17
absorption_ice, absorption_water = 0.01, 0.5
solarenergy = 273
Cp_water, Cp_ice, m_ice, m_water, Lf = 4186, 2108, 917, 1000, 334*10^3

# Create the window for animation
window = Toplevel()
window.geometry(str(int(n*1.2)) + "x" + str(int(n*1.2)))
window.configure(background='white')
canvas = Canvas(window, width=500, height=500, bd=2)
canvas.pack()

# Create initial T and ice block
partial_icebit = 4
percent_ice = 4

T, ice = initialize_system(N,percent_ice)
#ice = initialize_ice(N,partial_icebit)

# Set up colormap
cmap = LinearSegmentedColormap.from_list('rg', ["w", "blue"], N=256)

# Create image arrays for visualization
ice_im = np.zeros((N, N, 3))

# Animate the simulation
animate_simulation()

# Main loop for Tkinter window
window.mainloop()
#%%

