# =============================================================================   
# Importing
# =============================================================================   

import numpy as np
import model_params as model
import scrubbers_canistors as sg
import scipy.sparse as sp
from workspace import build_matrix
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from itertools import combinations

from alphashape import alphashape
from shapely.geometry import Point  

parameters = model.parameters()
  
# =============================================================================   
# Function
# =============================================================================   
 
def max_distance(parameters: model.parameters, gas_canisters: list, build_matrix):
    """
    Since each gas canister in the room is assumed to be identical, 
    (i.e idential initial concentration, gas, dimensions)
    the function finds the maximum distance from each canister a detector can be placed
    such that it detects a concentration value of 1% of the initial concentration released.
    
    It then calculates the optimal location the detector should be place such that 
    it is within this distance for each gas canister in the room.
    
    If it is not possible to find an optimal location for every canister, 
    the function finds the next best location to cover as many as it can.

    Args:
        parameters        Uses parameters defined in the model to compute the large block Matrix
            Nx_points     number of x spatial points
            Ny_points     number of y spatial points
            etc...
        
            U             array storing concentration values of the spatial grid at each time step
        
    Outputs:
        det_loc           a tuple array showing the location the detector should be placed in the room.
                          (length, width)
                            
    """

    # =========================================================================  
    # Finding convex/concave hull for each canister
    # =========================================================================  
    
    # Initialise storage lists for plotting
    gas_hulls = []
    canister_hulls = []
    satisfied_points = []
    canister_centres = []
    
    for canister in gas_canisters:
        
    # ------------------------------------------------------------------------- 
    # Running the model for each canister
    # -------------------------------------------------------------------------  
    # Perhaps turn model in to a function so we can just call the function here?      
    
        x = np.linspace(0, parameters.len_x, parameters.Nx_points)
        y = np.linspace(0, parameters.len_y, parameters.Ny_points)
        X, Y = np.meshgrid(x, y)

        A = build_matrix(parameters=parameters, scrubbers=scrubbers, x=x, y=y)
        U = np.zeros((parameters.Nt_points, parameters.Nx_points * parameters.Ny_points))

        u_init = canister.get_initial_concentration(x, y).reshape((parameters.Nx_points * parameters.Ny_points),)
        u = u_init.copy()
        U[0, :] = u_init

        for n in range(1, parameters.Nt_points):
            u = sp.linalg.minres(A, u)[0]
            U[n, :] = u
        
        # Unpack the flattened U matrix into a 3D matrix
        # ( with dimensions Nx_points, Ny_points, Nt_points )
        U_3D = np.zeros((parameters.Nx_points, parameters.Ny_points, parameters.Nt_points))
        for i in range(parameters.Nx_points):
            U_3D[i] = U[:, i::parameters.Nx_points].T
        
        # ---------------------------------------------------------------------
        # Finding points which satisfy 1% of initial conc. 
        # ---------------------------------------------------------------------       
        
        # Get threshold value for detector 
        threshold = 0.01 * canister.concentration  # 1% threshold
        
        # Get canister location
        x_loc = canister.x_loc
        y_loc = canister.y_loc
         
        # Find indices of canister location in the spatial grid
        # (closest value in the x and y arrays to x_loc and y_loc)
        x_idx = (np.abs(x - x_loc)).argmin()  
        y_idx = (np.abs(y - y_loc)).argmin()  
        
        # Store this location
        canister_centres.append(np.array([x_idx, y_idx]))

        # Find points where concentration is at least as much as the threshold
        # this is across all timesteps to return a 2d array
        threshold_mask = np.any(U_3D >= threshold, axis=2) # Boolean
        points = np.argwhere(threshold_mask)  # Extract the "True" grid points
        
        # Store Boolean of points which satisfy this condition for the canister
        satisfied_points.append(threshold_mask)
        
        # ---------------------------------------------------------------------
        # Plot for checking 
        # ---------------------------------------------------------------------  
        
        # Plot to check if the grid points are correct
        # fig, ax = plt.subplots()
        # for point in points:
        #     x, y = point  # Unpack the tuple into x and y
        #     ax.plot(x, y, 'o', color='cyan', markersize=1)
        # Commented out by default as only used for checking
            
        # ---------------------------------------------------------------------
        # Compute concave/convex hulls for gas and cannister 
        # ---------------------------------------------------------------------  
            
        # Compute a hull for gas using alpha shapes
        alpha = 0.1  # Lower alpha = tighter hull
        if len(points) >= 4:  # Alpha shapes require at least 4 points
            gas_hull = alphashape(points, alpha)
            if gas_hull.is_empty:
                print(f"Hull could not be formed for canister at ({canister.x_loc}, {canister.y_loc})")
            else:
                # Store the computed hull in storage list
                gas_hulls.append(gas_hull)

        # Compute a hull for canister using alpha shapes
        caniser_mask = canister.get_initial_concentration(x, y)
        caniser_points = np.argwhere(caniser_mask)  
        caniser_shape = alphashape(caniser_points, alpha)
        canister_hulls.append(caniser_shape)    
                
    
    # =========================================================================  
    # Finding optimal location for detector(s)
    # =========================================================================  
              
    # Generate the index combinations of all our canisters in descending order
    num_canisters = len(gas_canisters)
    total_combs = []
    for r in range(num_canisters, 0, -1):
        total_combs.extend(combinations(range(num_canisters), r))
    
    # Now we figure out how many detectors we need
    # We start with needing 1 detector
    detectors_needed = 1
    
   # Generate the range from num_canisters to 1 (inclusive)
   # Iterate over each combination layer
   # Top layer = 1 detector needed
   # Second layer = 2 needed etc...
    for i in range(num_canisters, 0, -1):
        
        current_combs = [comb for comb in total_combs if len(comb) == i]
        satisfied_combs = [] # Storage for index combinations that have overlap
        
        # Check if any of the combinations at this level can be satisfied
        for comb in current_combs:
            combined_mask = np.logical_and.reduce([satisfied_points[idx] for idx in comb])
            if np.any(combined_mask):
                # If combination is satisfied, add to list
                satisfied_combs.append(comb)
            
        if satisfied_combs:
            break # We found the minimum number of detectors needed at this level
            
        # Otherwise no combination at this level is satisfied, go to next level
        else:
            detectors_needed += 1
            
    # We now need to find what canisters are not included in the overlapping indices

    # Flatten the tuples to find all unique indices that are included in satisfied_combs
    flattened_satisfied_combs = set(index for comb in satisfied_combs for index in comb)
    
    # Find missing elements by checking which indices are not in the satisfied combinations
    missing_cans = [can for can in range(num_canisters) if can not in flattened_satisfied_combs]
    
    # Add each missing element as a single-element tuple to satisfied_combs
    satisfied_combs.extend((elem,) for elem in missing_cans)

    # Now that we have found a layer aka number of detectors we need, then:     
    detector_locs = []
    
    for comb in satisfied_combs:
        
        # Create boolean mask for this index combination
        combined_mask = np.logical_and.reduce([satisfied_points[idx] for idx in comb])
        
        # Extract the "True" grid points
        candidate_points = np.argwhere(combined_mask)
    
        # Function to calculate distance between two points
        def distance(point1, point2):
            return np.linalg.norm(point1 - point2)
    
        # Find the candidate point that maximizes the total distance to the centers of the canisters in the combination
        best_points = []
        max_total_distance = -np.inf

        for candidate in candidate_points:
            total_distance = sum(distance(candidate, canister_centres[idx]) for idx in comb)
            if total_distance > max_total_distance:
                max_total_distance = total_distance
                best_points = [candidate]  # Start a new list
  
            elif total_distance == max_total_distance:
                best_points.append(candidate)  # Append to the list of best points
        
        # If there are multiple best locations, we use a tie-breaking condition:
        # Select the location closest to the other canisters not included
        # within this index combination
        
        if len(best_points) > 1:
            min_distance_sum = np.inf
            best_point = None
            
            # Get all indices not in comb 
            remaining_indices = [i for i in range(num_canisters) if i not in comb]
            
            # Get all centers of the other canisters 
            canister_centres_filtered = [canister_centres[i] for i in range(len(canister_centres)) if i not in remaining_indices]
            
            if not canister_centres_filtered:
                canister_centres_filtered = canister_centres
                
            # Find the best point by minimizing the sum of distances to all other canisters
            for candidate in best_points:
                distance_sum = sum(distance(candidate, centre) for centre in canister_centres_filtered)
                
                # We also check for equality to just choose random
                # point if minimized distance is same
                if distance_sum <= min_distance_sum:
                    min_distance_sum = distance_sum
                    best_point = candidate  

        else:
            best_point = best_points[0]

        
        detector_locs.append(best_point)   
        
    # =========================================================================  
    # Plotting and Output
    # ========================================================================= 
    
    # Choose a random canister to leak
    canister = np.random.choice(gas_canisters)
    
    # Running the model for random canister
    x = np.linspace(0, parameters.len_x, parameters.Nx_points)
    y = np.linspace(0, parameters.len_y, parameters.Ny_points)
    X, Y = np.meshgrid(x, y)

    A = build_matrix(parameters=parameters, scrubbers=scrubbers, x=x, y=y)
    U = np.zeros((parameters.Nt_points, parameters.Nx_points * parameters.Ny_points))

    u_init = canister.get_initial_concentration(x, y).reshape((parameters.Nx_points * parameters.Ny_points),)
    u = u_init.copy()
    U[0, :] = u_init

    for n in range(1, parameters.Nt_points):
        u = sp.linalg.minres(A, u)[0]
        U[n, :] = u
    
    # Unpack the flattened U matrix into a 3D matrix
    # ( with dimensions Nx_points, Ny_points, Nt_points )
    U_3D = np.zeros((parameters.Nx_points, parameters.Ny_points, parameters.Nt_points))
    for i in range(parameters.Nx_points):
        U_3D[i] = U[:, i::parameters.Nx_points].T
        
    # -------------------------------------------------------------------------
    
    # Initialising plot
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_aspect('equal')
    ax.set_xlim(0, parameters.Nx_points)
    ax.set_ylim(0, parameters.Ny_points)
    ax.set_xlabel("X Grid Index")
    ax.set_ylabel("Y Grid Index")
    ax.set_title("Concave Hulls and Canister Locations")

    # -------------------------------------------------------------------------
        
    # Plot the heatmap
    min_concentration = 0
    max_concentration = 0.2
    
    heatmap = ax.imshow(U_3D[:, :, 0], cmap='magma', vmin=min_concentration, vmax=max_concentration, interpolation='nearest')

    # Animation function
    def update(frame):
        heatmap.set_array(U_3D[:, :, frame])
        ax.set_title(f"Time Step: {frame}")
        return heatmap,
    
    global ani
    ani = FuncAnimation(fig, update, frames=parameters.Nt_points, interval=100, blit=False, repeat=True)
    plt.colorbar(heatmap, ax=ax)

    # -------------------------------------------------------------------------
    
    # Iterate through each gas hull and plot it
    for gas_hull in gas_hulls:
        if not gas_hull.is_empty:
            hull_x, hull_y = gas_hull.exterior.xy  
            ax.plot(hull_y, hull_x, color='red', linestyle='-', linewidth=2, label="Max Distance to detect 1%")
            
    # -------------------------------------------------------------------------
            
    # Iterate through each canister hull and plot it
    for canister_shape in canister_hulls:
        if not canister_shape.is_empty:
            can_x, can_y = canister_shape.exterior.xy  
            ax.plot(can_x, can_y, color='blue', linestyle='-', linewidth=2, label="Gas Canister Circumference")
            
    # -------------------------------------------------------------------------
    
    # Plot the locations of the gas canisters on the grid
    for canister in gas_canisters:
        # Convert canister's physical coordinates to grid indices
        x_idx = (np.abs(x - canister.x_loc)).argmin()  # Find the closest grid point in x
        y_idx = (np.abs(y - canister.y_loc)).argmin()  # Find the closest grid point in y
        
        # Plot the canister's grid position on the heatmap
        ax.plot(x_idx, y_idx, 'o', color='green', markersize=8, label="Gas Canister Center")
    
    # -------------------------------------------------------------------------
    
    # Plot optimal point if it exists
    for point in detector_locs:
        ax.plot(point[1], point[0], '*', color='yellow', markersize=15, label="Optimal Location")

    # -------------------------------------------------------------------------
    
    # After plotting everything, remove duplicate legend entries
    handles, labels = ax.get_legend_handles_labels()
    unique_labels = dict(zip(labels, handles))  # remove duplicate labels
    plt.legend(unique_labels.values(), unique_labels.keys())
    
    # Show plot
    plt.tight_layout()
    plt.show()
    
    # -------------------------------------------------------------------------
        
    no_detectors = len(detector_locs)
    print(f"{no_detectors} are needed to cover all the canisters.")
    
    for point in detector_locs:
        
        # Reverse the point's coordinates for grid intuition
        point = point[::-1]
    
        # Calculate the physical coordinates for the given indices
        length_det = point[0] * parameters.dx
        width_det = point[1] * parameters.dy
            
        print(f"Optimal detector placement is at grid point of {point}",
              f"This is at a length of {length_det:.4g} and width of {width_det:.4g}",
              sep="\n")

# -------------------------------------------------------------------------

gas_canisters = [ 
    sg.GasCan2D(x_loc=0.5, y_loc=0.5,radius=0.05, concentration=1.0),
    sg.GasCan2D(x_loc=0.1, y_loc=0.2,radius=0.05, concentration=1.0),
    sg.GasCan2D(x_loc=0.9, y_loc=0.4,radius=0.05, concentration=1.0),
    sg.GasCan2D(x_loc=0.1, y_loc=0.9,radius=0.05, concentration=1.0)
]

scrubbers = [sg.Scrubber2D(x_loc=0.3, y_loc=0.5, radius=0.1, efficiency=0),
             sg.Scrubber2D(x_loc=0.7, y_loc=0.5, radius=0.1, efficiency=0)]   

max_distance(parameters=parameters, gas_canisters=gas_canisters, build_matrix=build_matrix)

# For the meshgrid, since using (indexing='ij'), the rows correspond to Y and
# columns correspond to X.



       
