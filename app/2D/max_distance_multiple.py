# =============================================================================   
# Importing
# =============================================================================   

import numpy as np
import model_params as model
import scrubbers_canisters as sg
import matplotlib.pyplot as plt

from workspace import run_simulation
from matplotlib.animation import FuncAnimation
from itertools import combinations

# Modules for plotting gas hulls
# Note - need to install alphashape
from alphashape import alphashape
from shapely.geometry import Point  

parameters = model.parameters()

# =============================================================================   
# Notes
# =============================================================================   

# - Throughout this code matrices have been made to resemble the gridspace
#   as much as possible, mainly for better intuition of what's going on. However
#   the "y-direction" is sometimes flipped using np.flipud() to perform indexing
#   operations in order to conform with Python's default indexing system.

# - It is assumed that the origin for our gridspace is a particular corner of
#   the room, and this corner remains the origin throughout. Thus when we refer
#   to "X" length, where X is a number, we are meaning the length from this 
#   origin point, and likewise for the width. 

# =============================================================================   
# Function
# =============================================================================   
 
def max_distance(parameters: model.parameters, gas_canisters: list, scrubbers: list, run_simulation):
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
    # Finding points locations in grid that satisfy 1% of initial conc.
    # =========================================================================  

    # Initialise storage lists for plotting
    satisfied_points = []
    gas_hulls = []
    canister_hulls = []
    canister_centres = []
     
    # Running the model for each canister
    for canister in gas_canisters:
  
        U = run_simulation(gas_canister = canister, scrubbers = scrubbers)     
        
        # Unpack the flattened U matrix into a 3D matrix
        # ( with dimensions Nx_points, Ny_points, Nt_points )
        U_3D = np.zeros((parameters.Nt_points, parameters.Nx_points, parameters.Ny_points))

        for i in range(parameters.Nt_points):
            U_3D[i] = np.flipud(U[i].reshape((parameters.Nx_points, parameters.Ny_points)))
        
        # Define gridspace (discretized spatial points on x and y grid)
        x = np.linspace(0, parameters.len_x, parameters.Nx_points)
        y = np.linspace(0, parameters.len_y, parameters.Ny_points)        
        
        # Get canister location
        x_loc = canister.x_loc
        y_loc = canister.y_loc
         
        # Find indices of canister location in the spatial grid
        # (closest values in gridspace to x_loc and y_loc)
        x_idx = (np.abs(x - x_loc)).argmin()  
        y_idx = (np.abs(y - y_loc)).argmin()  
        
        # Store the center of container in grid format
        canister_centres.append(np.array([x_idx, y_idx]))
        
        # Get threshold value for detector 
        threshold = 0.01 * canister.concentration  # 1% threshold
        
        # Find points where concentration is at least as much as the threshold
        # on the spatial grid, across all timesteps
        threshold_mask = np.any(U_3D >= threshold, axis=0) 
        
        # Extract the grid points for these grid points that satisfy condition
        points = np.argwhere(np.flipud(threshold_mask))  
        
        # Reverse the coordinates to keep in grid format rather than row x column
        points = points[:, ::-1]
            
        # Store Boolean of points which satisfy this condition for the canister
        satisfied_points.append(threshold_mask)
        
        # ---------------------------------------------------------------------
        # Plot for checking 
        # ---------------------------------------------------------------------  
        
        # Plot to check if the grid points are correct
        
        # fig, ax = plt.subplots(figsize=(6, 6))
        # ax.set_xlabel('X')
        # ax.set_ylabel('Y')
        # ax.set_xlim(0, parameters.Nx_points)
        # ax.set_ylim(0, parameters.Ny_points)
        # ax.plot(points[:, 0], points[:, 1], 'o', color='cyan', markersize=5)
        # plt.show()
        
        # Commented out by default as only used for checking
            
        # =========================================================================  
        # Finding convex/concave hull for gas/canister/scrubber for plotting
        # =========================================================================    
            
        # Compute a hull for gas using alpha shapes to visualise in plot
        alpha = 0.1  # Lower alpha = tighter hull
        if len(points) >= 4:  # Alpha shapes require at least 4 points
            gas_hull = alphashape(points, alpha)
            if not gas_hull.is_empty:
                # Store the computed hull in storage list
                gas_hulls.append(gas_hull)
             
        # Compute a hull for canister using alpha shapes
        canister_mask = canister.get_initial_concentration(x, y)
        canister_points = np.argwhere(np.flipud(canister_mask))  
        canister_points = canister_points[:, ::-1]
        canister_shape = alphashape(canister_points, alpha)
        canister_hulls.append(canister_shape)    
    
    # Exit canister loop and compute hulls for scrubbers
    scrubber_hulls = []
    scrubber_centres = []
    
    for scrubber in scrubbers:
                 
        # Find and store centre of scrubbers, similar to before
        x_loc_scrub = scrubber.x_loc
        y_loc_scrub = scrubber.y_loc
        x_idx_scrub = (np.abs(x - x_loc_scrub)).argmin()  
        y_idx_scrub = (np.abs(y - y_loc_scrub)).argmin()  
        scrubber_centres.append(np.array([x_idx_scrub, y_idx_scrub]))
        
        # Find and store hull of scrubbers, similar to that of canister
        scrubber_mask = scrubber.get_affected_indices(x,y)
        scrubber_points = np.argwhere(np.flipud(scrubber_mask))  
        scrubber_points = scrubber_points[:, ::-1]
        scrubber_shape = alphashape(scrubber_points, alpha)
        scrubber_hulls.append(scrubber_shape)    
                     
    # =========================================================================  
    # Finding optimal location for detector(s)
    # ========================================================================= 
     
    # Generate the index combinations of all our canisters in descending order
    num_canisters = len(gas_canisters)
    total_combs = []
    
    for r in range(num_canisters, 0, -1):
        total_combs.extend(combinations(range(num_canisters), r))
      
    # Iterate over each combination layer and add successful combinations to storage list
    satisfied_combs = []  # Storage for index combinations that have overlap
    
    # Don't iterate over last layer yet
    for i in range(num_canisters, 1, -1):
        
        # Check if any of the combinations at this level can be satisfied
        current_combs = [comb for comb in total_combs if len(comb) == i]
        
        # First check if any indices have already been covered:
        # Flatten the tuples to find all unique indices that are already in satisfied_combs
        flattened_satisfied_combs = set(index for comb in satisfied_combs for index in comb)
        
        for comb in current_combs:
    
            if all(idx in flattened_satisfied_combs for idx in comb):
                # Skip combinations that only include already covered canisters
                continue
            
            combined_mask = np.logical_and.reduce([satisfied_points[idx] for idx in comb])
            
            if np.any(combined_mask):
                # If combination is satisfied, add to list
                satisfied_combs.append(comb)

        
    # Flatten the tuples to find all unique indices that are included in satisfied_combs
    flattened_satisfied_combs = set(index for comb in satisfied_combs for index in comb)
    
    # Find missing elements by checking which indices are not in the satisfied combinations
    missing_cans = [can for can in range(num_canisters) if can not in flattened_satisfied_combs]
    
    # Add each missing element as a single-element tuple to satisfied_combs
    satisfied_combs.extend((elem,) for elem in missing_cans)

    # Now that we have found number of detectors we need and for what gas containers, then:     
    detector_locs = []
    
    # For each of our satisfied combinations
    for comb in satisfied_combs:
        
        # Create boolean mask for this index combination
        combined_mask = np.logical_and.reduce([satisfied_points[idx] for idx in comb])

        # Extract the "True" grid points
        candidate_points = np.argwhere(np.flipud(combined_mask))
                
        # Reverse the coordinates to keep in grid format rather than row x column
        candidate_points = candidate_points[:, ::-1]
    
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
    
    U = run_simulation(gas_canister = canister, scrubbers = scrubbers)     
    
    # Unpack the flattened U matrix into a 3D matrix
    # ( with dimensions Nx_points, Ny_points, Nt_points )
    U_3D = np.zeros((parameters.Nt_points, parameters.Nx_points, parameters.Ny_points))

    # Don't flip as just using for plotting
    for k in range(parameters.Nt_points):
        U_3D[k] = U[k].reshape((parameters.Nx_points, parameters.Ny_points))
        
    # -------------------------------------------------------------------------
    
    # Initialising plot
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_aspect('equal')
    ax.set_xlim(0, parameters.Nx_points)
    ax.set_ylim(0, parameters.Ny_points)
    ax.set_xlabel("X Grid Index")
    ax.set_ylabel("Y Grid Index")
    ax.set_title("1% Gas Spread Radius with Container and Optimal Detector Locations")

    # -------------------------------------------------------------------------
        
    # Plot the heatmap
    min_concentration = 0
    max_concentration = 0.1
    
    heatmap = ax.imshow(U_3D[0], cmap='magma', vmin=min_concentration, vmax=max_concentration, interpolation='nearest', origin='lower')

    # Animation function
    def update(frame):
        heatmap.set_array(U_3D[frame])
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
            ax.plot(hull_x, hull_y, color='red', linestyle='-', linewidth=2, label="Max Distance to detect 1%")
            
    # -------------------------------------------------------------------------
            
    # Iterate through each canister hull and plot it
    for canister_shape in canister_hulls:
        if not canister_shape.is_empty:
            can_x, can_y = canister_shape.exterior.xy  
            ax.plot(can_x, can_y, color='blue', linestyle='-', linewidth=2, label="Gas Canister Circumference")
            
    # -------------------------------------------------------------------------
    
    # Plot the locations of the gas canisters on the grid
    for centre in canister_centres:
        ax.plot(centre[0], centre[1], 'o', color='green', markersize=8, label="Gas Canister Center")

    # -------------------------------------------------------------------------
    
    # Plot optimal point if it exists
    for point in detector_locs:
        ax.plot(point[0], point[1], '*', color='yellow', markersize=15, label="Optimal Location")

    # -------------------------------------------------------------------------
    
    # Plot scrubbers
    for scrubber_centre in scrubber_centres:
        ax.plot(scrubber_centre[0], scrubber_centre[1], 'o', color='cyan', markersize=8, label="Scrubber Center")
        
    for scrubber_shape in scrubber_hulls:
        if not scrubber_shape.is_empty:
            scrub_x, scrub_y = scrubber_shape.exterior.xy  
            ax.plot(scrub_x, scrub_y, color='cyan', linestyle='-', linewidth=2, label="Scrubber Circumference")
            
    # -------------------------------------------------------------------------
    
    # After plotting everything, remove duplicate legend entries
    handles, labels = ax.get_legend_handles_labels()
    unique_labels = dict(zip(labels, handles))  # remove duplicate labels
    plt.legend(unique_labels.values(), unique_labels.keys())
    
    # Show plot
    plt.tight_layout()
    plt.show()
    
    # -------------------------------------------------------------------------
    
    # Print statements
    
    # Number of detectors needed
    no_detectors = len(detector_locs)
    print(f"\n{no_detectors} detectors are needed to cover all of the gas canisters.\n")
    
    # Print header of the table
    print(f"{'Detector':<10} {'Grid Point (x, y)':<20} {'Length (m)':<15} {'Width (m)':<15}")
    print("=" * 60)
    
    # Iterate through detector locations and print their placement details
    for idx, point in enumerate(detector_locs, start=1):
        # Calculate the physical coordinates for the given indices
        length_det = point[0] * parameters.dx
        width_det = point[1] * parameters.dy
        
        # Print row in tabular format
        print(f"{idx:<10} {str(tuple(point)):<20} {length_det:<15.4g} {width_det:<15.4g}")
                
    return U_3D

# =============================================================================   
# Running the algorithm
# =============================================================================  

gas_canisters = [
    sg.GasCan2D(x_loc=3.2, y_loc=3.2, radius=0.23, concentration=1.0),
    sg.GasCan2D(x_loc=0.4, y_loc=3.6, radius=0.23, concentration=1.0),
    sg.GasCan2D(x_loc=1.6, y_loc=2.4, radius=0.23, concentration=1.0),
    sg.GasCan2D(x_loc=2.8, y_loc=0.8, radius=0.23, concentration=1.0),
    sg.GasCan2D(x_loc=3.6, y_loc=1.6, radius=0.23, concentration=1.0)
]

# Scrubbers turned off in this run
scrubbers = [sg.Scrubber2D(x_loc=1, y_loc=2, radius=0.46, efficiency=0),
                sg.Scrubber2D(x_loc=3, y_loc=2, radius=0.46, efficiency=0)]

U_3D = max_distance(parameters=parameters, 
                    gas_canisters=gas_canisters, 
                    scrubbers=scrubbers, 
                    run_simulation=run_simulation)
