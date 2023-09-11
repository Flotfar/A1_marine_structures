
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


###################################
#### Scour profile plot
###################################

def Scour_profile_plot(D, W1, W2, S_eq, condition):
    # Coordinates of three points (W1, W2, and S)
    p1 = (-W1, 0)
    p2 = (W2, 0)
    p3 = (0, -S_eq)

    # Calculate the axis limits based on W1 and W2
    x_min = -5
    x_max = 4
    y_min = -3
    y_max = 2

    # Create a fixed-size figure and axis
    fig, ax = plt.subplots(figsize=(6, 6))  # Adjust the size as needed

    # Create a circle patch
    circle = Circle((0, D / 2), radius=D / 2, fill=False, color='black')  # (0, D/2) is the center of the circle

    # Add the circle patch to the axis
    ax.add_patch(circle)

    # Plot the three input points
    ax.plot(*p1, marker='o', markersize=3, color='black', label='W1')
    ax.plot(*p2, marker='o', markersize=3, color='black', label='W2')
    ax.plot(*p3, marker='o', markersize=3, color='black', label='S')

    # Set the aspect ratio to 'equal' to ensure the circle looks circular
    ax.set_aspect('equal', adjustable='box')

    # Set axis limits to include the circle
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # Add a grid
    ax.set_axisbelow(True)
    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7, zorder=0)

    # Plot a solid line for y=0
    plt.axhline(0, color='black', linestyle='dashed', linewidth=1)

    # Optional: Add labels or customize the plot further
    ax.set_xlabel('Lenght [m]')
    ax.set_ylabel('Height [m]')
    ax.set_title('Scour profile in ' + condition)

    # Display the plot
    plt.show()


