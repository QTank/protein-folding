import matplotlib.pyplot as plt

def plot_grid_with_arrows_and_text(grid_size, marked_locations):

    # Mark the locations on the grid
    for idx, location in enumerate(marked_locations):
        x, y, direction = location

        # Add arrow to indicate direction
        plt.arrow(x, y, *get_arrow_direction(direction),
                  head_width=0.1, head_length=0.1, fc='red', ec='red')

        # Add text with number
        plt.text(x+0.1, y+0.1, str(idx + 1), color='black',
                 ha='center', va='center', fontsize=13)

    # Mark the locations with red dots
    for location in marked_locations:
        plt.plot(location[0], location[1], 'ro')  # +0.5 to center the marker

    # Customize the plot
    plt.title('Protein folding')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.xticks(range(grid_size + 1))
    plt.yticks(range(grid_size + 1))
    for i in range(0,5):
        plt.scatter(range(0,5), [i]*5, c='b')
    # Show the plot
    plt.show()

def get_arrow_direction(direction):
    # Define arrow direction based on encoding
    d = 0.15
    if direction == '01':
        return (1-d, 0)
    elif direction == '10':
        return (-1+d, 0)
    elif direction == '11':
        return (0, 1-d)
    elif direction == '00':
        return (0, -1+d)

# Example usage
grid_size = 5
marked_locations = [(0, 0, '01'), (1, 0, '00'), (1, -1, '00')]
plot_grid_with_arrows_and_text(grid_size, marked_locations)