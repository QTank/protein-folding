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


def add_direction(location):
    start = (0,0)
    loc_with_direction = []

    for idx, loc in enumerate(location):

        if idx != 0:
            start = location[idx-1]
        direction = ""
        x = loc[0] - start[0]
        y = loc[1] - start[1]
        if x == 1 and y == 0:
            direction = "01"
        elif x == -1 and y == 0:
            direction = "10"
        elif x == 0 and y == 1:
            direction = "11"
        else:
            direction = "00"

        loc_with_direction.append((start[0], start[1], direction))

    return loc_with_direction

# Example usage
marked_locations = [(0, 0, '01'), (1, 0, '00'), (1, -1, '00')]
#loc = [[(-1, 0), (-2, 0), (-2, 1), (-1, 1), (0, 1), (1, 1), (1, 0)], [(-1, 0), (-2, 0), (-2, -1), (-1, -1), (0, -1), (1, -1), (1, 0)], [(-1, 0), (-1, 1), (0, 1), (1, 1), (2, 1), (2, 0), (1, 0)], [(-1, 0), (-1, 1), (0, 1), (0, 2), (1, 2), (1, 1), (1, 0)], [(-1, 0), (-1, 1), (-1, 2), (0, 2), (1, 2), (1, 1), (1, 0)], [(-1, 0), (-1, 1), (-1, 2), (0, 2), (0, 1), (1, 1), (1, 0)], [(-1, 0), (-1, -1), (0, -1), (1, -1), (2, -1), (2, 0), (1, 0)], [(-1, 0), (-1, -1), (0, -1), (0, -2), (1, -2), (1, -1), (1, 0)], [(-1, 0), (-1, -1), (-1, -2), (0, -2), (1, -2), (1, -1), (1, 0)], [(-1, 0), (-1, -1), (-1, -2), (0, -2), (0, -1), (1, -1), (1, 0)], [(0, 1), (1, 1), (2, 1), (2, 0), (2, -1), (1, -1), (1, 0)], [(0, 1), (1, 1), (1, 2), (2, 2), (2, 1), (2, 0), (1, 0)], [(0, 1), (-1, 1), (-1, 2), (0, 2), (1, 2), (1, 1), (1, 0)], [(0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1), (1, 0)], [(0, 1), (0, 2), (1, 2), (2, 2), (2, 1), (1, 1), (1, 0)], [(0, 1), (0, 2), (1, 2), (2, 2), (2, 1), (2, 0), (1, 0)], [(0, 1), (0, 2), (1, 2), (1, 1), (2, 1), (2, 0), (1, 0)], [(0, -1), (1, -1), (2, -1), (2, 0), (2, 1), (1, 1), (1, 0)], [(0, -1), (1, -1), (1, -2), (2, -2), (2, -1), (2, 0), (1, 0)], [(0, -1), (-1, -1), (-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0)], [(0, -1), (-1, -1), (-1, -2), (0, -2), (1, -2), (1, -1), (1, 0)], [(0, -1), (0, -2), (1, -2), (2, -2), (2, -1), (1, -1), (1, 0)], [(0, -1), (0, -2), (1, -2), (2, -2), (2, -1), (2, 0), (1, 0)], [(0, -1), (0, -2), (1, -2), (1, -1), (2, -1), (2, 0), (1, 0)]]for idx, location in enumerate(loc):plot_grid_with_arrows_and_text(grid_size, add_direction(location))
a = []

for i in range(0, 7+1):
    for j in range(0, 7+1):
        if i + j <= 7:
            a.append(i**2 + j**2)
a.sort()
print(list(dict.fromkeys(a)))
