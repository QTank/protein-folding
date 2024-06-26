from matplotlib import pyplot as plt


def plot_grid_with_arrows_and_text(amino_sequence, marked_locations, range_distance, interact):
    # Mark the locations on the grid
    for idx, location in enumerate(marked_locations):
        x, y, direction = location

        # Add arrow to indicate direction
        if idx != len(marked_locations) - 1:
            plt.arrow(x, y, *get_arrow_direction(direction, True),
                  head_width=0.1, head_length=0.1, fc='red', ec='red')

        # Add text with number
        plt.text(x+0.1, y+0.1, str(amino_sequence[idx]), color='black',
                 ha='center', va='center', fontsize=13)

    # Mark the locations with red dots
    # Customize the plot
    plt.title('Protein folding')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')

    for each in interact:
        interact_x = [marked_locations[each[0]][0], marked_locations[each[1]][0]]
        interact_y = [marked_locations[each[0]][1], marked_locations[each[1]][1]]
        plt.plot(interact_x, interact_y, '--', c='y')

    for i in range(range_distance[2], range_distance[3]):
        x_range = range(range_distance[0], range_distance[1])
        plt.scatter(x_range, [i]*len(x_range), c='b')

    plt.legend()
    # Show the plot
    plt.show()


def get_arrow_direction(direction, b=False):
    # Define arrow direction based on encoding
    d = 0.15 if b else 0
    if direction == '01':
        return (1-d, 0)
    elif direction == '10':
        return (-1+d, 0)
    elif direction == '11':
        return (0, 1-d)
    elif direction == '00':
        return (0, -1+d)


def get_figure(amino, fold, interact):
    p = [0, 0]
    x, y = [0], [0]

    for i in range(len(amino) - 1):
        if fold[2 * i] == fold[2 * i + 1]:
            p[1] += 2 * int(fold[2 * i + 1]) - 1
        else:
            p[0] += 1 - 2 * int(fold[2 * i])

        x.append(p[0])
        y.append(p[1])
    plt.plot(x, y, '-')
    plt.scatter(x, y, s=200, zorder=3)

    for each in interact:
        interact_x = [x[each[0]-1], x[each[1]-1]]
        interact_y = [y[each[0]-1], y[each[1]-1]]
        plt.plot(interact_x, interact_y, '--', 'r')
    plt.grid()
    plt.axis('off')
    plt.margins(0.1)  # enough margin so that the large scatter dots don't touch the borders
    plt.gca().set_aspect('equal')
    for i in range(len(amino)):
        plt.text(x[i], y[i] + 0.1, amino[i], fontsize=15)
    plt.show()


def plot_structure(amino_sequence, encode_sequence, interact):
    p = [0, 0, ""]
    location = []

    for i in range(len(encode_sequence)//2):

        direction = get_arrow_direction(encode_sequence[2*i:2*i+2])
        p[2] = encode_sequence[2*i:2*i+2]
        location.append(p.copy())
        p[0] = p[0] + direction[0]
        p[1] = p[1] + direction[1]

    location.append(p.copy())

    x_range = [row[0] for row in location]
    second_range = [row[1] for row in location]

    range_distance = [min(x_range)-1, max(x_range)+2, min(second_range)-1, max(second_range)+2]
    plot_grid_with_arrows_and_text(amino_sequence, location, range_distance, interact)