import String_Tool
from ProteinFolding import ProteinFold
from matplotlib import pyplot as plt
from qiskit.algorithms.optimizers import SLSQP
from qiskit.circuit.library import EfficientSU2
import numpy as np


def plot_grid_with_arrows_and_text(marked_locations):
    # Mark the locations on the grid
    for idx, location in enumerate(marked_locations):
        x, y, direction = location

        # Add arrow to indicate direction
        if idx != len(marked_locations) - 1:
            plt.arrow(x, y, *get_arrow_direction(direction, True),
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
    #plt.xticks(range(1 + 1))
    #plt.yticks(range(1 + 1))
    for i in range(0,5):
        plt.scatter(range(0,5), [i]*5, c='b')
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


def plot_structure(amino_sequence, encode_sequence):
    p = [0, 0, ""]
    location = []
    for i in range(len(encode_sequence)//2):

        direction = get_arrow_direction(encode_sequence[2*i:2*i+2])
        p[2] = encode_sequence[2*i:2*i+2]
        location.append(p.copy())
        p[0] = p[0] + direction[0]
        p[1] = p[1] + direction[1]

    location.append(p.copy())

    plot_grid_with_arrows_and_text(location)


def test(amino_sequence):
    proteinObj = ProteinFold(amino_sequence)
    num_qubits = proteinObj.qubits_used
    ansatz = EfficientSU2(num_qubits, reps=2, entanglement='circular', insert_barriers=True)
    # ansatz = RealAmplitudes(3, entanglement='circular', reps=2, insert_barriers=True)
    np.random.seed(100)  # seed for reproducibility
    initial_point = np.random.random(ansatz.num_parameters)
    slsqp = SLSQP()
    sequence = proteinObj.fold(ansatz, slsqp, initial_point)
    print("The sequence: ", sequence)
    interact_pair = get_interact_pair(proteinObj, sequence[3:], len(amino_sequence))

    return sequence, interact_pair
    #energy = proteinObj.test_energy()
    #for aa, v in energy.items():
    #    print(f"int: {int(aa, 2)}, sequence: {aa}, val: {v}")


def get_interact_pair(proteinObj, sequence, amino_len):
    interact = []
    for i in range(1, amino_len-2):
        for j in range(i+3, amino_len+1):
            d = proteinObj.get_distance(i, j, sequence)
            if d == 1:
                interact.append((i, j))
    return interact



if __name__ == "__main__":
    protein_sequence = "DAYAQW"
    #encode_sequence, interact_pair = test(protein_sequence)
    encode_sequence = "0100001011"
    plot_structure(protein_sequence, encode_sequence)
    #get_figure(protein_sequence, encode_sequence, [])