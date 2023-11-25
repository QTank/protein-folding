import String_Tool
from ProteinFolding import ProteinFold
from matplotlib import pyplot as plt
from qiskit.algorithms.optimizers import SLSQP
from qiskit.circuit.library import EfficientSU2
import numpy as np


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

    plt.axis('off')
    plt.margins(0.1)  # enough margin so that the large scatter dots don't touch the borders
    plt.gca().set_aspect('equal')
    for i in range(len(amino)):
        plt.text(x[i], y[i] + 0.1, amino[i], fontsize=15)
    plt.show()


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
    protein_sequence = "CYVIRV"
    protein_sequence = "DAYAQW"
    encode_sequence, interact_pair = test(protein_sequence)
    get_figure(protein_sequence, encode_sequence, interact_pair)
