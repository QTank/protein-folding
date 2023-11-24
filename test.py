from ProteinFolding import ProteinFold
from matplotlib import pyplot as plt
from qiskit.algorithms.optimizers import SLSQP
from qiskit.algorithms.optimizers import SPSA
from qiskit.utils import algorithm_globals
from qiskit.circuit.library import EfficientSU2
import numpy as np

def get_figure(amino, fold):
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
    plt.axis('off')
    plt.margins(0.1)  # enough margin so that the large scatter dots don't touch the borders
    plt.gca().set_aspect('equal')
    for i in range(len(amino)):
        plt.text(x[i], y[i] + 0.1, amino[i], fontsize=15)
    plt.show()


if __name__ == "__main__":
    amino_sequence = "PSVKMA"
    proteinObj = ProteinFold(amino_sequence)
    num_qubits = proteinObj.qubits_used
    ansatz = EfficientSU2(num_qubits, reps=2, entanglement='circular', insert_barriers=True)
    # ansatz = RealAmplitudes(3, entanglement='circular', reps=2, insert_barriers=True)
    np.random.seed(100)  # seed for reproducibility
    initial_point = np.random.random(ansatz.num_parameters)
    slsqp = SLSQP()
    sequence = proteinObj.fold(ansatz, slsqp, initial_point)
    print("The sequence: ", sequence)
    get_figure(amino_sequence, sequence)

    energy = proteinObj.test_energy()
    for aa, v in energy.items():
        print(f"int: {int(aa, 2)}, sequence: {aa}, val: {v}")