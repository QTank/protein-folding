import sympy
import numpy as np
from math import *
import String_Tool
import EnergyTable
from VQE_Qiskit import VQE
from qiskit.quantum_info import SparsePauliOp


class ProteinFold:
    def __init__(self, amino_acid_chain):
        self.N = len(amino_acid_chain)
        self.amino_acid_chain = amino_acid_chain
        self.overlap_lambda = 311  # sympy.Symbol('l_overlap')
        self.redundant_direction = 100
        self.qubits_used = 3 * self.N - 3 * 2
        self.ancilla_len = 0
        self.qubits = [sympy.Symbol("q%d" % i) for i in range(self.qubits_used+10*2)]
        self.ancilla = [sympy.Symbol("q%d" % i) for i in range(self.qubits_used, self.qubits_used+self.N)]
        self.z_list = [sympy.Symbol("z%d" % i) for i in range(len(self.qubits) + len(self.ancilla))]
        self.energy_table, self.index_dict = EnergyTable.read_data()
        self.distance_all = [1, 2, 5, 11, 13]#[1, 2, 4, 5, 8, 9]#, 10, 13, 16, 17, 18, 20, 25, 26, 29, 36, 37, 49]
        self._hamiltonian = self.get_hamiltonian()

        # direction: left:101, right:110, up:001, down:010, in: 100, out:111, redutant: 000, 011


    def left(self, q0, q1, q2):
        return q0 * (1-q1) * q2

    def right(self, q0, q1, q2):
        return q0 * q1 * (1-q2)

    def up(self, q0, q1, q2):
        return (1 - q0) * (1 - q1) * q2

    def down(self, q0, q1, q2):
        return (1-q0) * q1 * (1-q2)

    def in_direction(self, q0, q1, q2):
        return q0 * (1-q1) * (1-q2)

    def out_direction(self, q0, q1, q2):
        return q0 * q1 * q2

    def back_each(self, index):
        q_current0, q_current1, q_current2 = self.qubits[3*index], self.qubits[3*index+1], self.qubits[3*index+2]
        q_latter0, q_latter1, q_latter2 = self.qubits[3*index+3], self.qubits[3*index+1+3], self.qubits[3*index+2+3]
        return self.right(q_current0, q_current1, q_current2) * self.left(q_latter0, q_latter1, q_latter2) \
            + self.right(q_current0, q_current1, q_current2) * self.left(q_latter0, q_latter1, q_latter2) \
            + self.up(q_current0, q_current1, q_current2) * self.down(q_latter0, q_latter1, q_latter2) \
            + self.down(q_current0, q_current1, q_current2) * self.up(q_latter0, q_latter1, q_latter2) \
            + self.in_direction(q_current0, q_current1, q_current2) * self.out_direction(q_latter0, q_latter1, q_latter2) \
            + self.out_direction(q_current0, q_current1, q_current2) * self.in_direction(q_latter0, q_latter1, q_latter2)

    def back(self):
        val = 0
        for i in range(self.N-3):
            val += self.back_each(i)
        return val

    def get_location(self, index):
        if index == 0:
            return 0, 0, 0

        if index == 1:
            return 1, 0, 0

        qubit_start = 3 * (index - 2)
        pre_x, pre_y, pre_z = self.get_location(index - 1)
        x = pre_x + self.qubits[qubit_start] * self.qubits[qubit_start + 1] * (1 - self.qubits[qubit_start + 2]) \
            - self.qubits[qubit_start] * (1 - self.qubits[qubit_start + 1]) * self.qubits[qubit_start + 2]

        y = pre_y + self.qubits[qubit_start] * (1 - self.qubits[qubit_start + 1]) * (
                1 - self.qubits[qubit_start + 2]) \
            - self.qubits[qubit_start] * self.qubits[qubit_start + 1] * self.qubits[qubit_start + 2]

        z = pre_z + (1 - self.qubits[qubit_start]) * (1 - self.qubits[qubit_start + 1]) * self.qubits[
            qubit_start + 2] \
            - (1 - self.qubits[qubit_start]) * self.qubits[qubit_start + 1] * (1 - self.qubits[qubit_start + 2])

        return x, y, z

    def redundant(self):
        val = 0
        for i in range(0, self.qubits_used, 3):
            val += (1 - self.qubits[i]) * (1 - self.qubits[i + 1]) * (1 - self.qubits[i + 2]) \
                   + (1 - self.qubits[i]) * self.qubits[i + 1] * self.qubits[i + 2]

        return self.redundant_direction * val

    def get_distance(self, i, j):
        i_loc = self.get_location(i)
        j_loc = self.get_location(j)
        a = (i_loc[0] - j_loc[0]) ** 2 + (i_loc[1] - j_loc[1]) ** 2 + (i_loc[2] - j_loc[2]) ** 2
        return self.reduce(a.expand())#a.expand()(i_loc[0] - j_loc[0]) ** 2 + (i_loc[1] - j_loc[1]) ** 2 + (i_loc[2] - j_loc[2]) ** 2

    def overlap(self):
        val_all = 0
        for i in range(0, self.N - 4):
            for j in range(i + 4, self.N, 2):
                val_all += self.overlap_item(i, j)
                print("test", i, j)
        return self.overlap_lambda * val_all.as_poly()

    def overlap_item(self, i, j):
        n = ceil(log2(abs(i-j)**2))
        return (2**n - self.get_distance(i, j) - self.ancilla_item(i, j)) ** 2

    def ancilla_item(self, i, j):
        n = ceil(log2(abs(i-j)**2))
        val = 0
        for i in range(0, n):
            val += self.qubits[self.qubits_used + self.ancilla_len+i] * 2 ** i

        self.ancilla_len += n
        return val

    def interaction_energy(self):
        energy = 0
        for i in range(0, self.N - 1):
            for j in range(i + 1, self.N):
                print("i, j", i, j)
                distance = self.get_distance(i, j)
                val, real = 1, 1
                for d_index in range(1, len(self.distance_all)):
                    val *= (self.distance_all[d_index] - distance) ** 2
                    real *= 1 / (self.distance_all[d_index] - 1) ** 2
                val = 0
                energy += val * real
        return energy

    def build_hamiltonian(self):
        h = self.overlap()# + self.interaction_energy() + self.redundant_direction
        return h#h.expand()

    def reduce(self, h):
        for step in range(6, 1, -1):
            for i in range(self.qubits_used):
                h = h.subs(self.qubits[i] ** (step), self.qubits[i])
        return h

    def get_hamiltonian(self):
        h = self.build_hamiltonian()
        #h = 1
        return self.reduce(h)

    def to_z(self, h):
        for i in range(self.qubits_used):
            h = h.subs(self.qubits[i], 0.5 * (1 - self.z_list[i]))
        return h

    def get_coefficient(self, h):
        arr = h.split(" ")
        index_coefficient = {}
        for i in range(0, len(arr), 2):
            item = arr[i].split("*")
            coefficient = float(item[0])
            if i > 0 and arr[i - 1] == '-':
                coefficient *= -1.0

            tmp = ["0" for w in range(self.qubits_used)]
            tmp_string = "".join(tmp)
            if 'z' not in arr[i] and 'q' not in arr[i]:
                index_coefficient["".join(tmp)] = coefficient
            else:
                for w in range(1, len(item)):
                    tmp[int(item[w][1:]) - 1] = "1"
                index_coefficient["".join(tmp)] = coefficient

        index = []
        coefficient = [index_coefficient[tmp_string]]

        for k in index_coefficient:
            if k != tmp_string:
                index.append(k)
                coefficient.append(index_coefficient[k])
        return index, coefficient


test = ProteinFold("PSVAVVA")
print(test._hamiltonian)
#print(test.get_distance(0, 3))
