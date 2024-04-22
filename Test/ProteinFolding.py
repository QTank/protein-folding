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
        self.lamda_overlap = 311  # sympy.Symbol('l_overlap')
        self.qubits_used = 2 * self.N - 5
        self.qubits = [sympy.Symbol("q%d" % i) for i in range(self.qubits_used + 1)]
        self.z_list = [sympy.Symbol("z%d" % i) for i in range(self.qubits_used + 1)]
        self.energy_table, self.index_dict = EnergyTable.read_data()
        self._hamiltonian = self.get_hamiltonian()

    # find the x location of nth amino acid
    def get_x_position(self, n):
        if n == 1 or n == 2:
            return n - 1

        res = 1 + self.qubits[1]
        for i in range(3, n):
            res += self.qubits[2 * i - 3] - self.qubits[2 * i - 4]

        return res

    # find the y location of nth amino acid
    def get_y_position(self, n):
        if n == 1 or n == 2:
            return 0

        res = self.qubits[1] - 1
        for i in range(3, n):
            res += self.qubits[2 * i - 4] + self.qubits[2 * i - 3] - 1

        return res

    # encode the right direction in x-axis on the step
    def direction_x_right(self, step):
        if step == 1:
            return 1

        if step == 2:
            return self.qubits[1]

        return (1 - self.qubits[2 * step - 4]) * self.qubits[2 * step - 3]

    # encode the left direction in x-axis on the step
    def direction_x_left(self, step):
        if step == 1 or step == 2:
            return 1

        return self.qubits[2 * step - 4] * (1 - self.qubits[2 * step - 3])

    # encode the up direction in y-axis on the step
    def direction_y_up(self, step):
        if step == 1 or step == 2:
            return 1

        return self.qubits[2 * step - 4] * self.qubits[2 * step - 3]

    # encode the down direction in y-axis on the step
    def direction_y_down(self, step):
        if step == 1:
            return 1

        if step == 2:
            return 1 - self.qubits[1]

        return (1 - self.qubits[2 * step - 4]) * (1 - self.qubits[2 * step - 3])

    # find the square of distance between ith and jth amino acid
    def distance(self, i, j):
        x_i = self.get_x_position(i)
        x_j = self.get_x_position(j)
        y_i = self.get_y_position(i)
        y_j = self.get_y_position(j)
        return (x_i - x_j) ** 2 + (y_i - y_j) ** 2

    def A_item(self, i, j):
        n = ceil(log2(abs(j - i) ** 2))

        A = self.qubits[self.qubits_used + 1]
        for i in range(1, n):
            A += self.qubits[self.qubits_used + 1 + i] * (2 ** i)
        self.qubits_used += n
        return A

    # encode
    def overlap_item(self, i, j):
        n = ceil(log2(abs(j - i) ** 2))

        overlap = (2 ** n - self.distance(i, j) - self.A_item(i, j)) ** 2
        return overlap

    def overlap_val(self):
        res = 0
        for i in range(1, self.N - 3):
            for j in range(i + 4, self.N + 1):
                if abs(1 - i + j) % 2 == 1:
                    res += self.overlap_cir(i, 3)  # self.overlap_item(i, j)
        return res * self.lamda_overlap

    def back_item(self, i):  # i step
        return self.direction_x_right(i) * self.direction_x_left(i + 1) + self.direction_x_left(
            i) * self.direction_x_right(i + 1) + self.direction_y_up(i) * self.direction_y_down(
            i + 1) + self.direction_y_down(i) * self.direction_y_up(i + 1)

    def overlap_cir(self, index, length):
        lap = self.direction_x_right(index) * self.direction_y_up(index + 1) * self.direction_x_left(
            index + 2) * self.direction_y_down(index + 3)
        lap += self.direction_x_right(index) * self.direction_y_down(index + 1) * self.direction_x_left(
            index + 2) * self.direction_y_up(index + 3)
        lap += self.direction_x_left(index) * self.direction_y_up(index + 1) * self.direction_x_right(
            index + 2) * self.direction_y_down(index + 3)
        lap += self.direction_x_left(index) * self.direction_y_down(index + 1) * self.direction_x_right(
            index + 2) * self.direction_y_up(index + 3)
        lap += self.direction_y_up(index) * self.direction_x_right(index + 1) * self.direction_y_down(
            index + 2) * self.direction_x_left(index + 3)
        lap += self.direction_y_up(index) * self.direction_x_left(index + 1) * self.direction_y_down(
            index + 2) * self.direction_x_right(index + 3)
        lap += self.direction_y_down(index) * self.direction_x_right(index + 1) * self.direction_y_up(
            index + 2) * self.direction_x_left(index + 3)
        lap += self.direction_y_down(index) * self.direction_x_left(index + 1) * self.direction_y_up(
            index + 2) * self.direction_x_right(index + 3)
        return lap

    def back_val(self):
        res = 0
        if self.N < 3:
            return res

        energy_back = self.qubits[1] * self.qubits[2] + self.qubits[2] * self.qubits[3] - 2 * self.qubits[1] * \
                      self.qubits[2] * self.qubits[3]
        for i in range(3, self.N - 1):
            energy_back += self.back_item(i)

        return energy_back * self.lamda_overlap

    def pair_energy(self):
        energy = 0
        if self.N < 4:
            return energy

        for i in range(0, self.N - 3):
            for j in range(i + 3, self.N):
                amino_index_i = self.index_dict[self.amino_acid_chain[i]]
                amino_index_j = self.index_dict[self.amino_acid_chain[j]]
                val = self.energy_table[amino_index_i][amino_index_j] * 10
                #energy += val * (2 - self.distance(i + 1, j + 1))
                d = self.distance(i+1, j+1)
                dd = (10-d)*(5-d)*(2-d)*0.11*0.25 + (10-d)*(5-d)*(d-1)*0.125*0.33 * 0.25 \
                    + (10-d)*(d-2)*(d-1)*0.5*0.33*0.25 * 0.125 * 0.125 + (d-5)*(d-2)*(d-1)*0.5*0.128*0.11 * 0.1 **2

                energy += val * (dd*d)
        return energy

    def hamiltonian_function(self):
        #h = self.back_val() + self.overlap_val() + self.pair_energy()
        h = self.back_val() + self.overlap_val() + self.pair_energy_update()
        return self.reduce(h.expand(), self.qubits_used)

    def get_hamiltonian(self):
        hamiltonian = self.hamiltonian_function()
        return self.to_ising(hamiltonian).expand()

    def reduce(self, h, power):
        for step in range(power, 1, -1):
            for i in range(self.qubits_used + 1):
                h = h.subs(self.qubits[i] ** (step), self.qubits[i])
        return h

    def get_index_arr(self, hamiltonian):
        arr = hamiltonian.split(" ")
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

    # convert the hamiltonian to ising model
    def to_ising(self, hamiltonian):
        for w in range(self.qubits_used + 1):
            hamiltonian = hamiltonian.subs(self.qubits[w], 0.5 * (1 - self.z_list[w]))
        return hamiltonian

    # parse the operator and coefficient on the hamiltonian of the protein chain
    def get_index_coefficient(self):
        h = self._hamiltonian
        index, coefficient = self.get_index_arr(str(h))
        return index, coefficient

    # iteration the whole binary string in order to verify the result of VQE
    def test_energy(self):
        energy_dict = {}
        energy = self._hamiltonian
        for i in range(2 ** self.qubits_used):
            tp = energy
            encode = String_Tool.int_to_binary(i, self.qubits_used)
            for j in range(self.qubits_used):
                v = 1 if int(encode[j]) == 0 else -1
                tp = tp.subs(self.z_list[j + 1], v)

            energy_dict[encode] = tp

        return dict(sorted(energy_dict.items(), key=lambda item: item[1]))

    # compute hamiltonian on the binary string
    def binary_to_hval(self, binary_string):
        h = self._hamiltonian
        for i in range(len(binary_string)):
            val = 1 if int(binary_string[i]) == 0 else -1
            h = h.subs(self.z_list[i+1], val)

        return h

    def convert_op(self, op):
        for i in range(len(op)):
            op[i] = op[i].replace("1", "Z")
            op[i] = op[i].replace("0", "I")
        return op

    def fold(self, ansatz, optimizer, initial_point):
        op, co = self.get_index_coefficient()
        op = self.convert_op(op)
        co = co[1:] # delete the constant

        ob = SparsePauliOp(op, np.array(co))
        vqe = VQE(ansatz, ob, optimizer, initial_point)
        theta = vqe.solve().x
        res = vqe.get_expection(theta)
        prob = res.quasi_dists[0]
        top = max(prob.items(), key=lambda x: x[1])
        opt_sequence = String_Tool.int_to_binary(top[0], self.qubits_used)
        return "010%s" % opt_sequence

    def select_opt_sequence(self, top):
        min_val = 10**7
        opt_sequence = ""
        for each in top:
            sequence = String_Tool.int_to_binary(each[0], self.qubits_used)
            val_hamiltonian = self.binary_to_hval(sequence)
            if val_hamiltonian < min_val:
                min_val = val_hamiltonian
                opt_sequence = sequence
        return "010%s" % opt_sequence

    # compute the val of distance between ith and jth amino acid
    def get_distance(self, i, j, binary_string):
        d = self.distance(i, j)
        for i in range(len(binary_string)):
            d = d.subs(self.qubits[i+1], binary_string[i])
        return d

    # update the method of computing interaction energy between the pair of amino acids, but it does not work now
    def pair_energy_update(self):
        energy = 0
        if self.N < 4:
            return energy

        for i in range(0, self.N - 3):
            for j in range(i + 3, self.N, 2):
                # print("interaction ", self.interaction[i][j])
                amino_index_i = self.index_dict[self.amino_acid_chain[i]]
                amino_index_j = self.index_dict[self.amino_acid_chain[j]]
                val = self.energy_table[amino_index_i][amino_index_j] * 10 * 2
                if j - i + 1 == 4:
                    energy += val * self.tag_4(i + 1)
                elif j - i + 1 == 6:

                    energy += val * self.tag_6(i + 1)

        return energy

    # manually encode the overlap in a square with the length of 4, but it will not be good.
    def tag_4(self, index):
        if index + 3 > self.N:
            return 0

        lap = self.direction_x_right(index) * self.direction_y_up(index + 1) * self.direction_x_left(index + 2)
        lap += self.direction_x_right(index) * self.direction_y_down(index + 1) * self.direction_x_left(index + 2)
        lap += self.direction_x_left(index) * self.direction_y_up(index + 1) * self.direction_x_right(index + 2)
        lap += self.direction_x_left(index) * self.direction_y_down(index + 1) * self.direction_x_right(index + 2)
        lap += self.direction_y_up(index) * self.direction_x_right(index + 1) * self.direction_y_down(index + 2)
        lap += self.direction_y_up(index) * self.direction_x_left(index + 1) * self.direction_y_down(index + 2)
        lap += self.direction_y_down(index) * self.direction_x_right(index + 1) * self.direction_y_up(index + 2)
        lap += self.direction_y_down(index) * self.direction_x_left(index + 1) * self.direction_y_up(index + 2)
        return lap

    # manually encode the overlap in a rectangle th the length of 6, but it will not be good.
    # it will be merged to one method with tag_4
    def tag_6(self, index):
        if index + 5 > self.N:
            return 0

        lap = self.direction_x_right(index) * self.direction_y_down(index + 1) * self.direction_y_down(
            index + 2) * self.direction_x_left(index + 3) * self.direction_y_up(index + 4)
        lap += self.direction_x_right(index) * self.direction_x_right(index + 1) * self.direction_y_down(
            index + 2) * self.direction_x_left(index + 3) * self.direction_x_left(index + 4)
        lap += self.direction_x_right(index) * self.direction_y_up(index + 1) * self.direction_y_up(
            index + 2) * self.direction_x_left(index + 3) * self.direction_y_down(index + 4)
        lap += self.direction_x_right(index) * self.direction_x_right(index + 1) * self.direction_y_up(
            index + 2) * self.direction_x_left(index + 3) * self.direction_x_left(index + 4)
        lap += self.direction_x_right(index) * self.direction_y_down(index + 1) * self.direction_x_left(
            index + 2) * self.direction_x_left(index + 3) * self.direction_y_up(index + 4)
        lap += self.direction_x_right(index) * self.direction_y_up(index + 1) * self.direction_x_left(
            index + 2) * self.direction_x_left(index + 3) * self.direction_y_down(index + 4)

        lap += self.direction_x_left(index) * self.direction_y_down(index + 1) * self.direction_y_down(
            index + 2) * self.direction_x_right(index + 3) * self.direction_y_up(index + 4)
        lap += self.direction_x_left(index) * self.direction_x_left(index + 1) * self.direction_y_down(
            index + 2) * self.direction_x_right(index + 3) * self.direction_x_right(index + 4)
        lap += self.direction_x_left(index) * self.direction_y_up(index + 1) * self.direction_y_up(
            index + 2) * self.direction_x_right(index + 3) * self.direction_y_down(index + 4)
        lap += self.direction_x_left(index) * self.direction_x_left(index + 1) * self.direction_y_up(
            index + 2) * self.direction_x_right(index + 3) * self.direction_x_right(index + 4)
        lap += self.direction_x_left(index) * self.direction_y_down(index + 1) * self.direction_x_right(
            index + 2) * self.direction_x_right(index + 3) * self.direction_y_up(index + 4)
        lap += self.direction_x_left(index) * self.direction_y_up(index + 1) * self.direction_x_right(
            index + 2) * self.direction_x_right(index + 3) * self.direction_y_down(index + 4)

        lap += self.direction_y_down(index) * self.direction_x_right(index + 1) * self.direction_x_right(
            index + 2) * self.direction_y_up(index + 3) * self.direction_x_left(index + 4)
        lap += self.direction_y_down(index) * self.direction_y_down(index + 1) * self.direction_x_right(
            index + 2) * self.direction_y_up(index + 3) * self.direction_y_up(index + 4)
        lap += self.direction_y_down(index) * self.direction_x_left(index + 1) * self.direction_x_left(
            index + 2) * self.direction_y_up(index + 3) * self.direction_x_right(index + 4)
        lap += self.direction_y_down(index) * self.direction_y_down(index + 1) * self.direction_x_left(
            index + 2) * self.direction_y_up(index + 3) * self.direction_y_up(index + 4)
        lap += self.direction_y_down(index) * self.direction_x_right(index + 1) * self.direction_y_up(
            index + 2) * self.direction_y_up(index + 3) * self.direction_x_left(index + 4)
        lap += self.direction_y_down(index) * self.direction_x_left(index + 1) * self.direction_y_up(
            index + 2) * self.direction_y_up(index + 3) * self.direction_x_right(index + 4)

        lap += self.direction_y_up(index) * self.direction_x_right(index + 1) * self.direction_x_right(
            index + 2) * self.direction_y_down(index + 3) * self.direction_x_left(index + 4)
        lap += self.direction_y_up(index) * self.direction_y_up(index + 1) * self.direction_x_right(
            index + 2) * self.direction_y_down(index + 3) * self.direction_y_down(index + 4)
        lap += self.direction_y_up(index) * self.direction_x_left(index + 1) * self.direction_x_left(
            index + 2) * self.direction_y_down(index + 3) * self.direction_x_right(index + 4)
        lap += self.direction_y_up(index) * self.direction_y_up(index + 1) * self.direction_x_left(
            index + 2) * self.direction_y_down(index + 3) * self.direction_y_down(index + 4)
        lap += self.direction_y_up(index) * self.direction_x_right(index + 1) * self.direction_y_down(
            index + 2) * self.direction_y_down(index + 3) * self.direction_x_left(index + 4)
        lap += self.direction_y_up(index) * self.direction_x_left(index + 1) * self.direction_y_down(
            index + 2) * self.direction_y_down(index + 3) * self.direction_x_right(index + 4)

        return lap


if __name__ == '__main__':
    test = ProteinFold("YYPDETET")
    print(test.tag_8(1))

