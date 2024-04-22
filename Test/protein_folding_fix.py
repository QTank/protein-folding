import time
import util, EnergyTable
from bead import Bead


class ProteinFolding:

    def __init__(self, protein_sequence):
        self.protein_sequence = protein_sequence
        self.protein_len = len(protein_sequence)
        self.used_qubit_len = 2*self.protein_len - 2
        self.lamda_overlap = 311000
        self.lamda_energy = 100*2
        self.energy_table, self.index_dict = EnergyTable.read_data()
        self.bead_list = self.set_bead_list()
        self.identity = util.build_full_identity(self.used_qubit_len)

    def set_bead_list(self):
        return [Bead(self.protein_sequence[index], index, self.used_qubit_len) for index in range(self.protein_len)]

    def get_x_position(self, index):
        if index < 2:
            return index
        x = 1

        for loc in range(2, index):
            bead = self.bead_list[index]
            turn = bead.get_turn()
            x += turn[1] - turn[2]

        return x

    def get_y_position(self, index):
        if index < 2:
            return 0
        y = 0

        for loc in range(2, index):
            bead = self.bead_list[index]
            turn = bead.get_turn()
            y += turn[3] - turn[0]

        return y

    def get_distance(self, i, j):
        x_i = self.get_x_position(i)
        y_i = self.get_y_position(i)

        x_j = self.get_x_position(j)
        y_j = self.get_y_position(j)

        return ((x_i - x_j) @ (x_i - x_j) + (y_i - y_j) @ (y_i - y_j)).reduce()

    def back_item(self, i):

        bead_cur = self.bead_list[i]
        bead_next = self.bead_list[i+1]

        cur_turn = bead_cur.get_turn()
        next_turn = bead_next.get_turn()

        h = cur_turn[0] @ next_turn[3]
        h += cur_turn[1] @ next_turn[2]
        h += cur_turn[2] @ next_turn[1]
        h += cur_turn[3] @ next_turn[0]


        return h.reduce()

    def get_h_back(self):
        h_back = 0
        for i in range(2, self.protein_len-1):
            h_back += self.back_item(i)

        return h_back.reduce() * self.lamda_overlap

    def overlap_cir(self, index):
        h_overlap = 0
        turn_clockwise = [3, 1, 0, 2]
        turn_counterclockwise = [1, 3, 2, 0]

        for i in range(len(turn_clockwise)):
            h_clockwise, h_counterclockwise = self.identity, self.identity

            for bead_loc in range(4):
                bead_turn = self.bead_list[index + bead_loc].get_turn()[turn_clockwise[(i + bead_loc) % 4]]

                h_clockwise @= bead_turn

                bead_turn = self.bead_list[index + bead_loc].get_turn()[turn_counterclockwise[(i+bead_loc) % 4]]

                h_counterclockwise @= bead_turn

            h_overlap += h_clockwise + h_counterclockwise

        return h_overlap.reduce()

    def get_h_overlap(self):
        h_overlap = 0
        for i in range(self.protein_len - 3):
            h_overlap += self.overlap_cir(i)

        for i in range(self.protein_len - 5):

            h_overlap += self.overlap_cir_6(i)
        return h_overlap * self.lamda_overlap

    def overlap_cir_6(self, index):
        h = 0
        turn_counterclockwise = [1, 1, 3, 2, 2, 0]
        turn_clockwise = [3, 1, 1, 0, 2, 2]
        turn_counterclockwise_2 = [1, 3, 3, 2, 0, 0]
        turn_clockwise_2 = [3, 3, 1, 0, 0, 2]
        for i in range(len(turn_clockwise)):

            h_clockwise = self.identity
            h_counterclockwise = self.identity
            h_counterclockwise_2 = self.identity
            h_clockwise_2 = self.identity
            for bead_loc in range(6):
                bead_turn = self.bead_list[index + bead_loc].get_turn()[turn_counterclockwise[(i+bead_loc) % 4]]
                h_clockwise @= bead_turn

                bead_turn = self.bead_list[index + bead_loc].get_turn()[turn_counterclockwise[(i+bead_loc) % 6]]
                h_counterclockwise @= bead_turn

                bead_turn = self.bead_list[index + bead_loc].get_turn()[turn_clockwise_2[(i+bead_loc) % 6]]
                h_clockwise_2 @= bead_turn

                bead_turn = self.bead_list[index + bead_loc].get_turn()[turn_counterclockwise_2[(i+bead_loc) % 6]]

                h_counterclockwise_2 @= bead_turn

            h += h_clockwise + h_counterclockwise + h_clockwise_2 + h_counterclockwise_2

        return h.reduce()

    def get_pair_interaction(self):
        h_energy = 0
        for i in range(self.protein_len - 3):
            for j in range(i+3, self.protein_len):
                amino_index_i = self.index_dict[self.protein_sequence[i]]
                amino_index_j = self.index_dict[self.protein_sequence[j]]
                val = self.energy_table[amino_index_i][amino_index_j]
                if j - i + 1 == 4:
                    h_energy += self.tag_4(i) * val
                if j - i + 1 == 6:
                    h_energy += self.tag_6(i) * val

        return h_energy.reduce() * self.lamda_energy

    def tag_4(self, index):
        if index + 3 > self.protein_len:
            return 0

        turn_counterclockwise = [1, 3, 2, 0]
        turn_clockwise = [3, 1, 0, 2]
        h = 0

        for i in range(len(turn_clockwise)):

            h_counterclockwise = self.identity
            h_clockwise = self.identity
            for bead_loc in range(3):
                bead_turn = self.bead_list[index + bead_loc].get_turn()[turn_counterclockwise[(i+bead_loc) % 4]]
                h_counterclockwise @= bead_turn

                bead_turn = self.bead_list[index + bead_loc].get_turn()[turn_clockwise[(i+bead_loc) % 4]]
                h_clockwise @= bead_turn

            h += h_counterclockwise + h_counterclockwise
        return h.reduce()

    def tag_6(self, index):
        if index + 5 > self.protein_len:
            return 0
        h = 0
        turn_counterclockwise = [1, 1, 3, 2, 2, 0]
        turn_clockwise = [3, 1, 1, 0, 2, 2]
        turn_counterclockwise_2 = [1, 3, 3, 2, 0, 0]
        turn_clockwise_2 = [3, 3, 1, 0, 0, 2]
        for i in range(len(turn_clockwise)):
            h_clockwise = self.identity
            h_counterclockwise = self.identity
            h_counterclockwise_2 = self.identity
            h_clockwise_2 = self.identity
            for bead_loc in range(5):
                bead_turn = self.bead_list[index + bead_loc].get_turn()[turn_clockwise[(i + bead_loc) % 6]]
                h_clockwise @= bead_turn

                bead_turn = self.bead_list[index + bead_loc].get_turn()[turn_counterclockwise[(i + bead_loc) % 6]]
                h_counterclockwise @= bead_turn

                bead_turn = self.bead_list[index + bead_loc].get_turn()[turn_clockwise_2[(i + bead_loc) % 6]]
                h_clockwise_2 @= bead_turn

                bead_turn = self.bead_list[index + bead_loc].get_turn()[turn_counterclockwise_2[(i + bead_loc) % 6]]
                h_counterclockwise_2 @= bead_turn

            h += h_clockwise + h_counterclockwise + h_clockwise_2 + h_counterclockwise_2

        return h.reduce()

    def get_qubit_op(self):
        print(f"construct the energy Hamiltonian")
        start_time = time.time()
        h = self.get_h_back() + self.get_h_overlap() + self.get_pair_interaction()
        print(f"the constructing time: {time.time() - start_time}")
        return h.reduce()

    # def get_qubit_fix_op(self):
    #     print(f"construct the energy Hamiltonian with fixing qubits")
    #     start_time = time.time()
    #     h = self.get_h_back() + self.get_h_overlap() + self.get_pair_interaction()
    #     print(f"the constructing time: {time.time() - start_time}")
    #     h = h.reduce()
    #     h = qubit_fixing._fix_qubits(h.reduce())
    #     h, unused_qubits = qubit_number_reducer.remove_unused_qubits(h)
    #     return h.reduce()