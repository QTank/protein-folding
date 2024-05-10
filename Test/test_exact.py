import numpy as np
import EnergyTable

class ExactProtein:

    def __init__(self, protein_sequence, bitstring):
        self.protein_sequence = protein_sequence
        self.protein_len = len(protein_sequence)
        self.qubit_len = 2 * (len(protein_sequence) - 1)
        self.bitstring = bitstring
        self.turn_list = self.set_turn_list()

        self.energy_table, self.index_dict = EnergyTable.read_energy_data("../EnergyData.csv")

    def set_turn_list(self):
        turn_list = []
        for i in range(self.protein_len):

            turn_list.append(self.bitstring[2*i:2*i+2])
        return turn_list


    def get_position(self, index):
        return self.get_x_position(index), self.get_y_position(index)

    def get_x_position(self, index):
        if index == 0:
            return 0

        x = 0
        for i in range(index):
            turn = self.turn_list[i]
            if turn == "01":
                x += 1
            if turn == "10":
                x -= 1

        return x

    def get_y_position(self, index):
        if index == 0:
            return 0
        y = 0

        for i in range(index):
            turn = self.turn_list[i]
            if turn == "00":
                y -= 1

            if turn == "11":
                y += 1

        return y

    def get_distance(self, i, j):
        x_i, y_i = self.get_position(i)
        x_j, y_j = self.get_position(j)

        return (x_i - x_j) ** 2 + (y_i - y_j) **2

    def judge(self):
        for i in range(self.protein_len):
            for j in range(i+2, self.protein_len):
                if (j - i) % 2 == 0:
                    d = self.get_distance(i, j)
                    if d == 0:
                        return False
        return True

    def interaction(self):
        h = 0

        for i in range(self.protein_len):
            for j in range(i+2, self.protein_len):
                d = self.get_distance(i, j)

                if d == 1:
                    amino_i_index = self.index_dict[self.protein_sequence[i]]
                    amino_j_index = self.index_dict[self.protein_sequence[j]]
                    energy = self.energy_table[amino_i_index][amino_j_index]
                    print(f"{i}, {j}, {self.protein_sequence[i]}, {self.protein_sequence[j]}, {energy}")
                    h += energy
        return h
