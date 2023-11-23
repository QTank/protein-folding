import csv
import numpy as np


def read_data():
    energy = np.zeros((20, 20))
    index_dict = {}

    with open("EnergyData.csv", "r") as f:
        csv_reader = csv.reader(f, delimiter=",")
        line = 0
        for each in csv_reader:
            for i in range(1, len(each)):
                energy[line][i-1] = float(each[i])

            index_dict[each[0]] = line
            line += 1

        for i in range(20):
            for j in range(i, 20):
                energy[i][j] = min(energy[i][j], energy[j][i])
                energy[j][i] = min(energy[i][j], energy[j][i])
    return energy, index_dict

