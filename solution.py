import time, VQE, util
from protein_folding import ProteinFolding
import plot_protein_folding


def solution(protein_sequence):
    start_time = time.time()
    folding = ProteinFolding(protein_sequence)
    qubit_op = folding.get_qubit_op()
    print("start to execute VQE ... ")
    opt_bitstring, opt_value = VQE.get_min(qubit_op)
    executing_time = round(time.time() - start_time, 2)
    print(f"bitstring: {opt_bitstring}, value:{opt_value}, executing time: {executing_time} s")
    return opt_bitstring


def get_interaction_pair(sequence):
    interact = []
    point = util.calculate_position(sequence)
    amino_len = len(point)

    for i in range(0, amino_len):
        for j in range(i + 3, amino_len):
            i_loc = point[i]
            j_loc = point[j]
            distance = (i_loc[0]- j_loc[0]) ** 2 + (i_loc[1] - j_loc[1]) **2
            if distance == 1:
                interact.append((i, j))
    return interact


if __name__ == '__main__':
    protein_sequence = "YYDPETGT" # 00 01 11 01 11 10 10     00 01 11 11 10 10 00 01 00
    encoding_string = solution(protein_sequence)
    interaction = get_interaction_pair("00011101111010")

    plot_protein_folding.plot_structure(protein_sequence, encoding_string, interaction)




