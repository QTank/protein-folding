import sympy
import numpy as np
from math import *
import String_Tool
import EnergyTable
from VQE_Qiskit import VQE
from qiskit.quantum_info import SparsePauliOp

class Protein:
    def __init__(self, amino_acid_chain):
        self.N = len(amino_acid_chain)
        self.amino_acid_chain = amino_acid_chain
        self.lamda_overlap = 311  # sympy.Symbol('l_overlap')
        self.qubits_used = 2 * self.N - 6
        self.qubits = [sympy.Symbol("q%d" % i) for i in range(self.qubits_used + 1)]
        self.z_list = [sympy.Symbol("z%d" % i) for i in range(self.qubits_used + 1)]
        #self.energy_table, self.index_dict = EnergyTable.read_data()
        #self._hamiltonian = self.get_hamiltonian()

    def axis_lattice(self, index):
        if index == 0:
            return [0, 1, 0, 0]


        if index == 1:
            return [1, 0, 0, 0]

        start_index = 2*index - 4
        q0 = self.qubits[start_index]
        q1 = self.qubits[start_index+1]

        return [(1-q0)*(1-q1), (1-q0)*q1, q0*(1-q1), q0*q1]

    def distance(self, i, j):
        d = [0, 0, 0, 0]

        for index in range(i, j):

            if j % 2 == 0:
                d += self.axis_lattice(index)
            else:
                d -= self.axis_lattice(index)
        return d


protein = Protein("ASD")
print(protein.distance(0, 2))

