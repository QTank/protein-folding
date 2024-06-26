# Quantum protein-folding on a planar lattice model with VQE
The method is from the paper https://www.nature.com/articles/srep00571, which provides the technique for mapping protein construction onto qubits and folding proteins on a planar lattice model on quantum annealing. For a detailed explanation of the approach, you can refer to the paper. My code constructed the Hamiltonian according to the paper. However, instead of using quantum annealing, I implemented the method with SapmpleVQE on Qiskit. This method required fewer qubits compared to quantum annealing. 


You can run the solution.py and get the structure of the protein sequence YTDPETGT. Of course, you also can get the structure of the protein that you want to fold. 
