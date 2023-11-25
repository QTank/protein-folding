from qiskit.primitives import Estimator
from qiskit.primitives import Sampler


class VQE:
    def __init__(self, ansatz, ob, optimizer, init_point):
        self.estimator = Estimator()
        self.ob = ob  # SparsePauliOp(ob[0], ob[1])
        self.ansatz = ansatz
        self.optimizer = optimizer
        self.init_point = init_point
        self.sampler = Sampler()

    def solve(self):
        return self.optimizer.minimize(self.get_energy, self.init_point)

    def get_expection(self, theta):
        sampler = Sampler()
        self.ansatz.measure_all()
        task = sampler.run(self.ansatz, theta)
        return task.result()

    def get_energy(self, theta):
        return self.estimator.run(self.ansatz, self.ob, theta).result().values[0]
