import util


class Bead:

    def __init__(self, amino, index, used_qubit_len):
        self.amino = amino
        self.used_qubit_len = used_qubit_len
        self.identity = util.build_full_identity(used_qubit_len)
        self.index = index
        self.indicator_0, self.indicator_1 = self.set_indicator()
        self._turn_0 = self.set_turn_0()
        self._turn_1 = self.set_turn_1()
        self._turn_2 = self.set_turn_2()
        self._turn_3 = self.set_turn_3()

    def set_indicator(self):
        qubit_index = 2 * self.index
        return util.build_indicator_qubit(self.used_qubit_len, qubit_index), \
            util.build_indicator_qubit(self.used_qubit_len, qubit_index + 1)

    def set_turn_0(self):
        return (self.identity - self.indicator_0) @ (self.identity - self.indicator_1)

    def set_turn_1(self):
        return (self.identity - self.indicator_0) @ self.indicator_1

    def set_turn_2(self):
        return self.indicator_0 @ (self.identity - self.indicator_1)

    def set_turn_3(self):
        return self.indicator_0 @ self.indicator_1

    def get_turn(self):
        return [self._turn_0, self._turn_1, self._turn_2, self._turn_3]