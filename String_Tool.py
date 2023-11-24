def int_to_binary(index, length):
    binary = bin(index)[2:]
    bin_len = len(binary)

    if length > bin_len:
        pre = "".join(["0" for i in range(length - bin_len)])
        return "%s%s" % (pre, binary)

    return binary


def convert_op(op):
    for i in range(len(op)):
        op[i] = op[i].replace("1", "Z")
        op[i] = op[i].replace("0", "I")
    return op


def rand_sequence(length):
    import random
    chain = []
    table = ['C','M','F','I','L','V','W','Y','A','G','T','S','N','Q','D','E','H','R','K','P']
    for i in range(length):
        chain.append(table[random.randint(0, len(table)-1)])
    return "".join([table[random.randint(0, len(table)-1)] for i in range(length)])


def rand_only(length):
    import random
    chain = []
    table = ['C','M','F','I','L','V','W','Y','A','G','T','S','N','Q','D','E','H','R','K','P']
    for i in range(length):
        a = random.choice(table)
        while a in chain:
            a = random.choice(table)
        chain.append(a)

    return "".join(chain)
