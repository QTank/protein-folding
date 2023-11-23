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


def rand(length):
    import random
    chain = []
    table = ["A", "K", "M", "S", "V", "P"]
    for i in range(length):
        chain.append(table[random.randint(0, 5)])
    return "".join(chain)


def rand_only(length):
    import random
    chain = []
    table = ["A", "K", "M", "S", "V", "P"]
    for i in range(length):
        a = random.choice(table)
        while a in chain:
            a = random.choice(table)
        chain.append(a)

    return "".join(chain)



