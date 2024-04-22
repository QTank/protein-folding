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

def is_valid_move(x, y, n, path):
    return -n <= x <= n and -n <= y <= n and (x, y) not in path and (x, y) not in [(0,0)]

def find_paths(x, y, dest_x, dest_y, length, path=[], paths=[]):
    if length == 0 and x == dest_x and y == dest_y:
        paths.append(path.copy())
        return
    n = 2

    for nx, ny in [(x+1, y), (x-1, y), (x, y+1), (x, y-1)]:
        if is_valid_move(nx, ny, n, path):
            path.append((nx, ny))
            find_paths(nx, ny, dest_x, dest_y, length - 1, path, paths)
            path.pop()


def add_direction(location):
    start = (0,0)
    loc_with_direction = []

    for idx, loc in enumerate(location):

        if idx != 0:
            start = location[idx-1]
        direction = ""
        x = loc[0] - start[0]
        y = loc[1] - start[1]
        if x == 1 and y == 0:
            direction = "01"
        elif x == -1 and y == 0:
            direction = "10"
        elif x == 0 and y == 1:
            direction = "11"
        else:
            direction = "00"

        loc_with_direction.append((start[0], start[1], direction))

    return loc_with_direction