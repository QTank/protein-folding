from ProteinFolding import ProteinFold


def get_figure(amino, sequence):
    from matplotlib import pyplot as plt
    p = [0, 0]
    fold = sequence
    length = 6
    x = [0]
    y = [0]

    for i in range(length - 1):
        if fold[2 * i] == fold[2 * i + 1]:
            p[1] += 2 * int(fold[2 * i + 1]) - 1
        else:
            p[0] += 1 - 2 * int(fold[2 * i])

        x.append(p[0])
        y.append(p[1])

    plt.plot(x, y, '-')
    plt.scatter(x, y, s=200, zorder=3)
    plt.axis('off')
    plt.margins(0.1)  # enough margin so that the large scatter dots don't touch the borders
    plt.gca().set_aspect('equal')
    for i in range(length):
        plt.text(x[i], y[i] + 0.1, amino[i], fontsize=15)
    plt.show()


if __name__ == "__main__":
    amino_sequence = "PSVKMA"
    proteinObj = ProteinFold(amino_sequence)
    sequence = proteinObj.fold()
    get_figure(amino_sequence, sequence)
