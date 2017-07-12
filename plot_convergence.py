import matplotlib.pyplot as plt
import numpy as np
import os

def plot(f):
    if not os.path.exists(f): return None

    data = np.loadtxt(f)
    dofs, data = data[:, 0], data[:, 1:]

    plt.figure()
    if data.shape[1] == 1:
        plt.semilogx(dofs, data, '-o', basex=10.)
    else:
        with open(f, 'r') as foo:
            header = foo.readline().strip()
            lmbdas = header.split(' ')[1:]
            print lmbdas
        
        data = data.T
        for lmbda, niters in zip(lmbdas, data):
            plt.semilogx(dofs, niters, '-o', basex=10, label=r'$\lambda=%1.0E$' % float(lmbda))
        plt.legend(loc='best')

    plt.xlabel('dofs')
    plt.ylabel('iterations')
    plt.title(os.path.splitext(f)[0][1:])

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    for f in ['.lagrange_primal.txt', '.lagrange_mixed.txt', '.energy.txt']:
        plot(f)
    plt.show()
