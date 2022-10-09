import matplotlib.pyplot as plt
import numpy as np
import sys

def plot2DPCA(df):
    plt.figure(figsize=[10,10])
    plt.scatter(df[:,0],df[:,1],alpha=.5)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA on cell genotypes: {}'.format(file))
    plt.show()


if __name__ == "__main__":
    file = sys.argv[1]
    df = np.loadtxt(file)
    df = df[:,2:]
    plot2DPCA(df)
