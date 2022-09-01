#!/usr/bin/env python
import matplotlib.pyplot as plt 
import numpy as np 
import sys 

'''

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

n = 100

# For each set of style and range settings, plot n random points in the box
# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
for c, m, zlow, zhigh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:
    xs = randrange(n, 23, 32)
    ys = randrange(n, 0, 100)
    zs = randrange(n, zlow, zhigh)
    ax.scatter(xs, ys, zs, c=c, marker=m)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
'''


def visualize3d(data):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for d in data:
        ax.scatter(d[0], d[1], d[2],c='blue',alpha=.05)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_zlabel('PC3')
    plt.show()

def visualize(data):
    fig, axs = plt.subplots(nrows=2,sharex=True,dpi=800)
    fig.suptitle('PCA chr21.GRCh38.20181129.phased')
    for i,ax in enumerate(axs):
        ax.set_box_aspect(1)
        if i==0:
            ax.scatter(data[:,0], data[:,1],alpha=.1,s=1)
            ax.set_ylabel('PC2')
        if i==1:
            ax.scatter(data[:,0], data[:,2],alpha=.1, s=1)
            ax.set_xlabel('PC1')
            ax.set_ylabel('PC3')
    plt.savefig('pca.png')


if __name__ == "__main__":
    eigvecs_file = sys.argv[1] 
    eigvec = np.genfromtxt(eigvecs_file)
    eigvec = eigvec[:,2:]
    visualize(eigvec)