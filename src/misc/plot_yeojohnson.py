#! /usr/bin/env python2.7
'''
A small program to plot the Yeo-Johnson transformation [1] and the corresponding
probability density function.

[1] I.-K. Yeo and R. A. Johnson, "A new family of power transformations to improve
normality or symmetry," Biometrika, vol. 87, no. 4, pp. 954-959, Dec. 2000.

@author:     Alexander Schmidt-Richberg
@copyright:  2014 Imperial College London. All rights reserved.
@contact:    a.schmidt-richberg@imperial.ac.uk
'''
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx


def main():
    lambdas = [0.0, 0.5, 1.0, 1.5, 2.0]
    mu = 1
    sigma = 1

    plt.subplot(1, 2, 1)
    # plt.xlim(-2, 6)
    plt.ylim(-10, 30)
    plt.xlabel('Response $y$')
    plt.ylabel('Yeo-Johnson transformation $\psi(y,\lambda)$')
    plt.axhline(0, color='k')
    plt.axvline(0, color='k')
    Y = np.arange(-5, 10, 0.01)
    plot(Y, yeojohnson, lambdas, mu, sigma)
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.xlim(-2, 6)
    # plt.ylim(-0.01, 1)
    plt.xlabel('Response $y$')
    plt.ylabel('Probability $g_Y(y)$')
    plt.axhline(0, color='k')
    plt.axvline(0, color='k')
    Y = np.arange(-5, 10, 0.01)
    plot(Y, yeojohnson_density, lambdas, mu, sigma)
    plt.legend()

    plt.tight_layout()
    plt.show()


def yeojohnson(y, lmbda, mu, sigma):
    ymu = float(y) / mu
    if ymu < 0:
        if lmbda == 2:
            return -math.log(-ymu + 1) / sigma
        else:
            return -(math.pow(-ymu + 1, 2 - lmbda) - 1) / ((2 - lmbda) * sigma)
    else:
        if lmbda == 0:
            return math.log(ymu + 1) / sigma
        else:
            return (math.pow(ymu + 1, lmbda) - 1) / (lmbda * sigma)


def yeojohnson_inv(phi, lmbda, mu, sigma):
    if phi < 0:
        if lmbda == 2:
            return mu * (1 - math.exp(-phi * sigma))
        else:
            return mu * (1 - math.pow(-(2 - lmbda) * sigma * phi + 1, 1 / ((2 - lmbda) * sigma)))
    else:
        if lmbda == 0:
            return mu * (math.exp(phi * sigma) - 1)
        else:
            return mu * (math.pow(phi * sigma * lmbda + 1, 1 / (lmbda * sigma)) - 1)


def std_normal_dist(x):
    return np.exp(-0.5 * x * x) / np.sqrt(2 * np.pi)


def yeojohnson_density(y, lmbda, mu, sigma):
    return (1 / sigma) * \
        std_normal_dist((yeojohnson(y, lmbda, 1, 1) - mu) / sigma) * \
        math.pow(math.fabs(y) + 1, math.copysign(1, y) * (lmbda - 1))


def boxcox(y, lmbda, mu, sigma):
    ymu = float(y) / mu
    if lmbda == 0:
        return math.log(ymu) / sigma
    else:
        return (math.pow(ymu, lmbda) - 1) / (lmbda * sigma)


def plot(X, transformation, lambdas, mu, sigma):
    yjn_cmap = cmx.ScalarMappable(
        norm=colors.Normalize(vmin=np.min(lambdas), vmax=np.max(lambdas)),
        cmap=plt.get_cmap('cool'))

    for lmbda in lambdas:
        Y = [transformation(x, lmbda, mu, sigma) for x in X]
        plt.plot(X, Y, label='$\lambda=' + str(lmbda) + '$', color=yjn_cmap.to_rgba(lmbda))

if __name__ == "__main__":
    main()
