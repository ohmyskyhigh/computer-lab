import pylab as p
import numpy as np


def choose(n, i):
    if 0 <= i <= n:
        p = 1
        for t in range(0, min(i, n - i), 1):
            p = (p * (n - t)) // (t + 1)
        return p
    else:
        return 0


def bezier1(a):
    n = len(a) - 1
    t = np.linspace(0, 1, 101)
    c = np.zeros([101, 2])
    for i in range(0, 101):
        b = []
        for z in range(n + 1):
            b.append(a[z, :] * choose(n, z) * t[i] ** z * (1 - t[i]) **
            (n - z))
        c[i, :] = sum(b)
    p.plot(c[:, 0], c[:, 1])
    p.plot(a[:, 0], a[:, 1], 'ko')
    p.plot(a[:, 0], a[:, 1], 'k')
    return c[:, 0:2]


def rational_bezier(points, weights):
    n = len(points) - 1
    t = np.linspace(0, 1, 101)
    c = np.zeros([101, 2])
    for i in range(0, 101):
        b = []
        d = []
        for a in range(n + 1):
            b.append(points[a] * weights[a] * choose(n, a) * (t[i] ** a)\
            * ((1 - t[i]) ** (n - a)))
            d.append(weights[a] * choose(n, a) * (t[i] ** a) * ((1 - t[i]) \
            ** (n - a)))
        c[i, :] = (sum(b) / sum(d))
    p.plot(c[:, 0], c[:, 1])
    p.plot(points[:, 0], points[:, 1], 'ko')
    p.plot(points[:, 0], points[:, 1], 'k')
    return c[:, 0:2]