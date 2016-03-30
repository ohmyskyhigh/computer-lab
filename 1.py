import pylab as p
import numpy as np

def rational_bezier(points, weights):
    n = len(points) - 1
    m = len(weight) - 1
    t = np.linspace(0, 1, 101)
    c = np.zeros([101, 2])
    for i in range(0, 101):
        b = []
        d = []
        for a in range(n + 1):
            b.append(points[a, :] * weights[e, :] * choose(n, a) * t[i] ** a *
            (1 - t[i]) ** (n - a)
            