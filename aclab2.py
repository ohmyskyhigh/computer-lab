import pylab as p
import numpy as np
from scipy.optimize import fsolve


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
    

def bezier_spline_aerofoil(file_path):
    pi = np.array([[1.0, 0.0], [0.5, 0.08], [0.0, -0.05]])
    qu = np.array([[0.0, 0.1], [0.4, 0.2], [1.0, 0.0]])
    zp = np.array([1, 1, 1, 1])
    zq = np.array([1, 1, 1, 1])
    n = np.float(pi.shape[0])
    m = np.float(pi.shape[0])
    q_start = p_end = (n / (m + n) * pi[-1, :] + (m / (m + n)) * qu[0, :])
    pp = np.vstack([pi, p_end])
    qq = np.vstack([q_start, qu])
    p.plot(pp[:, 0], pp[:, 1], 'ko')
    p.plot(qq[:, 0], qq[:, 1], 'ro')
    lower = rational_bezier(pp, zp)
    upper = rational_bezier(qq, zq)
    p.show()
    data_file = open(file_path + 'aerofoil.dat', 'w')
    data_file.write('OurFoil\n')
    for i in range(0, 100):
        data_file.write("%f %f\n" %(lower[i, 0], lower[i, 1]))
    for i in range(0, 101):
        data_file.write("%f %f\n" %(upper[i, 0], upper[i, 1])) 
    data_file.close()


def run_xfoil(file_path, xfoil_path):
    Re = 1000
    M = 0.5
    command_file=open(file_path + 'commands.in', 'w')
    command_file.write('load ' + file_path + 'aerofoil.dat\n\
    panel\n\
    oper\n\
    visc ' + str(Re) + '\n\
    M ' + str(M) + '\n\
    type 1\n\
    pacc\n'\
    + file_path + 'polar.dat\n\
    \n\
    iter\n 1000\n\
    cl 1.2\n\
    \n\
    \n\
    quit\n')
    command_file.close()
    import os
    run_xfoil_command = '\"' + xfoil_path + 'xfoil.exe\" < ' + file_path + 'commands.in > '+ file_path + 'dump.out'
    os.system(run_xfoil_command)
    print(run_xfoil_command)
    aero_data_file = open(file_path + 'polar.dat', 'r')
    oin = aero_data_file.readlines()
    aero_data_file.close()
    os.system('del ' + file_path + 'polar.dat')
    cl = np.float(oin[-1][11: 17])
    cd = np.float(oin[-1][20: 27])
    print cl, cd

