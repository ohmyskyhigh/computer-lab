import pylab as p
import numpy as np
from scipy.optimize import fsolve
file_path = 'C:\\applycomputing\\'
xfoil_path = 'C:\\applycomputing\\'


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


def parametric_aerofoil(w, file_path):
    pi = np.array([[1.0, 0.0], [0.5, 0.08], [0.0, -0.05]])
    qu = np.array([[0.0, 0.1], [0.4, 0.2], [1.0, 0.0]])
    zp = np.array([1, 1, 1, 1])
    zq = np.array([1, 1, w, 1])
    print zq
    n = np.float(pi.shape[0])
    m = np.float(pi.shape[0])
    q_start = p_end = (n / (m + n) * pi[-1, :] + (m / (m + n)) * qu[0, :])
    pp = np.vstack([pi, p_end])
    qq = np.vstack([q_start, qu])
    p.plot(pp[:, 0], pp[:, 1], 'ko')
    p.plot(qq[:, 0], qq[:, 1], 'ro')
    lower = rational_bezier(pp, zp)
    upper = rational_bezier(qq, zq)
    data_file = open(file_path + 'aerofoil.dat', 'w')
    data_file.write('OurFoil\n')
    for i in range(0, 100):
        data_file.write("%.18f %.18f\n" %(lower[i, 0], lower[i, 1]))
    for i in range(0, 101):
        data_file.write("%.18f %.18f\n" %(upper[i, 0], upper[i, 1])) 
    data_file.close()
    

def run_xfoil_wcl(w, cl, file_path, xfoil_path):
    a = 340.3
    v = 12.5
    nu = 0.00001461
    L = 0.698
    Re = v * L / nu
    M = v / a
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
    iter\n 5000\n\
    cl ' + str(cl) +  '\n\
    \n\
    \n\
    quit\n')
    command_file.close()
    import os
    parametric_aerofoil(w, file_path)
    run_xfoil_command = '\"' + xfoil_path + 'xfoil.exe\" < ' + file_path + 'commands.in > '+ file_path + 'dump.out'
    os.system(run_xfoil_command)
    aero_data_file = open(file_path + 'polar.dat', 'r')
    oin = aero_data_file.readlines()
    aero_data_file.close()
    os.system('del ' + file_path + 'polar.dat')
    cd = np.float(oin[-1][20: 27])
    cl = np.float(oin[-1][11: 17])
    print cd
    return cl, cd


def parameter_sweep(w_array, cl, file_path, xfoil_path):
    cd = np.zeros(len(w_array))
    for i in range(len(w_array)):
        [cl, cd[i]] = (run_xfoil_wcl(w_array[i], cl, file_path, xfoil_path))
    return cd

    
def mls(x, X, y, sigma):
    N = max(np.shape(y))
    weights = np.zeros(N)
    A = np.zeros([N, 3])
    A[:, 0] = X ** 2
    A[:, 1] = X
    A[:, 2] = np.ones([1, N])
    for i in range(1, N):
        weights[i] = np.exp(-(np.sum((x - X[i]) ** 2)) / (2 * sigma))
    W = np.diag(weights)
    a = np.linalg.lstsq(np.dot(np.dot(A.conj().T, W), A), np.dot(np.dot(A.conj().T, W), y))
    f = a[0][0] * x ** 2 + a[0][1] * x + a[0][2]
    return f


def mls_error(sigma, X, y):
    y_test = np.zeros(11)
    error = np.zeros(11)
    for i in range(0, 11):
        y_test[i] = mls(X[i], np.append(X[0: i], X[i + 1: -1]), np.append(y[0: i], y[i + 1: -1]), sigma)
        error[i] = (y[i] - y_test[i]) ** 2
    sum_error = np.sum(error)
    return sum_error


def mls_curve_fit(w_array, cd, x):
    sigma_best = fsolve(mls_error, 0.5, args = (w_array, cd))
    w_fine = np.linspace(0.6, 1.2, 101)
    y_pred = np.zeros(101)
    for i in range(0, 101):
        y_pred[i] = mls(w_fine[i], w_array, cd, sigma_best)
    p.plot(w_array, cd, 'o', label = 'XFOIL_data')
    p.plot(w_fine, y_pred, label = 'MLS_fit')
    p.legend()
    f = mls(x, w_array, cd, sigma_best) 
    return f 
    

def linspace(start, end, interval):
    a = np.linspace(start, end, interval)
    return a.tolist()