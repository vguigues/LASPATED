#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import numpy as np

C = 3
D = 7
T = 48
R = 297
J = 4


def read_file(path):
    arq = open(path, "r")
    arq.readline()
    x = np.zeros((C, D, T, J))
    grad = np.zeros((C, D, T, J))
    xmz = np.zeros((C, D, T, J))
    for r in range(C * D * T * J):
        tokens = arq.readline().split()
        c = int(tokens[0][1])
        d = int(tokens[1][1])
        t = int(tokens[2][1])
        j = int(tokens[3][1])
        x[c, d, t, j] = float(tokens[7])
        grad[c, d, t, j] = float(tokens[10])
        xmz[c, d, t, j] = float(tokens[13])

    rhs = float(arq.readline().split()[2])
    arq.close()

    return x, grad, xmz, rhs


x1, g1, xmz1, rhs1 = read_file("derivatives_matlab.txt")
x2, g2, xmz2, rhs2 = read_file("derivatives_cpp.txt")

print("1 = matlab, 2 = cpp")
print("norm (x1-x2) =", np.linalg.norm(x1 - x2))
print("norm (g1-g2) =", np.linalg.norm(g1 - g2))
print("norm (xmz1-xmz2) =", np.linalg.norm(xmz1 - xmz2))
print("rhs1 =", rhs1)
print("rhs2 =", rhs2)
print("rhs1 - rhs2 =", abs(rhs1 - rhs2))


py_rhs1 = 0
py_rhs2 = 0
diff_sum = 0
precision = 10**-3
for c in range(C):
    for d in range(D):
        for t in range(T):
            for j in range(J):
                py_rhs1 += g1[c, d, t, j] * xmz1[c, d, t, j]
                py_rhs2 += g2[c, d, t, j] * xmz2[c, d, t, j]
                if abs(g1[c, d, t, j] - g2[c, d, t, j]) > precision:
                    diff_sum += abs(g1[c, d, t, j] - g2[c, d, t, j])
                    print(f"{c},{d},{t},{j}: g1 = {g1[c,d,t,j]}, g2={g2[c, d, t, j]}")

print(f"diff_sum  greater than {precision}  = {diff_sum}")
