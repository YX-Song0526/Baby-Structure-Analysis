import numpy as np
from solver.skyline import skyline_storage


def ldl_raw(K):
    n = len(K)
    L = np.eye(n)
    D = np.zeros(n)

    for j in range(n):
        for i in range(j):
            L[j, i] = (K[j, i] - np.sum(L[i, :i] * D[:i] * L[j, :i])) / D[i]
        D[j] = K[j, j] - np.sum(L[j, :j] * D[:j] * L[j, :j])

    return L, D


def ldl_modified(K):
    n = len(K)
    L = np.eye(n)
    D = np.zeros(n)
    g = np.zeros((n, n))

    for j in range(n):
        for i in range(j):
            g[j, i] = K[j, i] - np.sum(L[i, :i] * g[j, :i])
        D[j] = K[j, j] - np.sum(g[j, :j] ** 2 / D[:j])
        for i in range(j):
            L[j, i] = g[j, i] / D[i]
    return L, D


def ldl_with_skyline(K):
    n = len(K)
    L = np.eye(n)
    D = np.zeros(n)

    K_sky, pointer, edge = skyline_storage(K)

    for j in range(n):
        for i in range(edge[j], j):
            L[j, i] = (K[j, i] - np.sum(L[i, :i] * D[:i] * L[j, :i])) / D[i]
        D[j] = K[j, j] - np.sum(L[j, :j] * D[:j] * L[j, :j])

    return L, D


def ldl_modified_with_skyline(K):
    n = len(K)

    L = np.eye(n)
    D = np.zeros(n)
    g = np.zeros((n, n))

    K_sky, pointer, edge = skyline_storage(K)

    for j in range(n):
        # 求g_ij
        for i in range(edge[j], j):
            g[j, i] = K[j, i] - np.sum(L[i, :i] * g[j, :i])
            for k in range(max(edge[i], edge[j]), i):
                K_sky[pointer[j] - (j - i)] -= K_sky[pointer[i] - (i - k)] * K_sky[pointer[j] - (j - k)]

        # 求D_ij
        D[j] = K[j, j] - np.sum(g[j, edge[j]:j] ** 2 / D[edge[j]:j])

        for k in range(edge[j], j):
            K_sky[pointer[j]] -= K_sky[pointer[j] - (j - k)] ** 2 / K_sky[pointer[k]]

        # 求L_ij
        for i in range(edge[j], j):
            L[j, i] = g[j, i] / D[i]
            K_sky[pointer[j] - (j - i)] = K_sky[pointer[j] - (j - i)] / K_sky[pointer[i]]

    return L, D, K_sky
