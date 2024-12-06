import numpy as np
from numpy import cos, sin


def little_mat(phi):
    return np.array([[cos(phi) ** 2, cos(phi) * sin(phi), -cos(phi) ** 2, -cos(phi) * sin(phi)],
                     [cos(phi) * sin(phi), sin(phi) ** 2, -cos(phi) * sin(phi), -sin(phi) ** 2],
                     [-cos(phi) ** 2, -cos(phi) * sin(phi), cos(phi) ** 2, cos(phi) * sin(phi)],
                     [-cos(phi) * sin(phi), -sin(phi) ** 2, cos(phi) * sin(phi), sin(phi) ** 2]])


def tiny_mat(L):
    return np.array([[12, 6 * L, -12, 6 * L],
                     [6 * L, 4 * L ** 2, -6 * L, 2 * L ** 2],
                     [-12, -6 * L, 12, -6 * L],
                     [6 * L, 2 * L ** 2, -6 * L, 4 * L ** 2]])


def big_mat(E, A, I, L):
    return np.array([[E * A / L, 0, 0, -E * A / L, 0, 0],
                     [0, 12 * E * I / L ** 3, 6 * E * I / L ** 2, 0, -12 * E * I / L ** 3, 6 * E * I / L ** 2],
                     [0, 6 * E * I / L ** 2, 4 * E * I / L, 0, -6 * E * I / L ** 2, 2 * E * I / L],
                     [-E * A / L, 0, 0, E * A / L, 0, 0],
                     [0, -12 * E * I / L ** 3, -6 * E * I / L ** 2, 0, 12 * E * I / L ** 3, -6 * E * I / L ** 2],
                     [0, 6 * E * I / L ** 2, 2 * E * I / L, 0, -6 * E * I / L ** 2, 4 * E * I / L]])


def trans_mat(phi):
    return np.array([[cos(phi), sin(phi), 0, 0, 0, 0],
                     [-sin(phi), cos(phi), 0, 0, 0, 0],
                     [0, 0, 1, 0, 0, 0],
                     [0, 0, 0, cos(phi), sin(phi), 0],
                     [0, 0, 0, -sin(phi), cos(phi), 0],
                     [0, 0, 0, 0, 0, 1]])


def augment(E: float, A: float, G: float, J: float, I: list, L: float):
    Iy, Iz = I
    K_bar = (E * A / L) * np.array([[1, -1],
                                    [-1, 1]])

    K_G = (G * J / L) * np.array([[1, -1],
                                  [-1, 1]])

    K_Bxy = (E * Iz / L ** 3) * np.array([[12, 6 * L, -12, 6 * L],
                                          [6 * L, 4 * L ** 2, -6 * L, 2 * L ** 2],
                                          [-12, -6 * L, 12, -6 * L],
                                          [6 * L, 2 * L ** 2, -6 * L, 4 * L ** 2]])

    K_Bxz = (E * Iy / L ** 3) * np.array([[12, 6 * L, -12, 6 * L],
                                          [6 * L, 4 * L ** 2, -6 * L, 2 * L ** 2],
                                          [-12, -6 * L, 12, -6 * L],
                                          [6 * L, 2 * L ** 2, -6 * L, 4 * L ** 2]])

    Ke = np.zeros((12, 12))

    for i in range(2):
        for j in range(2):
            Ke[6 * i, 6 * j] += K_bar[i, j]
            Ke[6 * i + 3, 6 * j + 3] += K_G[i, j]

    for i in range(4):
        for j in range(4):
            dof = [1, 5, 7, 11]
            Ke[dof[i], dof[j]] += K_Bxy[i, j]
            dof = [2, 4, 8, 10]
            Ke[dof[i], dof[j]] += K_Bxz[i, j]

    return Ke


def T_mat(dx, dy, dz):
    L = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    lx, mx, nx = dx / L, dy / L, dz / L
    local_x = np.array([lx, mx, nx])

    # print(lx, mx, nx)

    if np.allclose([lx, mx, nx], [0.0, 0.0, 1.0]) or np.allclose([lx, mx, nx], [0.0, 0.0, -1.0]):
        # Handle case where direction vector is aligned with global y-axis
        local_y = np.array([0, 1, 0])
        # print("Hello")
    else:
        global_z = np.array([0, 0, 1])
        local_y = np.cross(global_z, local_x)
        local_y /= np.linalg.norm(local_y)

    # Compute local z-axis
    local_z = np.cross(local_x, local_y)

    R = np.array([local_x,
                  local_y,
                  local_z])

    T = np.zeros((12, 12))
    for i in range(4):
        T[i * 3:(i + 1) * 3, i * 3:(i + 1) * 3] = R

    return T
