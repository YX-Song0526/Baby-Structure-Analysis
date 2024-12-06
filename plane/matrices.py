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


def trans_mat_for_bar(phi):
    return np.array([[cos(phi), sin(phi), 0, 0],
                     [-sin(phi), cos(phi), 0, 0],
                     [0, 0, cos(phi), sin(phi)],
                     [0, 0, -sin(phi), cos(phi)]])


def trans_mat_for_frame(phi):
    return np.array([[cos(phi), sin(phi), 0, 0, 0, 0],
                     [-sin(phi), cos(phi), 0, 0, 0, 0],
                     [0, 0, 1, 0, 0, 0],
                     [0, 0, 0, cos(phi), sin(phi), 0],
                     [0, 0, 0, -sin(phi), cos(phi), 0],
                     [0, 0, 0, 0, 0, 1]])
