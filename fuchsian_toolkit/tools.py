import time
from collections import deque
from copy import deepcopy
from functools import reduce
from typing import Any, List, Optional

import sympy as sp
from sympy import I, Matrix, S, eye, zeros

z = sp.symbols("z", complex=True)


class TMatrix:
    def __init__(self, size: int, index: int, s: complex, delta: int):
        self._s = s

        self._index = index
        self._size = size

        self._mat = sp.eye(size)
        self._mat_inv = sp.eye(size)

        self._mat[index, index] = (z - s) ** delta
        self._mat_inv[index, index] = (z - s) ** (-delta)

    def get_mat(self):
        return self._mat[:, :]

    def _zero_arg(self):
        m = sp.eye(self._size)
        m[self._index, self._index] = 0
        return m

    def __call__(self, x: complex):
        if x == self._s:
            return self._zero_arg()

        return self._mat.subs({z: x})

    def inv(self, x: complex):
        if x == self._s:
            return self._zero_arg()

        return self._mat_inv.subs({z: x})


class FuchKit:
    def __init__(self, mu_list, s_list, A_list):
        self.size = A_list[0].shape[0]

        self.steps = []

        self.T_list = []

        self.A_list = deepcopy(A_list)
        self.s_list = s_list[:]
        self.mu_list = mu_list = [sorted(m) for m in mu_list]

        self.C = zeros(self.size)

        self.steps.append((eye(self.size), self.A_list, self.C))

    def rotate(self):
        self.A_list = self.A_list[1:] + [self.A_list[0]]
        self.s_list = self.s_list[1:] + [self.s_list[0]]
        self.mu_list = self.mu_list[1:] + [self.mu_list[0]]

    def get_transform_matrix(self) -> Matrix:
        return sp.prod(self.T_list[::-1]).applyfunc(sp.simplify)

    @staticmethod
    def get_floor_diag(A: Matrix):
        return [sp.floor(sp.re(a)) for i, a in enumerate(A.diagonal())]

    @staticmethod
    def res(A, p):
        return A.applyfunc(lambda e: sp.residue(e, z, p))

    def transform_general(self, T: Matrix):
        T_inv = T.inv()

        size = T.shape[0]

        C_hat = T * self.C * T_inv + sp.diff(T, z) * T_inv

        B = sum(
            ((T * A * T_inv) / (z - s) for s, A in zip(self.s_list, self.A_list)), C_hat
        ).applyfunc(sp.simplify)

        A_list = [self.res(B, s) for s in self.s_list]
        C = sum((-A / (z - s) for s, A in zip(self.s_list, A_list)), B).applyfunc(
            sp.simplify
        )

        self.A_list = A_list
        self.C = C

    def transform_const(self, H: Matrix):
        H_inv = H.inv()

        self.T_list.append(H)

        self.A_list = [(H * A * H_inv).applyfunc(sp.simplify) for A in self.A_list]
        self.C = (H * self.C * H_inv).applyfunc(sp.simplify)

    def transform_jordan_form(self):
        A = self.A_list[0]
        mu = self.mu_list[0]

        size = self.size

        T_inv, J = A.jordan_form()
        T = T_inv.inv()

        eigenvals = self.get_floor_diag(J)
        indices = sorted(enumerate(eigenvals), key=lambda x: x[1])

        change_action = 0

        for m, (row_index, e) in zip(mu, indices):
            if m != e:
                change_action = 1 if m > e else -1
                break

        if change_action == 1:
            col = J[:, row_index].flat()
            col[row_index] = 0

            if any(x != 0 for x in col):
                H = eye(size)
                H.row_swap(row_index - 1, row_index)
                T = H * T

        elif change_action == -1:
            row = J[row_index, :].flat()
            row[row_index] = 0

            if any(x != 0 for x in row):
                H = eye(size)
                H.row_swap(row_index, row_index + 1)
                T = H * T

        self.transform_const(T)

    def transform_add_one(self, row_ind: int):
        size = self.size

        s = self.s_list[0]

        T = TMatrix(size, row_ind, s, 1)

        A_hat_list = [T(s) * A * T.inv(s) for s, A in zip(self.s_list, self.A_list)]

        A_hat = A_hat_list[0]

        A_hat[row_ind, row_ind] = self.A_list[0][row_ind, row_ind] + 1

        C_hat = self.C[:, :]

        for i in range(size):
            if i != row_ind:
                c_bar, c_tilde = sp.div(self.C[i, row_ind], z - s)

                A_hat[i, row_ind] += sum(
                    (
                        -A_j[i, row_ind] / (s_j - s)
                        for s_j, A_j in zip(self.s_list[1:], self.A_list[1:])
                    ),
                    c_tilde,
                )

                C_hat[i, row_ind] = c_bar
                C_hat[row_ind, i] *= z - s

                C_hat[row_ind, i] += sum(A_j[row_ind, i] for A_j in self.A_list)

        self.A_list = [A.applyfunc(sp.simplify) for A in A_hat_list]
        self.C = C_hat

        self.T_list.append(T.get_mat())

    def transform_sub_one(self, row_ind: int):
        size = self.size

        s = self.s_list[0]

        T = TMatrix(size, row_ind, s, -1)

        A_hat_list = [T(s) * A * T.inv(s) for s, A in zip(self.s_list, self.A_list)]

        A_hat = A_hat_list[0]

        A_hat[row_ind, row_ind] = self.A_list[0][row_ind, row_ind] - 1

        C_hat = self.C[:, :]

        for i in range(size):
            if i != row_ind:
                c_bar, c_tilde = sp.div(self.C[row_ind, i], z - s)

                A_hat[row_ind, i] += sum(
                    (
                        -A_j[row_ind, i] / (s_j - s)
                        for s_j, A_j in zip(self.s_list[1:], self.A_list[1:])
                    ),
                    c_tilde,
                )

                C_hat[row_ind, i] = c_bar
                C_hat[i, row_ind] *= z - s

                C_hat[i, row_ind] += sum(A_j[i, row_ind] for A_j in self.A_list)

        self.A_list = A_hat_list
        self.C = C_hat

        self.T_list.append(T.get_mat())

    def make_transformations(self):
        count_steps = 0

        st = time.perf_counter()

        def increase_step():
            nonlocal count_steps
            count_steps += 1
            print(f"step: {count_steps}")

        for ind, mu in enumerate(self.mu_list):
            print("cycle:", ind + 1)

            increase_step()

            self.transform_jordan_form()

            eigenvals = self.get_floor_diag(self.A_list[0])

            indices = sorted(enumerate(eigenvals), key=lambda x: x[1])

            for i, m in enumerate(mu):
                while indices[i][1] > m:
                    increase_step()

                    self.transform_sub_one(indices[i][0])

                    eigenvals = self.get_floor_diag(self.A_list[0])
                    indices = sorted(enumerate(eigenvals), key=lambda x: x[1])

                while indices[i][1] < m:
                    increase_step()

                    self.transform_add_one(indices[i][0])

                    eigenvals = self.get_floor_diag(self.A_list[0])
                    indices = sorted(enumerate(eigenvals), key=lambda x: x[1])

            self.rotate()

        ed = time.perf_counter()

        print("total time:", ed - st)


def get_fuchsian_transform(
    A_list: List[Matrix], s_list: List[complex], mu_list: List[List[int]]
):
    fuch = FuchKit(mu_list, s_list, A_list)

    fuch.make_transformations()

    T = fuch.get_transform_matrix()

    return T, fuch


# Examples


def example_0():
    A_list = [
        Matrix([[S("1"), S("0")], [S("-3/4"), S("2")]]),
        Matrix([[S("-3/2"), S("0")], [S("1/2"), S("-3/2")]]),
        Matrix([[S("1/2"), S("0")], [S("1/4"), S("-1/2")]]),
    ]

    s_list = [2, 0, 1]

    mu_list = [
        [0, 0],
        [0, 0],
        [0, 0],
    ]

    T, fuch = get_fuchsian_transform(A_list, s_list, mu_list)

    for A in fuch.A_list:
        print(A)

    print()

    print(fuch.C)

    print()

    print(T)
