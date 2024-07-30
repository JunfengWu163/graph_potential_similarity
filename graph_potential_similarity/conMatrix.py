import numpy as np
from scipy.linalg import svd

class ConMatrix:
    def __init__(self):
        self.sig = None  # Singular values
        self.U = None  # Left singular vectors
        self.V = None  # Right singular vectors
        self.sum = 0  # Sum of potential
        self.sum_sq = 0  # Sum of squares of potential
        self.sum_cu = 0  # Sum of cubes of potential
        self.map1 = {}  # Map to store indices of same singular values
        self.map2 = {}  # Map to store counts of points in stages S5-S8

    def cal_potential(self, m1, decimal_precision = 6):
        A = np.array(m1)

        # Perform Singular Value Decomposition
        U, s, Vt = svd(A)
        self.U = U
        self.V = Vt.T
        self.sig = np.round(s, decimal_precision)

        # Store indices of same singular values
        for i, value in enumerate(self.sig):
            if value not in self.map1:
                self.map1[value] = [i]
            else:
                self.map1[value].append(i)

        # Remove entries with less than 2 indices
        self.map1 = {k: v for k, v in self.map1.items() if len(v) >= 2}

        # Calculate potential vector S1
        S1 = []
        hm1 = {}
        for key, value in self.map1.items():
            nums = len(value)
            if nums > 1:
                if nums not in hm1:
                    hm1[nums] = nums
                else:
                    hm1[nums] += nums

        for key, value in hm1.items():
            for i in range(1, value + 1):
                S1.append(int(f"{key}{i}"))

        self.cal_rank_sum(S1)

    def cal_rank_sum(self, S):
        S1 = np.array(S, dtype=np.float64)
        ps = self.get_ps(S1)
        self.sum, self.sum_sq, self.sum_cu = ps

    @staticmethod
    def get_ps(sum_rows):
        sorted_rows = np.sort(sum_rows)
        ps = [0, 0, 0]
        for x in sorted_rows:
            x_sqr = x * x
            x_cubic = x_sqr * x
            ps[0] += x
            ps[1] += x_sqr
            ps[2] += x_cubic
        return ps

