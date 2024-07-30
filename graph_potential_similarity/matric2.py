import csv
import numpy as np

class Matric2:
    def __init__(self, v_matrix_csv_file_name:str=None, v_matrix:np.ndarray=None):
        self.VComplete = None
        self.vMatLen = 0
        self.sum = 0
        self.sum_Square = 0
        self.sum_Cubic = 0

        self.EComplete = None
        self.eMatLen = 0
        self.m = 0
        self.m_Square = 0
        self.m_Cubic = 0

        if v_matrix_csv_file_name:
            self.VComplete = self.read_csv(v_matrix_csv_file_name)
            self.vMatLen = len(self.VComplete)
        elif not v_matrix is None:
            self.VComplete = v_matrix.tolist()
            self.vMatLen = len(self.VComplete)
        
        if not self.VComplete is None:
            v1s, v2s = [], []
            for i in range(self.vMatLen):
                for j in range(self.vMatLen):
                    # this is how we deal with general complex graphs: a complex graph allows more than one edges between two nodes
                    for k in range(self.VComplete[i][j]):
                        v1s.append(i)
                        v2s.append(j)

            self.eMatLen = len(v1s)
            self.EComplete = np.zeros((self.eMatLen, self.eMatLen), dtype=int)

            for k1 in range(len(v1s)):
                i1, j1 = v1s[k1], v2s[k1]
                for k2 in range(len(v1s)):
                    if k2 == k1:
                        continue
                    i2, j2 = v2s[k2], v2s[k2]
                    if i1 == i2 or i1 == j2 or j1 == i2 or j1 == j2:
                        self.EComplete[k1][k2] = 1

            v_sum_row = self.get_sum_row(self.VComplete)
            self.sum, self.sum_Square, self.sum_Cubic = self.get_ps(v_sum_row)

            e_sum_row = self.get_sum_row(self.EComplete)
            self.m, self.m_Square, self.m_Cubic = self.get_ps(e_sum_row)

    @staticmethod
    def get_sum_row(A):
        return [sum(1 for j in row if j > 0) for row in A]

    @staticmethod
    def get_ps(sum_rows):
        sorted_rows = sorted(sum_rows)
        ps = [0, 0, 0]
        for x in sorted_rows:
            x_sqr = x * x
            x_cubic = x_sqr * x
            ps[0] += x
            ps[1] += x_sqr
            ps[2] += x_cubic
        return ps

    @staticmethod
    def read_csv(file_path):
        with open(file_path, newline='') as csvfile:
            reader = csv.reader(csvfile)
            return [[int(value) for value in row] for row in reader]    

    def get_v_complete(self):
        return np.array(self.VComplete, dtype=float)

    def get_e_complete(self):
        return np.array(self.EComplete, dtype=float)

    def get_m(self):
        return self.m

    def set_m(self, m):
        self.m = m

    def set_m_square(self, m_square):
        self.m_Square = m_square

    def set_m_cubic(self, m_cubic):
        self.m_Cubic = m_cubic

    def get_m_square(self):
        return self.m_Square

    def get_m_cubic(self):
        return self.m_Cubic

    def get_sum(self):
        return self.sum

    def set_sum(self, sum):
        self.sum = sum

    def get_sum_square(self):
        return self.sum_Square

    def set_sum_square(self, sum_square):
        self.sum_Square = sum_square

    def get_sum_cubic(self):
        return self.sum_Cubic

    def set_sum_cubic(self, sum_cubic):
        self.sum_Cubic = sum_cubic

    def get_v_mat_len(self):
        return self.vMatLen

    def set_v_mat_len(self, v_mat_len):
        self.vMatLen = v_mat_len

    def get_e_mat_len(self):
        return self.eMatLen

    def set_e_mat_len(self, e_mat_len):
        self.eMatLen = e_mat_len

