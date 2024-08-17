import unittest
import sys
import os
# Add the parent directory to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import graph_potential_similarity as gps

seqs = [
    "GGACCGACAGGAATTCGCTCCTTAGGACGTTATAGTTACGGCCGCCGTTTACTGGGGCTTCAATTCGCAGCTTCGC",
    "AGCGAACGCTGGCGGCATGCTTAACACATGCAAGTCGCACGAAGGCTTCGGCCTTAGTGGCGGACGGGTGAGTAAC",
    "AGCGAACGCTGGCGGCATGCTTAACACATGCAAGTCGCACGAAGGCTTCGGCCTTAGTGGCGGACGGGTGAGTAAC",
    "TAGTGGCGGACGGGTGAGTAACGCGTAGGAATCTATCCATGGGTGGGGGATAACTCCGGGAAACTGGAGCTAATAC",
    "GTGGCGGACGGGTGAGTAACGCGTAGGAATCTATCCATGGGTGGGGGATAACTCCGGGAAACTGGAGCTAATACCG"
]
dists = [
    [0, 0.059118503933023835, 0.059118503933023835, 0.045407486149992286, 0.04649854198482426], 
    [0.059118503933023835, 0, 0.0, 0.0628405446069266, 0.05430110467820465], 
    [0.059118503933023835, 0.0, 0, 0.0628405446069266, 0.05430110467820465], 
    [0.045407486149992286, 0.0628405446069266, 0.0628405446069266, 0, 0.05023724897819748], 
    [0.04649854198482426, 0.05430110467820465, 0.05430110467820465, 0.05023724897819748, 0]
    ]
class TestDoubleHelix(unittest.TestCase):
    def main(self):
        n = len(seqs)
        output_matrix = [["" for _ in range(n)] for _ in range(n)]
        for i in range(n):
            output_matrix[i][i] = 0
            for j in range(i + 1, n):
                dist = gps.double_helix_distance(seqs[i],seqs[j],dist_type='hop_tree')
                output_matrix[i][j] = dist
                output_matrix[j][i] = dist
                print(f'{i},{j}:{dist}')
                # diff = abs(dist - dists[i][j])
                # self.assertLessEqual(diff, 1e-10)
        for i in range(n):
            print(f'{i}: {output_matrix[i]}')

if __name__ == '__main__':
    test = TestDoubleHelix()
    test.main()