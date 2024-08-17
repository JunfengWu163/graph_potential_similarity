import unittest
import sys
import os
# Add the parent directory to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import graph_potential_similarity as gps

compound_smiles = {
    # 'CFCl2-CF2-O-CF2Cl': 'FC(Cl)(Cl)C(F)(F)OC(F)(F)Cl',
    # 'CF3-CCl2-O-CF2Cl': 'FC(F)(F)C(Cl)(Cl)OC(F)(F)Cl',
    # 'CF2H-CF2-O-CFH2': 'FC(F)C(F)(F)OCF',
    # 'CF3-CH(CF3)-O-CH3': 'FC(F)(F)C(C(F)(F)F)OC',
    # 'CF3-CH(CF3)-O-CF2H': 'FC(F)(F)C(C(F)(F)F)OC(F)',
    'Isoflurane (C3H2ClF5O)': 'C(C(F)(F)F)(OC(F)F)Cl',
    'Isoflurane R (C3H2ClF5O)': '[C@@H](C(OC(F)F)(F)F)(F)Cl',
    'Isoflurane S (C3H2ClF5O)': '[C@H](C(F)(F)F)(OC(F)F)Cl',
    'Tropacocaine (C15H19NO2)': 'CN1[C@@H]2CC[C@H]1CC(C2)OC(=O)C3=CC=CC=C3',
    'Cocaine (C17H21NO4)': 'CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)C3=CC=CC=C3)C(=O)OC',
    'Eucaine alpha (C15H21NO2)': 'C[C@H]1C[C@H](CC(N1)(C)C)OC(=O)C2=CC=CC=C2',
    'Eucaine beta (C15H21NO2)': 'CC1CC(CC(N1)(C)C)OC(=O)C2=CC=CC=C2'
}

class TestSMILES(unittest.TestCase):
    def main(self):
        input_names = sorted(compound_smiles.keys())
        n = len(input_names)
        output_matrix = [[0 for _ in range(n)] for _ in range(n)]
        for i in range(n):
            output_matrix[i][i] = 1.0
            for j in range(i + 1, n):
                dist = gps.smiles_distance(compound_smiles[input_names[i]],compound_smiles[input_names[j]], dist_type='hop_tree')
                sim = 1.0 - dist
                output_matrix[i][j] = sim
                output_matrix[j][i] = sim
                print(f'{i},{j}:{sim}')
                # diff = abs(sim - float(sims[i][j]))
                # self.assertLessEqual(diff, 1e-10)
        for i in range(len(input_names)):
            print(f'{input_names[i]}: {output_matrix[i]}')

if __name__ == '__main__':
    test = TestSMILES()
    test.main()