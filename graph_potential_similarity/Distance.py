from .matric2 import Matric2
from .conMatrix import ConMatrix
from .gauss import gauss
from .RrefResult import RrefResult
from .calShiNum import cal_shi_num
from .directed_potential_graph import directed_potential_graph
from rdkit import Chem
import numpy as np

class Distance:
    Vcm1 = None
    Vcm2 = None
    Ecm1 = None
    Ecm2 = None

    @staticmethod
    def cal(m1, m2, dist_type = 'original', decimal_precision = 6):
        if dist_type == 'hop_tree':
            return 1.0 - Distance.cal_hop_tree_sim(m1.get_v_complete(), m2.get_v_complete())
        
        Distance.Vcm1 = ConMatrix()
        Distance.Vcm2 = ConMatrix()
        Distance.Ecm1 = ConMatrix()
        Distance.Ecm2 = ConMatrix()
        Distance.Vcm1.cal_potential(m1.get_v_complete(), decimal_precision)
        Distance.Vcm2.cal_potential(m2.get_v_complete(), decimal_precision)
        Distance.Ecm1.cal_potential(m1.get_e_complete(), decimal_precision)
        Distance.Ecm2.cal_potential(m2.get_e_complete(), decimal_precision)

        if dist_type == 'continuous':
            return (Distance.calS1(m1, m2) + Distance.calS2() + Distance.calS3()) / 3.0
        else:
            distance = 0
            if (distance := Distance.calS1(m1, m2)) != 0.0:
                return distance
            elif (distance := Distance.calS2(m1, m2, decimal_precision)) != 0.0:
                return distance
            elif (distance := Distance.calS3()) != 0.0:
                return distance
            return 0.0
    
    @staticmethod
    def cal_hop_tree_sim(m1,m2):
        g1 = directed_potential_graph(m1,max_hop=2)
        g2 = directed_potential_graph(m2,max_hop=2)
        return g1.get_similarity(g2)
    
    @staticmethod
    def get_graph_matrix(mol, atom_refs, ref_ids):
        n1 = len(atom_refs)
        n2 = len(ref_ids)
        n = n1 + n2
        graph_matrix = np.zeros(shape=(n, n),dtype=np.int8)
        graph_matrix[:n1,:n1] = Chem.GetAdjacencyMatrix(mol)
        for i in range(n1):
            graph_matrix[i,n1+atom_refs[i]] = 1
            graph_matrix[n1+atom_refs[i],i] = 1
        for i in range(n1,n-1):
            graph_matrix[i,i+1] = 1
            graph_matrix[i+1,i] = 1
        return graph_matrix

    @staticmethod
    def from_smiles(smiles1, smiles2, dist_type = 'continuous', decimal_precision = 6):
        print('creating graphs')
        # get molecules from SMILES strings
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        
        # find the union of their atom sets
        atoms = {}
        for atom in mol1.GetAtoms():
            atoms[atom.GetAtomicNum()] = atom.GetSymbol()
        for atom in mol2.GetAtoms():
            atoms[atom.GetAtomicNum()] = atom.GetSymbol()
        
        # convert the union to ref_ids
        atom_list = sorted(atoms.items())
        ref_ids = {}
        for i in range(len(atom_list)):
            ref_ids[atom_list[i][1]] = i

        # get atom refs
        atom_refs_1 = []
        for atom in mol1.GetAtoms():
            atom_refs_1.append(ref_ids[atom.GetSymbol()])
        atom_refs_2 = []
        for atom in mol2.GetAtoms():
            atom_refs_2.append(ref_ids[atom.GetSymbol()])
        
        # get graph matrices:
        m1 = Matric2(v_matrix=Distance.get_graph_matrix(mol1, atom_refs_1, ref_ids))
        m2 = Matric2(v_matrix=Distance.get_graph_matrix(mol2, atom_refs_2, ref_ids))

        # compute distance
        print('computing distasnce')
        return Distance.cal(m1, m2, dist_type, decimal_precision)

    @staticmethod
    def getGCATN(seqIn):
        seqOut = ''
        for char in seqIn:
            if char == 'G' or char == 'g':
                seqOut = seqOut + 'G'
            elif char == 'C' or char == 'c':
                seqOut = seqOut + 'C'
            elif char == 'A' or char == 'a':
                seqOut = seqOut + 'A'
            elif char == 'T' or char == 't':
                seqOut = seqOut + 'T'
            elif char == ' ':
                continue
            else:
                seqOut = seqOut + 'N'
        return seqOut
    
    @staticmethod
    def get_double_helix_mat(seq):
        n = len(seq)
        zlmap = {'G': 'C', 'C': 'G', 'A': 'T', 'T': 'A', 'N': 'N'}

        jjmap = {'G': 2*n, 'C': 2*n + 1, 'A': 2*n + 2, 'T': 2*n + 3, 'N': 2*n + 4}

        edges = []
        # the left chain of seq (0 to n-1), connects n-1 edges
        for i in range(n-1):
            edges.append((i,i+1))
        # the right chain of seq (n to 2n-1), connects n-1 edges
        for i in range(n, 2*n - 1):
            edges.append((i,i+1))
        # G(2n-1) T(2n) A(2n+1) C(2n+2) N(2n+3) (2n to 2n + 4), connects 4 edges
        for i in range(2*n, 2*n + 4):
            edges.append((i,i+1))

        for i in range(n):
            # compute the connections of the left chain with GTACN
            target = jjmap[seq[i]]
            source = i
            edges.append((source, target))

            # compute the connections of the right chain with GTACN
            target = jjmap[zlmap[seq[i]]]
            source = n + i
            edges.append((source ,target))

            # compute the connections between the left and right chain
            source = i
            target = n + i
            edges.append((source ,target))
        
        mat = np.zeros(shape=(2*n+5,2*n+5),dtype=np.int8)
        for edge in edges:
            mat[edge[0],edge[1]] = 1
            mat[edge[1],edge[0]] = 1
        return mat


    @staticmethod
    def from_double_helix(seq1, seq2, dist_type = 'continuous', decimal_precision = 6):
        # get graph matrices:
        print('creating graphs')
        m1 = Matric2(v_matrix=Distance.get_double_helix_mat(seq1))
        m2 = Matric2(v_matrix=Distance.get_double_helix_mat(seq2))

        # compute distance
        print('computing distasnce')
        return Distance.cal(m1, m2, dist_type, decimal_precision)
    
    # A complex graph refers to a graph with cycles, multi-edges, and self-loops.
    # The adjacency matrix of a complex graph with n nodes is of size n * n, and the (i,j) entry is a non-negative integer representing the number of directed edges 
    # from node i to node j. The matrix needs not to be symetric. 
    @staticmethod
    def from_complex_graphs(adjMat1, adjMat2, dist_type = 'continuous', decimal_precision = 6):
        # get graph matrices:
        m1 = Matric2(v_matrix=adjMat1)
        m2 = Matric2(v_matrix=adjMat2)

        # compute distance
        return Distance.cal(m1, m2, dist_type, decimal_precision)

    @staticmethod
    def calS1(m1, m2):
        x1, y1, z1 = m1.get_sum(), m1.get_sum_square(), m1.get_sum_cubic()
        i1, j1, k1 = m1.get_m(), m1.get_m_square(), m1.get_m_cubic()
        x2, y2, z2 = m2.get_sum(), m2.get_sum_square(), m2.get_sum_cubic()
        i2, j2, k2 = m2.get_m(), m2.get_m_square(), m2.get_m_cubic()

        dx, dy, dz = x2 - x1, y2 - y1, z2 - z1
        di, dj, dk = i2 - i1, j2 - j1, k2 - k1

        c = dx * dx + dy * dy + dz * dz
        d = x1 * x1 + x2 * x2 + y1 * y1 + y2 * y2 + z1 * z1 + z2 * z2
        VD = np.sqrt(c) / np.sqrt(d)

        c = di * di + dj * dj + dk * dk
        d = i1 * i1 + i2 * i2 + j1 * j1 + j2 * j2 + k1 * k1 + k2 * k2
        ED = np.sqrt(c) / np.sqrt(d)
        return 0.5 * VD + 0.5 * ED

    @staticmethod
    def calS2():
        if not Distance.Vcm1.map1 or not Distance.Vcm2.map1:
            VDis = 0
        else:
            x1, y1, z1 = Distance.Vcm1.sum, Distance.Vcm1.sum_sq, Distance.Vcm1.sum_cu
            x2, y2, z2 = Distance.Vcm2.sum, Distance.Vcm2.sum_sq, Distance.Vcm2.sum_cu
            dx, dy, dz = x2 - x1, y2 - y1, z2 - z1
            c = dx * dx + dy * dy + dz * dz
            d = x1 * x1 + x2 * x2 + y1 * y1 + y2 * y2 + z1 * z1 + z2 * z2
            VDis = np.sqrt(c) / np.sqrt(d)

        if not Distance.Ecm1.map1 or not Distance.Ecm2.map1:
            EDis = 0
        else:
            x1, y1, z1 = Distance.Ecm1.sum, Distance.Ecm1.sum_sq, Distance.Ecm1.sum_cu
            x2, y2, z2 = Distance.Ecm2.sum, Distance.Ecm2.sum_sq, Distance.Ecm2.sum_cu
            dx, dy, dz = x2 - x1, y2 - y1, z2 - z1
            c = dx * dx + dy * dy + dz * dz
            d = x1 * x1 + x2 * x2 + y1 * y1 + y2 * y2 + z1 * z1 + z2 * z2
            EDis = np.sqrt(c) / np.sqrt(d)

        return 0.25 * VDis + 0.25 * EDis

    @staticmethod
    def calS3():
        res = Distance.cal_rank(Distance.Vcm1, Distance.Vcm2, "V") + Distance.cal_rank(Distance.Ecm1, Distance.Ecm2, "E")
        return res

    @staticmethod
    def cal_rank(cm1, cm2, VorE):
        flagU = True
        flagV = True

        map1 = iter(cm1.map1.items())
        map2 = iter(cm2.map1.items())

        while True:
            try:
                entry1 = next(map1)
                entry2 = next(map2)
            except StopIteration:
                break

            length = len(entry1[1])
            length2 = len(entry2[1])
            if length == length2:
                v1 = np.array([cm1.V[index] for index in entry1[1]])
                v2 = np.array([cm2.V[index] for index in entry2[1]])
                u1 = np.array([[cm1.U[m][index] for index in entry1[1]] for m in range(len(cm1.sig))])
                u2 = np.array([[cm2.U[m][index] for index in entry2[1]] for m in range(len(cm2.sig))])

                u1 = gauss(u1)
                u2 = gauss(u2)
                v1 = gauss(v1)
                v2 = gauss(v2)

                flagU = flagU and cal_shi_num(u1, u2, VorE + " u1u2")
                flagV = flagV and cal_shi_num(v1, v2, VorE + " v1v2")

        US = 0.0 if flagU else 1.0
        VS = 0.0 if flagV else 1.0

        return 0.125 * US + 0.125 * VS

def smiles_distance(smiles1:str, smiles2:str, dist_type:str = 'continuous', decimal_precision:int = 6):
    return Distance.from_smiles(smiles1, smiles2, dist_type, decimal_precision)

def double_helix_distance(seq1:str, seq2:str, dist_type:str = 'continuous', decimal_precision:int = 6):
    return Distance.from_double_helix(seq1, seq2, dist_type, decimal_precision)

def graph_distance(adjMat1:np.ndarray, adjMat2:np.ndarray, dist_type:str = 'continuous', decimal_precision:int = 6):
    return Distance.from_complex_graphs(adjMat1, adjMat2, dist_type, decimal_precision)