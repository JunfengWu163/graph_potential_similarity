import numpy as np
from .hop_tree import hop_tree, edit_similarity
class directed_potential_graph:
    def __init__(self, m:np.array):
        self.n:int = m.shape[0]
        self.edges = {}
        for i in range(self.n):
            for j in range(self.n):
                if m[i,j] > 0:
                    if not i in self.edges:
                        self.edges[i] = set()
                    self.edges[i].add(j)
        self.hop_trees = self.get_hop_trees()

    def get_subtree(self, i:int, visited:set = set()):
        assert i < self.n and i >= 0
        if not i in self.edges:
            return [1]
        total_weight = 1
        subtrees = []
        new_visited = visited.union(set(self.edges[i]))
        new_visited.add(i)
        for j in self.edges[i]:
            if not j in visited:
                subtree = self.get_subtree(j, new_visited)
                total_weight += subtree[0]
                subtrees.append(hop_tree(subtree))
        tree = [total_weight]
        # do not sort subtrees in get_subtree. sort them at get_hop_tree
        for j in range(len(subtrees)):
            tree.append(subtrees[j].tree)
        return tree
    
    def get_hop_tree(self, i:int):
        return hop_tree(self.get_subtree(i))
    
    def get_hop_trees(self):
        hop_trees = []
        for i in range(self.n):
            hop_trees.append(self.get_hop_tree(i))
        return sorted(hop_trees, reverse=True)
    
    def get_similarity(self, other):
        max_level = 3
        sim = 0.0
        n1 = len(self.hop_trees)
        n2 = len(other.hop_trees)
        for level in range(max_level):
            s1 = []
            for i in range(n1):
                s1 += self.hop_trees[i].to_seq(level)
            s2 = []
            for i in range(n2):
                s2 += other.hop_trees[i].to_seq(level)
            mat = np.zeros(shape=(len(s1),len(s2)))
            # print(f's1 = {s1}')
            # print(f's2 = {s2}')
            for i1 in range(len(s1)):
                for i2 in range(len(s2)):
                    mat[i1,i2] = min(s1[i1],s2[i2]) / max(s1[i1],s2[i2])
            sim += edit_similarity(mat) / max_level
        return sim