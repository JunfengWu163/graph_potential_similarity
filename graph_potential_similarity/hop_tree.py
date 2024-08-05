import numpy as np
def edit_similarity(similarity_matrix:np.ndarray):
    n1, n2 = similarity_matrix.shape[0], similarity_matrix.shape[1]
    dp = np.zeros(shape = (n1 + 1, n2 + 1))

    for i in range(1, n1 + 1):
        dp[i][0] = dp[i - 1][0] + 1 - similarity_matrix[i - 1][0]
    for j in range(1, n2 + 1):
        dp[0][j] = dp[0][j - 1] + 1 - similarity_matrix[0][j - 1]

    for i in range(1, n1 + 1):
        for j in range(1, n2 + 1):
            dp[i][j] = dp[i - 1][j - 1] + 1 - similarity_matrix[i - 1][j - 1]
            if j < n2:
                delete = dp[i - 1][j] + 1 - similarity_matrix[i - 1][j]
                if delete < dp[i][j]:
                    dp[i][j] = delete
            if i < n1:
                insert = dp[i][j - 1] + 1 - similarity_matrix[i][j - 1]
                if insert < dp[i][j]:
                    dp[i][j] = insert
    return 1 - dp[n1][n2] / (n1 + n2)

class hop_tree:
    def __init__(self, tree:list = [1]):
        self.tree = [tree[0]]
        children = []
        for i in range(1,len(tree)):
            children.append(hop_tree(tree[i]))
        children = sorted(children, reverse=True)
        for j in range(len(children)):
            self.tree.append(children[j].tree)
    
    @staticmethod
    def tree_cmp(tree1, tree2):
        if tree1[0] < tree2[0]:
            return -1
        elif tree1[0] > tree2[0]:
            return 1
        elif len(tree1) == 1:
            if len(tree2) > 1:
                return -1
            else:
                return 0
        elif len(tree2) == 1:
            return 1
        else:
            n = min(len(tree1), len(tree2))
            for i in range(1, n):
                child_cmp = hop_tree.tree_cmp(tree1[i], tree2[i])
                if child_cmp != 0:
                    return child_cmp
            if len(tree1) < len(tree2):
                return -1
            elif len(tree1) > len(tree2):
                return 1
            else:
                return 0

    def __lt__(self, other):
        return hop_tree.tree_cmp(self.tree, other.tree)
    
    @staticmethod
    def tree_sim(tree1, tree2):
        if len(tree1) == 1:
            if len(tree2) > 1:
                return 0.0
            else:
                return 1.0
        elif len(tree2) == 1:
            return 0.0
        
        # compute a sim matrix for the subtrees
        n1 = len(tree1) - 1
        n2 = len(tree2) - 1
        mat = np.zeros(shape=(n1,n2))
        for i1 in range(n1):
            for i2 in range(n2):
                mat[i1,i2] = hop_tree.tree_sim(tree1[i1 + 1], tree2[i2 + 1])
        return edit_similarity(mat)
    
    def get_similarity(self,other):
        return hop_tree.tree_sim(self.tree,other.tree)
        