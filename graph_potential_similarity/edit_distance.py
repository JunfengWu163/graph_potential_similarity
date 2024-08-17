import numpy as np
from numba import njit
@njit
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