import numpy as np
from numba import njit
@njit
def edit_similarity(s1:np.ndarray, s2:np.ndarray):
    n1, n2 = s1.shape[0], s2.shape[0]
    dp = np.zeros(shape = (n1 + 1, n2 + 1))

    for i in range(1, n1 + 1):
        dp[i][0] = dp[i - 1][0] + int(s1[i - 1] != s2[0])
    for j in range(1, n2 + 1):
        dp[0][j] = dp[0][j - 1] + int(s1[0] != s2[j - 1])

    for i in range(1, n1 + 1):
        for j in range(1, n2 + 1):
            dp[i][j] = dp[i - 1][j - 1] + int(s1[i - 1] != s2[j - 1])
            if j < n2:
                delete = dp[i - 1][j] + int(s1[i - 1] != s2[j])
                if delete < dp[i][j]:
                    dp[i][j] = delete
            if i < n1:
                insert = dp[i][j - 1] + int(s1[i] != s2[j - 1])
                if insert < dp[i][j]:
                    dp[i][j] = insert
    return 1 - dp[n1][n2] / (n1 + n2)