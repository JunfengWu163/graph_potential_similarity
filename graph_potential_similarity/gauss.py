import numpy as np

EXP = 1e-6
PRECISION = 1000000

def gauss(arr):
    if len(arr) > len(arr[0]):
        return None

    result = rref(arr)
    R = result['R']
    J = result['J']

    change_precision(R, PRECISION)

    J_sy = [num for num in range(len(R[0])) if num not in J]
    R_del = np.zeros((len(R), len(J_sy)))

    for i in range(len(R)):
        for j in range(len(J_sy)):
            R_del[i][j] = R[i][J_sy[j]]

    return R_del

def rref(arr):
    row, col = len(arr), len(arr[0])
    j_save = []
    num_j_save = 0

    i, j = 0, 0
    while i < row and j < col:
        max_row = max(range(i, row), key=lambda r: abs(arr[r][j]))
        if abs(arr[max_row][j]) > EXP:
            if max_row != i:
                arr[i], arr[max_row] = arr[max_row], arr[i]
            j_save.append(j)
            num_j_save += 1

            arr[i] = [x / arr[i][j] for x in arr[i]]

            for r in range(i + 1, row):
                k = arr[r][j]
                arr[r] = [arr[r][c] - arr[i][c] * k for c in range(col)]

            i += 1
        j += 1

    for r in range(num_j_save - 1, -1, -1):
        c = j_save[r]
        for r2 in range(r - 1, -1, -1):
            k = arr[r2][c]
            arr[r2] = [arr[r2][c2] - k * arr[r][c2] for c2 in range(col)]

    J = j_save[:num_j_save]
    return {'R': arr, 'J': J}

def change_precision(arr, precision):
    for i in range(len(arr)):
        for j in range(len(arr[0])):
            arr[i][j] = round(arr[i][j] * precision) / precision

def print_arr2(arr):
    for row in arr:
        print(" ".join(f"{x:25.13f}" for x in row))
    print()

def print_arr1wei(arr):
    print(" ".join(map(str, arr)))
    print()


