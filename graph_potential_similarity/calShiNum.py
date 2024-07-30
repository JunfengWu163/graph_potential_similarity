def cal_shi_num(u1, u2, VorE):
    if u1 is None and u2 is None:
        return True

    row1 = []
    row2 = []

    for i in range(len(u1)):
        m1 = {}
        for j in range(len(u1[0])):
            if u1[i][j] in m1:
                m1[u1[i][j]] += 1
            else:
                m1[u1[i][j]] = 1
        m1 = {k: v for k, v in m1.items() if v >= 2}
        row1.append(m1)

    for i in range(len(u2)):
        m2 = {}
        for j in range(len(u2[0])):
            if u2[i][j] in m2:
                m2[u2[i][j]] += 1
            else:
                m2[u2[i][j]] = 1
        m2 = {k: v for k, v in m2.items() if v >= 2}
        row2.append(m2)

    for map1 in row1:
        iter2 = iter(row2)
        while True:
            try:
                map2 = next(iter2)
            except StopIteration:
                break

            map1_val_num = {}
            for val in map1.values():
                if val in map1_val_num:
                    map1_val_num[val] += 1
                else:
                    map1_val_num[val] = 1

            map2_val_num = {}
            for val in map2.values():
                if val in map2_val_num:
                    map2_val_num[val] += 1
                else:
                    map2_val_num[val] = 1

            if map1_val_num == map2_val_num:
                row2.remove(map2)
                break

    return not row2
