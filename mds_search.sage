import itertools

p = 2^64 - 2415919103  # Crandall prime
k = GF(p)
crandall_MDS = matrix(k, 8, 8,
           [[16140901062381928449, 2635249153041947502, 3074457345215605419, 11068046442776179508, 13835058053470224385, 6148914690431210838, 9223372035646816257, 1],
            [2049638230143736946, 16140901062381928449, 2635249153041947502, 3074457345215605419, 11068046442776179508, 13835058053470224385, 6148914690431210838, 9223372035646816257],
            [5534023221388089754, 2049638230143736946, 16140901062381928449, 2635249153041947502, 3074457345215605419, 11068046442776179508, 13835058053470224385, 6148914690431210838],
            [16769767337539665921, 5534023221388089754, 2049638230143736946, 16140901062381928449, 2635249153041947502, 3074457345215605419, 11068046442776179508, 13835058053470224385],
            [10760600708254618966, 16769767337539665921, 5534023221388089754, 2049638230143736946, 16140901062381928449, 2635249153041947502, 3074457345215605419, 11068046442776179508],
            [5675921252705733081, 10760600708254618966, 16769767337539665921, 5534023221388089754, 2049638230143736946, 16140901062381928449, 2635249153041947502, 3074457345215605419],
            [1317624576520973751, 5675921252705733081, 10760600708254618966, 16769767337539665921, 5534023221388089754, 2049638230143736946, 16140901062381928449, 2635249153041947502],
            [15987178195121148178, 1317624576520973751, 5675921252705733081, 10760600708254618966, 16769767337539665921, 5534023221388089754, 2049638230143736946, 16140901062381928449]])

def is_mds_slow(M, noisy=False):
    Mbar = matrix.block([[1], [M]]) # Put the identity matrix on top of M
    indices = itertools.combinations(range(Mbar.nrows()), M.nrows())

    for row_indices in indices:
        if not Mbar.delete_rows(row_indices).is_invertible():
            if noisy: print(f'Failed at {row_indices}')
            return False
    return True

def is_mds_fast(A, noisy=False):
    # 1-minors are just the elements themselves
    if any(any(r == 0 for r in row) for row in A):
        if noisy: print('FAILURE: matrix has zero entry')
        return False

    N = A.nrows()
    assert A.is_square() and N >= 2

    det_cache = A

    # Calculate all the nxn minors of A:
    for n in range(2, N+1):
        new_det_cache = dict()
        for rows in itertools.combinations(range(N), n):
            for cols in itertools.combinations(range(N), n):
                i, *rs = rows

                # Laplace expansion along row i
                det = 0
                for j in range(n):
                    # pick out c = column j; the remaining columns are in cs
                    c = cols[j]
                    cs = cols[:j] + cols[j+1:]

                    # Look up the determinant from the previous iteration
                    # and multiply by -1 if j is odd
                    cofactor = det_cache[(*rs, *cs)]
                    if j % 2 == 1:
                        cofactor = -cofactor

                    # update the determinant with the j-th term
                    det += A[i, c] * cofactor

                if det == 0:
                    if noisy: print(f'FAILURE on {n}-minor: rows={rows}, cols={cols}')
                    return False
                new_det_cache[(*rows, *cols)] = det
        if noisy: print(f'matrix has no zero {n}-minors')
        det_cache = new_det_cache
    return True


def is_mds_circulant(row):

    # 2-minors
    new_det_cache = dict()
    for rows in itertools.combinations(range(N), 2):
        for cols in itertools.combinations(range(N), 2):
            r1, r2 = rows
            c1, c2 = cols

            a = (c1 - r1) % N # (r1, c1)
            b = (c2 - r1) % N # (r1, c2)
            c = (c1 - r2) % N # (r2, c1)
            d = (c2 - r2) % N # (r2, c2)

            det = C.matrix_from_rows_and_columns(rows, cols).det()

            assert det == row[a]*row[d] - row[b]*row[c]

            if det == 0:
                return f'FAILURE on 2-minor: rows={rows}, cols={cols}'
            new_det_cache[(a, b, c, d)] = det

    print('matrix has no zero 2-minors')

    # 3-minors
    new_det_cache = dict()
    for rows in itertools.combinations(range(N), 3):
        for cols in itertools.combinations(range(N), 3):
            r1, r2, r3 = rows
            c1, c2, c3 = cols

            a = (c1 - r1) % N # (r1, c1)
            b = (c2 - r1) % N # (r1, c2)
            c = (c3 - r1) % N # (r1, c3)
            d = (c1 - r2) % N # (r2, c1)
            e = (c2 - r2) % N # (r2, c2)
            f = (c3 - r2) % N # (r2, c3)
            g = (c1 - r3) % N # (r3, c1)
            h = (c2 - r3) % N # (r3, c2)
            i = (c3 - r3) % N # (r3, c3)

            det = C.matrix_from_rows_and_columns(rows, cols).det()

            # Laplace expansion along row r1
            det = C[r1, c1] * det_cache[(r2, r3, c2, c3)] \
                - C[r1, c2] * det_cache[(r2, r3, c1, c3)] \
                + C[r1, c3] * det_cache[(r2, r3, c1, c2)]

            if det == 0:
                return False
            new_det_cache[(*rows, *cols)] = det


def reduce_mds(MM):
    assert is_mds(MM)

    M = matrix(MM)
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            N = matrix(M)
            while N[i,j] > 16:
                N[i,j] = N[i,j] >> 4
                if is_mds(N):
                    M = matrix(N)
                else:
                    break
            print(f'After {i}, {j}:\n{M}')
    return M


def reduce_circulant_mds(init_row, rate_bits=16):
    C = Matrix.circulant(init_row)
    k = C.base_ring()
    assert is_mds_fast(C)
    N = len(init_row)

    # ones first:
    for cols in itertools.combinations(range(N), 4):
        row = copy(init_row)
        for c in cols:
            row[c] = k(1)
        if is_mds_fast(Matrix.circulant(row)):
            print(f'row: {row}')

    row = copy(init_row)
    for j in range(len(row)):
        while row[j] > (1 << rate_bits):
            rj = row[j]
            row[j] = row[j] >> rate_bits
            if not is_mds_fast(Matrix.circulant(row)):
                row[j] = rj  # restore row
                break
        print(f'After {j}:\n{row}')
    return row


row0 = [4, 1, 2, 9, 10, 5, 1, 1]
row0 = [4, 1, 2, 256, 16, 8, 1, 1]
C = matrix.circulant(row0)
C = C.change_ring(k)

#print('C is MDS?', is_mds(C))

# row0 = [4, 1, 2, 256, 16, 8, 1, 1, 12345, 6789, 1111, 4321]
# D = matrix.circulant(row0)
# D = D.change_ring(k)
# print('D is MDS?', is_mds(D))
# #print('reduction:', reduce_mds2(D))


def do_thing(mat_rows, submat_rows):
    N = mat_rows
    n = submat_rows
    assert n < N

    # A is a tuple of indeterminates a_0, ..., a_{N-1}
    A = var(','.join(map(lambda i: f'a{i}', range(N))))
    #A = primes_first_n(N)
    C = Matrix.circulant(A)

    M = [C.matrix_from_rows_and_columns(rows, cols)
        for rows in itertools.combinations(range(N), n)
        for cols in itertools.combinations(range(N), n)]
    print(f'total submatrices:        {len(M):5}')

    for m in M:
        m.set_immutable()
    MM = set(M)
    t = Integer(len(MM))
    print(f'unique submatrices:       {t:5}   {t.factor()}')

    D = set(m.det() for m in M)
    t = Integer(len(D))
    print(f'unique minors:            {t:5}   {t.factor()}')

    E = set()
    for d in D:
        if d not in E and -d not in E:
            E.add(d)
    t = Integer(len(E))
    print(f'unique minors up to sign: {t:5}   {t.factor()}')



rescue_MDS = matrix(k, 12, 12,
    [[10760600708254618966, 16769767337539665921, 5534023221388089754, 2049638230143736946, 16140901062381928449, 2635249153041947502, 3074457345215605419, 11068046442776179508, 13835058053470224385, 6148914690431210838, 9223372035646816257, 1,],
    [5675921252705733081, 10760600708254618966, 16769767337539665921, 5534023221388089754, 2049638230143736946, 16140901062381928449, 2635249153041947502, 3074457345215605419, 11068046442776179508, 13835058053470224385, 6148914690431210838, 9223372035646816257,],
    [1317624576520973751, 5675921252705733081, 10760600708254618966, 16769767337539665921, 5534023221388089754, 2049638230143736946, 16140901062381928449, 2635249153041947502, 3074457345215605419, 11068046442776179508, 13835058053470224385, 6148914690431210838,],
    [15987178195121148178, 1317624576520973751, 5675921252705733081, 10760600708254618966, 16769767337539665921, 5534023221388089754, 2049638230143736946, 16140901062381928449, 2635249153041947502, 3074457345215605419, 11068046442776179508, 13835058053470224385,],
    [17293822566837780481, 15987178195121148178, 1317624576520973751, 5675921252705733081, 10760600708254618966, 16769767337539665921, 5534023221388089754, 2049638230143736946, 16140901062381928449, 2635249153041947502, 3074457345215605419, 11068046442776179508,],
    [3255307777287111620, 17293822566837780481, 15987178195121148178, 1317624576520973751, 5675921252705733081, 10760600708254618966, 16769767337539665921, 5534023221388089754, 2049638230143736946, 16140901062381928449, 2635249153041947502, 3074457345215605419,],
    [1024819115071868473, 3255307777287111620, 17293822566837780481, 15987178195121148178, 1317624576520973751, 5675921252705733081, 10760600708254618966, 16769767337539665921, 5534023221388089754, 2049638230143736946, 16140901062381928449, 2635249153041947502,],
    [9708812669101911849, 1024819115071868473, 3255307777287111620, 17293822566837780481, 15987178195121148178, 1317624576520973751, 5675921252705733081, 10760600708254618966, 16769767337539665921, 5534023221388089754, 2049638230143736946, 16140901062381928449,],
    [2767011610694044877, 9708812669101911849, 1024819115071868473, 3255307777287111620, 17293822566837780481, 15987178195121148178, 1317624576520973751, 5675921252705733081, 10760600708254618966, 16769767337539665921, 5534023221388089754, 2049638230143736946,],
    [878416384347315834, 2767011610694044877, 9708812669101911849, 1024819115071868473, 3255307777287111620, 17293822566837780481, 15987178195121148178, 1317624576520973751, 5675921252705733081, 10760600708254618966, 16769767337539665921, 5534023221388089754,],
    [17608255704416649217, 878416384347315834, 2767011610694044877, 9708812669101911849, 1024819115071868473, 3255307777287111620, 17293822566837780481, 15987178195121148178, 1317624576520973751, 5675921252705733081, 10760600708254618966, 16769767337539665921,],
    [15238614667590392076, 17608255704416649217, 878416384347315834, 2767011610694044877, 9708812669101911849, 1024819115071868473, 3255307777287111620, 17293822566837780481, 15987178195121148178, 1317624576520973751, 5675921252705733081, 10760600708254618966]])

#is_mds_fast(M, noisy=True)

row12 = [9, 7, 4, 1, 16, 2, 256, 128, 3, 32, 1, 1]
row12 = [1024, 8192, 4, 1, 16, 2, 256, 128, 3, 32, 1, 1]

C = Matrix.circulant(row12)
C = C.change_ring(k)

