import itertools

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

def is_mds_circ(row):
    return is_mds_fast(Matrix.circulant(row))

def is_binary_power(x):
    x = int(x)
    return (x == 1) or (x & (x - 1) == 0)

def make_binary_powers(init_row):
    C = Matrix.circulant(init_row)
    k = C.base_ring()
    assert is_mds_fast(C)
    N = len(init_row)

    # ones first:
    for cols in itertools.combinations(range(N), 3):
        row = copy(init_row)
        for c in cols:
            row[c] = k(1)
        if is_mds_circ(row):
            print(f'row: {row}')
            break

    for j in range(N):
        if is_binary_power(row[j]):
            continue
        e = 1
        while True:
            row[j] = k(1 << e)
            if is_mds_circ(row):
                print(f'row: {row}')
                break
            e += 1

    return row

def random_circulant_mds(k, n):
    ntries = 0
    while True:
        random_row = [k.random_element() for _ in range(n)]
        C = Matrix.circulant(random_row)
        assert C.base_ring() == k
        ntries += 1
        if is_mds_fast(C):
            return C.row(0)
        if ntries % 100 == 0:
            print(f'Continuing after {ntries} tries')


###
### Crandall field
###

crandall_prime = 2^64 - 2415919103
crandall_field = GF(crandall_prime)

crandall_small_mds8 = vector(crandall_field, [4, 1, 2, 9, 10, 5, 1, 1])
crandall_binary_mds8 = vector(crandall_field, [4, 1, 2, 256, 16, 8, 1, 1])
# print('crandall_small_mds8 is MDS?', is_mds_circ(crandall_small_mds8))
# print('crandall_binary_mds8 is MDS?', is_mds_circ(crandall_binary_mds8))

crandall_small_mds12 = vector(crandall_field, [9, 7, 4, 1, 16, 2, 256, 128, 3, 32, 1, 1])
crandall_binary_mds12 = vector(crandall_field, [1024, 8192, 4, 1, 16, 2, 256, 128, 3, 32, 1, 1])
# print('crandall_small_mds12 is MDS?', is_mds_circ(crandall_small_mds12))
# print('crandall_binary_mds12 is MDS?', is_mds_circ(crandall_binary_mds12))


###
### Goldilocks field
###

goldilocks_prime = 2^64 - 2^32 + 1
goldilocks_field = GF(goldilocks_prime)

# Produced with make_binary_powers(random_circulant_mds(goldilocks_field, 8))
goldilocks_mds8 = vector(goldilocks_field, [1, 1, 2, 1, 8, 32, 4, 256])
# print('goldilocks_mds8 is MDS?', is_mds_circ(goldilocks_mds8))

# Produced with make_binary_powers(random_circulant_mds(goldilocks_field, 12))
goldilocks_mds12 = vector(goldilocks_field, [1, 1, 2, 1, 8, 32, 2, 256, 4096, 8, 65536, 1024])
# print('goldilocks_mds12 is MDS?', is_mds_circ(goldilocks_mds12))
