###
### This file contains utilities for finding MDS matrices whose entries
### of small powers of 2.
###

import itertools

def is_mds_slow(M, noisy=False):
    '''Return True iff M is an MDS matrix.

    Algorithm is the naive one: If M is nxn, create a block matrix by
    putting an nxn identity matrix on top of M, then check that
    deleting any n rows of that block matrix is still invertible.
    [Ref: 'https://en.wikipedia.org/wiki/MDS_matrix']

    NB: This function is included mainly as a way to double check
    results, since `is_mds_fast` is complicated and a bug could
    potentially be hiding in there somewhere.
    '''
    Mbar = matrix.block([[1], [M]]) # Put the identity matrix on top of M
    indices = itertools.combinations(range(Mbar.nrows()), M.nrows())

    for row_indices in indices:
        if not Mbar.delete_rows(row_indices).is_invertible():
            if noisy: print(f'Failed at {row_indices}')
            return False
    return True

def is_mds_fast(A, noisy=False):
    '''Return True iff A is an MDS matrix.

    This function uses a Faster algorithm to avoid lots of
    recomputation: A matrix is MDS if all its submatrices are
    invertible, i.e. all its minors are non-zero; so we use the
    Laplace expansion of the determinant to calculate the (m+1)x(m+1)
    minors in terms of the (already computed) mxm minors [Ref:
    'https://en.wikipedia.org/wiki/Laplace_expansion']. There is
    probably a smarter way to do this, but this algorithm basically
    just puts the mxm minors in a dictionary and looks them up when
    calculating the (m+1)x(m+1) minors.
    '''

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
    '''Return the circulant matrix whose first row is 'row'.'''
    return is_mds_fast(Matrix.circulant(row))

def is_binary_power(x):
    '''Return true iff x = 2^e for some integer e >= 0.'''
    x = int(x)
    return (x == 1) or (x & (x - 1) == 0)

def next_binary_power(x, inclusive=True):
    '''Return the smallest binary power larger than x.

    If inclusive is True, return x if x is already a binary power;
    otherwise, return the smallest binary power strictly larger than x.
    '''
    if inclusive and is_binary_power(x):
        return x
    return 1 << int(x).bit_length()

def make_binary_powers(init_row, noisy=True):
    '''Use given circulant MDS matrix to generate a circulant MDS matrix
    with entries that are powers of two.

    Given the initial row of a circulant matrix, we (i) replace three
    of entries with 1 (three is the maximum number of 1's we can hope
    for) and (ii) replace each of the other entries with the smallest
    power of two that maintains the MDS property.
    '''
    C = Matrix.circulant(init_row)
    field = C.base_ring()
    assert is_mds_fast(C)
    N = len(init_row)

    # ones first:
    for cols in itertools.combinations(range(N), 3):
        row = copy(init_row)
        for c in cols:
            row[c] = field(1)
        if is_mds_circ(row):
            if noisy: print(f'row: {row}')
            break

    # for each non-binary power, replace it with successively larger
    # binary powers until the MDS property holds.
    for j in range(N):
        if is_binary_power(row[j]):
            continue
        e = 1  # the exponent
        while True:
            row[j] = field(1 << e)
            if is_mds_circ(row):
                if noisy: print(f'row: {row}')
                break
            e += 1

    return row

def random_circulant_mds(k, n, noisy=True):
    '''Return a random circulant MDS matrix.'''
    ntries = 0
    while True:
        random_row = vector(k, [k.random_element() for _ in range(n)])
        ntries += 1
        if is_mds_circ(random_row):
            return random_row
        if noisy and ntries % 100 == 0:
            print(f'Continuing after {ntries} tries')


###
### Crandall field
###

crandall_prime = 2^64 - 2415919103
crandall_field = GF(crandall_prime)

## These four were 'manually crafted'; not necessary better, just from
## an earlier version.
crandall_small_manual_mds8 = vector(crandall_field, [4, 1, 2, 9, 10, 5, 1, 1])
crandall_binary_manual_mds8 = vector(crandall_field, [4, 1, 2, 256, 16, 8, 1, 1])
#print('crandall_small_manual_mds8 is MDS?', is_mds_circ(crandall_small_manual_mds8))
#print('crandall_binary_manual_mds8 is MDS?', is_mds_circ(crandall_binary_manual_mds8))

crandall_small_manual_mds12 = vector(crandall_field, [9, 7, 4, 1, 16, 2, 256, 128, 3, 32, 1, 1])
crandall_binary_manual_mds12 = vector(crandall_field, [1024, 8192, 4, 1, 16, 2, 256, 128, 3, 32, 1, 1])
#print('crandall_small_mds12 is MDS?', is_mds_circ(crandall_small_manual_mds12))
#print('crandall_binary_mds12 is MDS?', is_mds_circ(crandall_binary_manual_mds12))

## Produced with:
# make_binary_powers(random_circulant_mds(crandall_field, 8)) produces
crandall_binary_mds8 = vector(crandall_field, [1, 1, 2, 1, 8, 32, 4, 256])
crandall_binary_mds12 = vector(crandall_field, [1, 1, 2, 1, 8, 32, 2, 256, 4096, 8, 65536, 1024])
#print('crandall_small_mds8 is MDS?', is_mds_circ(crandall_small_mds8))
#print('crandall_binary_mds8 is MDS?', is_mds_circ(crandall_binary_mds8))


###
### Goldilocks field
###

goldilocks_prime = 2^64 - 2^32 + 1
goldilocks_field = GF(goldilocks_prime)

# Produced with make_binary_powers(random_circulant_mds(goldilocks_field, 8))
goldilocks_mds8 = vector(goldilocks_field, [1, 1, 2, 1, 8, 32, 4, 256])
#print('goldilocks_mds8 is MDS?', is_mds_circ(goldilocks_mds8))

# Produced with make_binary_powers(random_circulant_mds(goldilocks_field, 12))
goldilocks_mds12 = vector(goldilocks_field, [1, 1, 2, 1, 8, 32, 2, 256, 4096, 8, 65536, 1024])
#print('goldilocks_mds12 is MDS?', is_mds_circ(goldilocks_mds12))
