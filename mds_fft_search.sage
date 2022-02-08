import itertools

# For check_minpoly_condition, algorithm_1, algorithm_2, algorithm_3
#
# Annoyingly, it seems we have to comment out the `if __name__ ==
# '__main__'` clause and `def main(row):` function at the bottom of
# `mds_security.sage` for this to work; not sure why.
load('mds_security.sage')

goldilocks_prime = 2^64 - 2^32 + 1
goldilocks_field = GF(goldilocks_prime)

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
    '''Return True if the circulant matrix whose first row is 'row' is MDS.'''
    return is_mds_fast(Matrix.circulant(row))

def to_butterfly_tuples(M):
    '''Given a 12-vector M, return the three 4-vectors from its radix-3
    Fourier transform.'''
    t = [M[0], M[3], M[6], M[9]]
    u = [M[1], M[4], M[7], M[10]]
    v = [M[2], M[5], M[8], M[11]]
    return t, u, v

def from_butterfly_tuples(t, u, v):
    return [t[0], u[0], v[0],
            t[1], u[1], v[1],
            t[2], u[2], v[2],
            t[3], u[3], v[3]]

def FT_vals(a, b, c, d):
    '''Given the four elements of a length 4 vector, return the values
    that appear in its complex Fourier transform.

    The complex Fourier transform of a 4-vector (a, b, c, d) is

    [a+b+c+d, (a-c) + (b-d)*i, a-b+c-d, (a-c) - (b-d)*i]

    where 'i' is sqrt(-1). This function returns the first and third
    elements, as well as the real and imaginary parts of the second
    and fourth elements.
    '''
    return [a + b + c + d,
            a - b + c - d,
            a - c,
            b - d]

def FT(M):
    '''Given a 12-vector M, return the Fourier transform values (as
    defined in FT_vals) of the three "butterfly tuples" of M.'''
    t, u, v = to_butterfly_tuples(M)
    return FT_vals(*t), FT_vals(*u), FT_vals(*v)

# U = matrix(QQ, [[1,  1,  1,  1],
#                 [1, -1,  1, -1],
#                 [1,  0, -1,  0],
#                 [0,  1,  0, -1]])
Uinv = 1/4 * matrix(QQ, [[1,  1,  2,  0],
                         [1, -1,  0,  2],
                         [1,  1, -2,  0],
                         [1, -1,  0, -2]])

def update_vec(vec, target):
    delta = vector(t - f for t, f in zip(target, FT_vals(*vec)))
    update = Uinv * vector(QQ, delta)
    if any(u.denominator() != 1 for u in update):
        return None
    new_vec = vector(QQ, vec) + update
    if any(v < 0 for v in new_vec):
        return None
    return new_vec

class DoCheck:
    def __init__(self, K, tuples, idx):
        self.K = K
        self.tuples = tuples
        self.idx = idx

    def __call__(self, target):
        t = self.tuples[self.idx]
        new_t = update_vec(t, target)
        if new_t is None:
            return None
        new_tuples = list(self.tuples)
        new_tuples[self.idx] = new_t
        ok = is_mds_circ(vector(self.K, from_butterfly_tuples(*new_tuples)))
        return new_t if ok else None

def find_tuple_at_idx_par(K, tuples, idx, trials):
    from multiprocessing import Pool, cpu_count
    pool = Pool(cpu_count())
    check = DoCheck(K, tuples, idx)
    return pool.imap(check, trials)


def find_tuple_at_idx(K, tuples, idx, trials):
    cnt = 0
    t = tuples[idx]
    for target in trials:
        new_t = update_vec(t, target)
        if new_t is not None:
            new_tuples = list(tuples)
            new_tuples[idx] = new_t
            ok = is_mds_circ(vector(K, from_butterfly_tuples(*new_tuples)))
            if ok:
                print(f'  {t} -> {new_t} gives {target}')
                # new tuples[idx] tuple
                yield new_t
        cnt += 1
    print(f'count: {cnt}')

def trial_list(max_exp):
    basic_entries = [1 << i for i in range(max_exp + 1)]
    trial_entries = []
    for b in basic_entries:
        trial_entries.extend([b, -b])

    # TODO: Document choice of sort order
    trials = list(itertools.product(trial_entries, repeat=4))
    trials.sort(key=lambda t: list(map(abs, reversed(t))))
    return trials

def find_matrix_tree(K, max_exp=4, M=None):
    if M == None:
        M = [9, 20, 4, 1, 16, 2, 22, 27, 3, 32, 1, 1]
    tuples = list(to_butterfly_tuples(M))
    trials = trial_list(max_exp)
    print(f'Index 0:')
    for t0 in find_tuple_at_idx(K, tuples, 0, trials):
        print(f' * Index 1:')
        for t1 in find_tuple_at_idx(K, (t0, tuples[1], tuples[2]), 1, trials):
            print(f' ** Index 2:')
            last_tups = find_tuple_at_idx_par(K, (t0, t1, tuples[2]), 2, trials)
            for t2 in last_tups:
                if t2 is None:
                    continue
                res = from_butterfly_tuples(t0, t1, t2)
                C = matrix.circulant(res)
                C.change_ring(K)
                mp = check_minpoly_condition(C)
                sec1, _ = algorithm_1(C)
                sec2, _ = algorithm_2(C)
                sec3, _ = algorithm_3(C)
                if all((sec1, sec2, sec3)):
                    print(f'    found: {res} -> {FT(res)}')
                    print(f'    SAFE!! (minpoly: {mp})')
                    if mp:
                        last_tups.terminate()
                        return res
