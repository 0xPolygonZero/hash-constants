import itertools

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

def to_butterfly_tuples(M):
    t = [M[0], M[3], M[6], M[9]]
    u = [M[1], M[4], M[7], M[10]]
    v = [M[2], M[5], M[8], M[11]]
    return t, u, v

def from_butterfly_tuples(t, u, v):
    return [t[0], u[0], v[0],
            t[1], u[1], v[1],
            t[2], u[2], v[2],
            t[3], u[3], v[3]]

def fft_vals(a, b, c, d):
    return [a + b + c + d,
            a - b + c - d,
            a - c,
            b - d]

def fftv(M):
    t, u, v = to_butterfly_tuples(M)
    return fft_vals(*t), fft_vals(*u), fft_vals(*v)

def check_candidate(args):
    K, U, candidate_elts = args
    t, u, v = candidate_elts
    # invert the FFT butterfly
    row = vector(K, [t[0], u[0], v[0],
                     t[1], u[1], v[1],
                     t[2], u[2], v[2],
                     t[3], u[3], v[3]])
    if is_mds_circ(row):
        fft_elts = [U * r for r in candidate_elts]
        print('\n', row, ' => ', fft_elts, flush=True)
        return [(row, fft_elts)]
    return []

# U = matrix(QQ, [[1,  1,  1,  1],
#                 [1, -1,  1, -1],
#                 [1,  0, -1,  0],
#                 [0,  1,  0, -1]])
Uinv = 1/4 * matrix(QQ, [[1,  1,  2,  0],
                         [1, -1,  0,  2],
                         [1,  1, -2,  0],
                         [1, -1,  0, -2]])

def update_vec(vec, target):
    delta = vector(t - f for t, f in zip(target, fft_vals(*vec)))
    update = Uinv * vector(QQ, delta)
    if any(u.denominator() != 1 for u in update):
        return None
    new_vec = vector(QQ, vec) + update
    if any(v < 0 for v in new_vec):
        return None
    return new_vec


def find_matrix(K, max_exp=4):
    '''Find a matrix.

    The complex Fourier transform of x0, x1, x2, x3 is

    X = [x0 + x1 + x2 + x3,
         (x0 - x2) + (x1 - x3) i,
         x0 - x1 + x2 - x3,
         (x0 - x2) - (x1 - x3) i]
      = [ 1  1  1  1 ] [ x0 ]
        [ 1  i -1 -i ] [ x1 ]
        [ 1 -1  1 -1 ] [ x2 ]
        [ 1 -i -1  i ] [ x3 ]

    We want the entries of X to be small (incl. zero or one) or
    otherwise a power of two.

    '''
    U = matrix(QQ, [[1,  1,  1,  1],
                    [1, -1,  1, -1],
                    [1,  0, -1,  0],
                    [0,  1,  0, -1]])
    Uinv = U^-1
    trial_entries = [1 << i for i in range(max_exp + 1)]
    trials = itertools.product(trial_entries, repeat=4)
    candidate_tuples = []
    for trial in trials:
        v = Uinv * vector(QQ, trial)
        if all(vi != 0 for vi in v):      # and len(set(v)) >= len(v) - 1:
            # v has only non-zero entries #and no duplicates
            #print(trial, ' => ', v)
            candidate_tuples.append(v)

    print('found', len(candidate_tuples), 'candidate tuples')
    n = len(candidate_tuples)
    print('total number of candidate matrices is', n * (n-1) * (n-2))

    results = []
    tries = 0

    M = [9, 20, 4, 1, 16, 2, 22, 27, 3, 32, 1, 1]
    tuples = to_butterfly_tuples(M)
    tuples = list(tuples)

    basic_entries = [1 << i for i in range(max_exp + 1)]
    trial_entires = []
    for b in basic_entries:
        trial_entries.extend([b, -b])
    trial_entries = basic_entries

    for idx in range(0, 3):
        print(f'Index {idx}:')
        trials = list(itertools.product(trial_entries, repeat=4))
        trials.sort(key=lambda t: list(map(abs, reversed(t))))

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
                    tuples[idx] = new_t
                    break
            cnt += 1

        print(f'count: {cnt}')
    return

    # from multiprocessing import Pool, cpu_count
    # pool = Pool(cpu_count())
    # for result in pool.map(check_candidate,
    #                        zip(itertools.repeat(K),
    #                            itertools.repeat(U),
    #                            itertools.permutations(candidate_tuples, r=int(3)))):
    #     if tries % 100000 == 0:
    #        print('.', end='', flush=True)
    #     tries += 1
    #     results += result

    #return results

