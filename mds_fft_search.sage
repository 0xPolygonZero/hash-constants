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
                    [0,  1,  0, -1],
                    [1,  0, -1,  0]])
    Uinv = U^-1
    trial_entries = [0] + [1 << i for i in range(max_exp + 1)]
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

    # for candidate_elts in itertools.combinations(candidate_tuples, 3):
    #     t, u, v = candidate_elts
    #     # invert the FFT butterfly
    #     row = vector(K,
    #                  [t[0], u[0], v[0],
    #                   t[1], u[1], v[1],
    #                   t[2], u[2], v[2],
    #                   t[3], u[3], v[3]])
    #     if is_mds_circ(row):
    #         fft_elts = [U * r for r in candidate_elts]
    #         print('\n', row, ' => ', fft_elts, flush=True)
    #         results.append((row, fft_elts))
    #     elif tries % 100000 == 0:
    #         print(tries, flush=True)
    #     tries += 1

    from multiprocessing import Pool, cpu_count
    pool = Pool(cpu_count())
    for result in pool.map(check_candidate,
                           zip(itertools.repeat(K),
                               itertools.repeat(U),
                               itertools.permutations(candidate_tuples, r=int(3)))):
        if tries % 100000 == 0:
           print('.', end='', flush=True)
        tries += 1
        results += result

    return results

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


def find_complex_coeff(K, tuples, idx, f_idx):
    t_orig = tuples[idx]
    print('t =', t_orig)
    F = fft_vals(*t_orig)
    assert f_idx == 2 or f_idx == 3
    f_val = F[f_idx]
    print('F =', F, '(f_val is', f_val, ')')
    t = t_orig
    tuples_orig = list(list(t) for t in tuples)

    def check_target(tuples, idx, target, F, f_idx):
        f_val = F[f_idx]
        i, j = (0, 2) if f_idx == 2 else (1, 3)
        print('trying for', target)
        t = tuples[idx]
        print(f'{t[i]}, {t[j]} -> ', end='')
        ti = t[i] + (target - f_val) // 2
        tj = t[j] - (target - f_val) // 2
        new_t = list(t)
        new_t[i] = ti
        new_t[j] = tj
        print(f'{ti}, {tj} ... ', end='', flush=True)
        print(f'(fft: {fft_vals(*new_t)}):', end=' ', flush=True)
        new_tuples = list(tuples)
        new_tuples[idx] = new_t
        ok = is_mds_circ(vector(K, from_butterfly_tuples(*new_tuples)))
        print('OK' if ok else 'FAIL')
        return ok, new_tuples

    def flip_f_val_parity(tuples, idx):
        # Switch parity of t[0] - t[2]; this modifies f2 but not f1
        t_orig = list(tuples[idx])
        t = tuples[idx]
        t[0] += 1
        t[1] += 1
        t[3] -= 2
        while True:
            print('.', end='', flush=True)
            tuples[idx] = t
            ok = is_mds_circ(vector(K, from_butterfly_tuples(*tuples)))
            if ok:
                print(f'{t_orig} -> {t} (fft: {fft_vals(*t)})')
                return tuples
            t[0] += 2
            t[1] += 2
            t[3] -= 4

    # First try to make f_val +/- 1; for this f_val must be odd
    if f_val % 2 == 0:
        # f_val is even; make it odd
        print('find equivalent t with odd f_val')
        tuples = flip_f_val_parity(tuples, idx)
        F = fft_vals(*tuples[idx])

    ok, new_tuples = check_target(tuples, idx, 1, F, f_idx)
    if ok:
        return new_tuples
    ok, new_tuples = check_target(tuples, idx, -1, F, f_idx)
    if ok:
        return new_tuples

    # Can't make f_val +/- 1, so we make it a power of two; for this
    # f_val must be even.
    print('settle for power of two')
    if f_val % 2 == 1:
        # f_val is odd; make it even
        print('find equivalent t with even f_val')
        tuples = flip_f_val_parity(tuples_orig, idx)
    else:
        tuples = tuples_orig

    t = tuples[idx]
    F = fft_vals(*tuples[idx])

    target = 2
    while True:
        ok, new_tuples = check_target(tuples, idx, target, F, f_idx)
        if ok:
            break
        target *= 2
    return new_tuples

def find_first_val(t):
    pass

def matrix_search_t(K, M, idx):
    tuples = to_butterfly_tuples(M)
    t = tuples[idx]

    print(f'M = {M}  =>  {tuples}')
    print('t =', t)
    F = fft_vals(*t)
    print('F =', F)
    f1, f2, f3, f4 = F

    print('\n=== FFT value 1 ===')
    print(f'{t} -> ', end='')
    if idx == 0:
        t[2] = 22
    elif idx == 1:
        t[2] = 104
    elif idx == 2:
        t[3] = 7
    else:
        assert False
    print(f'{t} ... ', end='', flush=True)
    tuples = list(tuples)
    tuples[idx] = t
    ok = is_mds_circ(vector(K, from_butterfly_tuples(*tuples)))
    print('OK' if ok else 'fail')
    if not ok:
        return

    #target1 = 64
    # while True:
    #     #target1 = next_binary_power(f1)
    #     print(f'{t[2]} -> ', end='')
    #     #t[0] += target1 - f1
    #     t[2] += target1 - f1
    #     print(f'{t[2]} ... ', end='', flush=True)
    #     f1 = target1
    #     ok = is_mds_circ(vector(K, from_butterfly_tuples(t, u, v)))
    #     print('OK' if ok else 'fail')
    #     if ok: break

    print('t =', t)
    F = fft_vals(*t)
    print('F =', F)
    f1, f2, f3, f4 = F

    print('\n=== FFT value 2 ===')

    if not is_binary_power(f2):

    print('t =', t)
    F = fft_vals(*t)
    print('F =', F)
    f1, f2, f3, f4 = F

    ###
    ## This is to force a zero in the second FFT value; doesn't seem
    ## to find anything (probably because forcing a value to be zero
    ## forces a linear relation between the values, which goes
    ## contrary to the MDS condition). Might be worth trying to spread
    ## the increment across three or four values, rather than just
    ## two?
    ###
    # increment = f1
    # while True:
    #     print(f'incr: {increment}, f2: {f2}')
    #     print(f'{t[0]}, {t[1]} -> ', end='')
    #     t0 = t[0] + (increment - f2) // 2
    #     t1 = t[1] + (increment + f2) // 2
    #     t_bar = (t0, t1, t[2], t[3])
    #     print(f'{t0}, {t1} ... ', end='', flush=True)
    #     print('(fft:', fft_vals(*t_bar), ')')
    #     ok = is_mds_circ(vector(K, from_butterfly_tuples(t_bar, u, v)))
    #     print('OK' if ok else 'fail')
    #     if ok:
    #         t[0] = t0
    #         t[1] = t1
    #         break
    #     increment = 2*increment + f1

    print('\n=== FFT value 3 ===')

    tuples = find_complex_coeff(K, tuples, idx, 2)
    t = tuples[idx]

    print('t =', t)
    F = fft_vals(*t)
    print('F =', F)
    f1, f2, f3, f4 = F

    print('\n=== FFT value 4 ===')
    tuples = find_complex_coeff(K, tuples, idx, 3)
    t = tuples[idx]

    print('t =', t)
    F = fft_vals(*t)
    print('F =', F)
    f1, f2, f3, f4 = F

    return from_butterfly_tuples(*tuples)

def matrix_search(K):
    # M is the starting point, a small MDS matrix over Goldilocks
    M = [9, 7, 4, 1, 16, 2, 256, 128, 3, 32, 1, 1]

    print('\n=== FIRST TUPLE ===')
    M = matrix_search_t(K, M, 0)
    print('\n=== SECOND TUPLE ===')
    M = matrix_search_t(K, M, 1)
    print('\n=== THIRD TUPLE ===')
    M = matrix_search_t(K, M, 2)

    return M
