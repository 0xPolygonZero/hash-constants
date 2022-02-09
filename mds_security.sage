# TODO: Work out what the role of this global is.
s = 1

def check_minpoly_condition(M):
    t = M.ncols()
    max_period = 2*t
    M_pow_i = M
    for i in range(1, max_period + 1):
        minpol = M_pow_i.minimal_polynomial()
        if minpol.degree() != t or not minpol.is_irreducible():
            return [False, i]
        M_pow_i = M * M_pow_i
    return True

def isAllInvertible(M):
    t = M.nrows()
    # Test all square submatrices for invertibility
    all_invertible = True
    for i in range(2, t):
        choices_i = Combinations(range(0, t), i)
        for m in range(0, choices_i.cardinality()):
            for n in range(0, choices_i.cardinality()):
                M_sub = M[choices_i[m], choices_i[n]]
                is_inv = M_sub.is_invertible()
                all_invertible = all_invertible and is_inv
    return all_invertible

def subspace_times_matrix(subspace, M):
    F = M.base_ring()
    t = M.ncols()
    V = VectorSpace(F, t)
    subspace_basis = subspace.basis()
    new_basis = []
    for vec in subspace_basis:
        new_basis.append(M * vec)
    new_subspace = V.subspace(new_basis)
    return new_subspace

# Returns True if the matrix is considered secure, False otherwise
def algorithm_1(M):
    F = M.base_ring()
    t = M.ncols()
    R.<x> = PolynomialRing(F)
    decomposition = M.minimal_polynomial().squarefree_decomposition()

    # Get A's
    V = VectorSpace(F, t)
    unit_vector_space = V.subspace(V.basis()[s:])
    A_list = []
    basis_vectors = []
    for i in range(0, len(list(decomposition))):
        poly = list(decomposition)[i][0]
        exponent = list(decomposition)[i][1]
        A_i = R(poly^exponent)(M).right_kernel()
        A_list.append(A_i)

    basis_vectors = []
    for A_i in A_list:
        X_i = A_i.intersection(unit_vector_space)
        while X_i.dimension() > 0:
            X_i_new = X_i.intersection(subspace_times_matrix(X_i, M))
            if X_i == X_i_new:
                break
            X_i = X_i_new
        basis_vectors += X_i.basis()
    
    P_full_space = V.subspace(basis_vectors)
    if P_full_space.dimension() > 0:
        return [False, P_full_space]

    return [True, 0]

# Returns True if the matrix is considered secure, False otherwise
def algorithm_2(M):
    F = M.base_ring()
    t = M.ncols()
    V = VectorSpace(F, t)
    trail = [None, None]
    test_next = False
    I = list(range(0, s))
    I_powerset = list(sage.misc.misc.powerset(I))[1:]
    for I_s in I_powerset:
        test_next = False
        new_basis = []
        for l in I_s:
            new_basis.append(V.basis()[l])
        IS = V.subspace(new_basis)
        for i in range(s, t):
            new_basis.append(V.basis()[i])
        full_iota_space = V.subspace(new_basis)
        for l in I_s:
            v = V.basis()[l]
            while True:
                delta = IS.dimension()
                v = M * v
                IS = V.subspace(IS.basis() + [v])
                if IS.dimension() == t or IS.intersection(full_iota_space) != IS:
                    test_next = True
                    break
                if IS.dimension() <= delta:
                    break
            if test_next == True:
                break
        if test_next == True:
            continue
        return [False, [IS, I_s]]

    return [True, None]

# Returns True if the matrix is considered secure, False otherwise
def algorithm_3(M):
    F = M.base_ring()
    t = M.ncols()

    V = VectorSpace(F, t)
    l = 2*t
    r_limit = floor((t - s) / float(s))

    flag_secure = True
    subspaces_found = []

    # Generate round matrices
    M_round = []
    for j in range(0, l+1):
        M_round.append(M^(j+1))

    I = range(0, s)
    I_powerset = list(sage.misc.misc.powerset(I))[1:]

    for r in range(2, l+1):
        next_r = False
        for I_s in I_powerset:
            IS = None
            res_alg_2 = algorithm_2(M^r)
            if res_alg_2[1] == None:
                continue
            IS = res_alg_2[1][0]
            I_s = res_alg_2[1][1]

            if IS != None and IS.dimension() > 0:
                active_sbox_positions = [[] for _ in range(0, r)]
                active_sbox_positions[0] = I_s
                for j in range(1, r):
                    if IS == subspace_times_matrix(IS, M):
                        next_r = True
                        break
                    IS = subspace_times_matrix(IS, M)
                    for i in range(0, s):
                        # new_basis = [V.basis()[k] for k in range(0, t) if k != i]
                        new_basis = []
                        for k in range(0, t):
                            if k != i:
                                new_basis.append(V.basis()[k])
                        iota_space = V.subspace(new_basis)
                        if IS.intersection(iota_space) != IS:
                            single_iota_space = V.subspace([V.basis()[i]])
                            if IS.intersection(single_iota_space) == single_iota_space:
                                active_sbox_positions[j].append(i)
                            else:
                                next_r = True
                                break
                    if next_r == True:
                        break
                if next_r == True:
                    break
                if active_sbox_positions != [[] for _ in range(0, r)]:
                    flag_secure = False
                    subspaces_found.append([IS, r, active_sbox_positions])
        if next_r == True:
            continue
    
    return [flag_secure, subspaces_found]

def main(row):
    C = matrix.circulant(row)
    print('=== Invariant subspace tests for matrix ===')
    print(C)
    print('Prime: ', hex(C.base_ring().cardinality()))
    print('Is MDS: ', isAllInvertible(C))
    print('Satisfies minpoly condition (this is sufficient, but not necessary): ',
          check_minpoly_condition(C))
    print('Secure by algo 1:', algorithm_1(C))
    print('Secure by algo 2:', algorithm_2(C))
    print('Secure by algo 3:', algorithm_3(C))

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(f'usage: {sys.argv[0]} <prime> <first> <row> <of> <matrix>')
        print(f'  e.g: {sys.argv[0]} 0xffffffff70000001 1 1 2 1 8 32 4 256')
        print(" <prime> can also be 'goldilocks' or 'crandall' to pick those primes")
        exit()

    crandall_prime = 2^64 - 9*2^28 + 1
    goldilocks_prime = 2^64 - 2^32 + 1

    primes = {'crandall': crandall_prime,
              'goldilocks': goldilocks_prime}

    try:
        prime = Integer(sys.argv[1])
        if not prime.is_prime():
            raise ValueError(prime, 'first argument must be prime')
    except TypeError:
        try:
            prime = primes[sys.argv[1]]
        except KeyError:
            print(f'ERROR: first argument must be a prime or one of the '
                  + ' strings {list(primes.keys())}')
            exit()

    field = GF(prime)
    row = vector(field, map(int, sys.argv[2:]))

    main(row)
