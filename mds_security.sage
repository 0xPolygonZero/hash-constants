import random
import time

if len(sys.argv) != 5:
    print("Usage: <script> <n> <t> <sample_size> <seed>")
    exit()

n = int(sys.argv[1])
t = int(sys.argv[2])
sample_size = int(sys.argv[3])
random_seed = int(sys.argv[4])
s = 1
r = floor((t - s) / float(s))
use_seed = True
single_matrix_test = True
if single_matrix_test == True:
    use_seed = False
if use_seed == True:
    set_random_seed(random_seed)
    random.seed(random_seed)

def get_prime(bit_length):
    while True:
        p = random_prime((2^bit_length) - 1, false, 2^(bit_length - 1))
        # if gcd(p - 1, 3) == 1:
        #     return p
        return p

# use_seed = 1
# set_random_seed(use_seed)
# random.seed(use_seed)

if t == 8 or t == 12:
    p = 2^64 - 2415919103  # Crandall prime
    k = GF(p)
    prime = p
    F = k
    MS = MatrixSpace(F, t, t)
else:
    prime = get_prime(n)
    F = GF(prime)
    MS = MatrixSpace(F, t, t)

def print_matrix_format(M_int, n, t):
    print("n:", n)
    print("t:", t)
    print("N:", (n * t))
    hex_length = int(ceil(float(n) / 4)) + 2 # +2 for "0x"
    print("Prime number:", "0x" + hex(prime))
    print("MDS matrix (rows):")
    for i in range(0, t):
        print(["{0:#0{1}x}".format(entry, hex_length) for entry in M_int[i]])

def matrix_entries_to_int(M, t):
    M_int = []
    for i in range(0, t):
        M_int.append([])
        for j in range(0, t):
            M_int[i].append(int(M[i, j]))
    return M_int

def create_mds_p_spec(n, t):
    M = matrix(F, t, t)
    xs = []
    ys = []
    
    for i in range(0, t):
        xs.append(F(i))
        ys.append(F(t + i))
    
    for i in range(0, t):
        for j in range(0, t):
            entry = (xs[i] + ys[j])^(-1)
            M[i, j] = entry
    return M

def isAllInvertible(M, t):
    # Test all square submatrices for invertibility
    counter = 0
    all_invertible = True
    for i in range(2, t):
        choices_i = Combinations(range(0, t), i)
        for m in range(0, choices_i.cardinality()):
            for n in range(0, choices_i.cardinality()):
                M_sub = M[choices_i[m], choices_i[n]]
                is_inv = M_sub.is_invertible()
                all_invertible = all_invertible and is_inv
                #if is_inv == False:
                    #print("FALSE")
                    #print(M_sub)
                counter += 1
    #print("Submatrices checked:", counter)
    return all_invertible

def create_mds_p(n, t):
    M = matrix(F, t, t)

    # Sample random distinct indices and assign to xs and ys
    while True:
        flag = True
        rand_list = [F(ele) for ele in random.sample(range(0, prime), 2*t)]
        xs = rand_list[:t]
        ys = rand_list[t:]
        # xs = [F(ele) for ele in range(0, t)]
        # ys = [F(ele) for ele in range(t, 2*t)]
        for i in range(0, t):
            for j in range(0, t):
                if (flag == False) or ((xs[i] + ys[j]) == 0):
                    flag = False
                else:
                    entry = (xs[i] + ys[j])^(-1)
                    M[i, j] = entry
        if flag == False:
            continue
        return M

def get_eigenvalues(mat):
    #print(mat.charpoly().roots())
    eigenvalues = [mat.charpoly().roots()[i][0] for i in range(0, len(mat.charpoly().roots()))]
    return eigenvalues

def get_eigenvectors(mat):
    # Computing the right eigenvectors
    eigenvalues = get_eigenvalues(mat)
    eigenvectors = []
    for i in range(0, len(eigenvalues)):
        #print(eigenvalues[i])
        eig_basis = ((mat - eigenvalues[i] * matrix.identity(F, t)).right_kernel()).basis()
        for vec in eig_basis:
            eigenvectors.append(vec)
    return eigenvectors

def generate_vectorspace(round_num, M, M_round):
    V = VectorSpace(F, t)
    if round_num == 0:
        return V
    elif round_num == 1:
        return V.subspace(V.basis()[s:])
    else:
        mat_temp = matrix(F)
        for i in range(0, round_num-1):
            add_rows = []
            for j in range(0, s):
                add_rows.append(M_round[i].rows()[j][s:])
            mat_temp = matrix(mat_temp.rows() + add_rows)
        r_k = mat_temp.right_kernel()
        extended_basis_vectors = []
        for vec in r_k.basis():
            extended_basis_vectors.append(vector([0]*s + list(vec)))
        S = V.subspace(extended_basis_vectors)

        return S

def subspace_times_matrix(subspace, M):
    V = VectorSpace(F, t)
    subspace_basis = subspace.basis()
    new_basis = []
    for vec in subspace_basis:
        new_basis.append(M * vec)
    new_subspace = V.subspace(new_basis)
    return new_subspace

# Returns True if the matrix is considered secure, False otherwise
def algorithm_1(M):

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

def check_minpoly_condition(M):
    max_period = 2*t
    all_fulfilled = True
    M_temp = M
    for i in range(1, max_period + 1):
        if not ((M_temp.minimal_polynomial().degree() == t) and (M_temp.minimal_polynomial().is_irreducible() == True)):
            all_fulfilled = False
            break
        M_temp = M * M_temp
    return all_fulfilled

def construct_test_matrix_inv_active():

    M = None
    while M == None or M.is_invertible() == False:
        # New matrix
        M = matrix(F, t, t)

        # Random a, b, c, ...
        while M[1, 0] == F(0):
            M[1, 0] = F.random_element()
        for i in range(2, t):
            M[i, 0] = F.random_element()

        # Random M_{0,2}, M_{1,2}, ...
        for i in range(0, t):
            for j in range(2, t):
                M[i, j] = F.random_element()
        
        # Second column
        prod_sum = 1
        for i in range(2, t):
            prod_sum -= M[0, i]*M[i, 0]
        M[0, 1] = prod_sum/M[1, 0]
        for i in range(1, t):
            prod_sum = 0
            for j in range(2, t):
                prod_sum -= M[i, j]*M[j, 0]
            M[i, 1] = prod_sum/M[1, 0]

    return M

def construct_test_matrix_inv_active_gen():
    t_custom = 4
    M = None
    while M == None or M.is_invertible() == False:
        # New matrix
        M = matrix(F, t_custom, t_custom)

        # Random a, b, c, ...
        M[0, 0] = F(1)
        while M[1, 0] == F(0):
            M[1, 0] = F.random_element()
        for i in range(2, t_custom):
            M[i, 0] = F.random_element()

        # Random M_{0,2}, M_{1,2}, ...
        for i in range(0, t_custom):
            for j in range(2, t_custom):
                M[i, j] = F.random_element()
        
        # Second column
        prod_sum = 0
        for i in range(2, t_custom):
            prod_sum -= M[0, i]*M[i, 0]
        M[0, 1] = prod_sum/M[1, 0]
        for i in range(1, t_custom):
            prod_sum = -M[i, 0]
            for j in range(2, t_custom):
                prod_sum -= M[i, j]*M[j, 0]
            M[i, 1] = prod_sum/M[1, 0]

    return M

def construct_test_matrix_e2():
    t_custom = 8
    M = None
    while M == None or M.is_invertible() == False:
        # New matrix
        # Copy M_1 (CHANGE GLOBAL t)
        global t
        old_t = t
        t = 4
        M_1 = construct_test_matrix_521()
        t = old_t
        M = matrix(F, t_custom, t_custom)
        for i in range(0, 4):
            for j in range(0, 4):
                M[i,j] = M_1[i,j]

        # Construct M_2
        MS_reduced = MatrixSpace(F, 4, 4)
        M_2 = MS_reduced.random_element()
        for i in range(0, 4):
            M_2[i,0] = M_2[i,2]
            M_2[i,1] = M_2[i,3]
        for i in range(0, 4):
            for j in range(0, 4):
                M[i,j+4] = M_2[i,j]

        # Construct M_3
        M_3 = MS_reduced.random_element()
        for i in range(0, 4):
            M_3[i,0] = F(0)
            M_3[i,1] = -M_3[i,2] - M_3[i,3]
        for i in range(0, 4):
            for j in range(0, 4):
                M[i+4,j] = M_3[i,j]

        # Construct M_4 (circ(a, b, c, d) with two eigenvalues)
        init_vec = vector([F.random_element() for _ in range(0, 4)])
        M_4 = matrix.circulant(init_vec)
        while len(get_eigenvalues(M_4)) != 2:
            init_vec = vector([F.random_element() for _ in range(0, 4)])
            M_4 = matrix.circulant(init_vec)
        for i in range(0, 4):
            for j in range(0, 4):
                M[i+4,j+4] = M_4[i,j]
    
    return M

def custom_random_matrix():
    range_limit = 4
    while True:
        #random_list = random.sample(xrange(0, range_limit), t*t)
        random_list = [random.randint(1, range_limit) for _ in range(0, t*t)]
        #print(random_list)
        M = matrix(F, t, t)
        for i in range(0, t):
            for j in range(0, t):
                M[i, j] = random_list[i*t + j]
        if isAllInvertible(M, t) == True:
            return M

test_algorithms = [algorithm_1]
discard_by_algorithms = []

if single_matrix_test == True:
    # Test a single matrix
    if t == 8:
        p = 2^64 - 2415919103  # Crandall prime
        k = GF(p)
        row0 = [4, 1, 2, 256, 16, 8, 1, 1] # 5, 7, 256
        C = matrix.circulant(row0)
        C = C.change_ring(k)
        print("=== Small entries matrix ===")
        print(C)
        print("Prime:", p)
        print("Is MDS: ", isAllInvertible(C, t))
        print("Secure by algo 1:", algorithm_1(C))
        print("Secure by algo 2:", algorithm_2(C))
        print("Secure by algo 3:", algorithm_3(C))

    if t == 12:
        p = 2^64 - 2415919103  # Crandall prime
        k = GF(p)
        row0 = [1024, 8192, 4, 1, 16, 2, 256, 128, 3, 32, 1, 1]
        C = matrix.circulant(row0)
        C = C.change_ring(k)
        print("=== Small entries matrix ===")
        print(C)
        print("Prime:", p)
        print("Is MDS: (yep!)")
        print("Secure by algo 1:", algorithm_1(C))
        print("Secure by algo 2:", algorithm_2(C))
        print("Secure by algo 3:", algorithm_3(C))

    if t == 3:
        # Iterative subspace trails with active S-boxes, but no invariant ones
        M = matrix(F, [[F(0), F(1), F(-1)], [F(1), F(-2), F(1)], [F(1), F(-4), F(2)]])
        print("=== Iterative Subspace Trails with Active S-Boxes, but no Invariant Ones ===")
        print(M)
        print("Prime:", prime)
        print("Result (expected: True):", algorithm_1(M))
        print("Result (expected: True):", algorithm_2(M))
        print("Result (expected: False):", algorithm_3(M))
    if t == 4:
        # Invariant subspace trails with active S-boxes
        M = construct_test_matrix_inv_active()
        print("=== Invariant Subspace Trails with Active S-Boxes ===")
        print(M)
        print("Prime:", prime)
        print("Result (expected: False):", algorithm_2(M))
        # vs = algorithm_2(M)[1][0]
        # V = VectorSpace(F, t)
        # V_sub = V.subspace([vector([F(1), F(0), F(0), F(0)]), vector([F(0), F(M[1,0]), F(M[2,0]), F(M[3,0])])])
        # print("Spaces match:", vs == V_sub)
        # Invariant subspace trails with active S-boxes (generalization)
        M = construct_test_matrix_inv_active_gen()
        print("=== Invariant Subspace Trails with Active S-Boxes (Generalization) ===")
        print(M)
        print("Prime:", prime)
        print("Result (expected: False):", algorithm_2(M))
        # vs = algorithm_2(M)[1][0]
        # V = VectorSpace(F, t)
        # V_sub = V.subspace([vector([F(1), F(0), F(0), F(0)]), vector([F(1), F(M[1,0]), F(M[2,0]), F(M[3,0])])])
        # print("Spaces match:", vs == V_sub)
    # E.2 (ONLY FOR GF(11))
    # if t == 5:
    #     if prime != 11:
    #         print("Example E.2 works only in GF(11)! Set n=4 and run again.")
    #         exit()
    #     M = matrix(F, t, [[F(0),F(3),F(5),F(0),F(7)],
    #                     [F(0),F(4),F(5),F(0),F(6)],
    #                     [F(1),F(4),F(8),F(3),F(0)],
    #                     [F(1),F(4),F(3),F(7),F(1)],
    #                     [F(1),F(4),F(9),F(9),F(4)]])
    #     print("=== E.2 ===")
    #     print(M)
    #     print("Prime:", prime)
    #     print("Result (Algorithm 1) (expected: True):", algorithm_1(M))
    #     print("Result (Algorithm 2) (expected: False):", algorithm_2(M))
    # E.2
    # if t == 8:
    #     M = construct_test_matrix_e2()
    #     # Example matrix from paper (paper uses GF(131))
    #     # M = matrix(F, t, t, [[1, 112, 125, 116, 59, 19, 59, 19],
    #     #     [51, 128, 72, 36, 18, 34, 18, 34],
    #     #     [94, 54, 7,  22, 41, 30, 41, 30],
    #     #     [55, 46, 19, 12, 86, 70, 86, 70],
    #     #     [0, 38, 130, 94, 84, 127, 7, 42],
    #     #     [0, 35, 92, 4, 42, 84, 127, 7],
    #     #     [0, 128, 129, 5, 7, 42, 84, 127],
    #     #     [0, 105, 108, 49, 127, 7, 42, 84]])
    #     print("=== E.2 ===")
    #     print(M)
    #     print(latex(M))
    #     print("Prime:", prime)
    #     print("Result (Algorithm 1) (expected: False):", algorithm_1(M))
    #     print("Result (Algorithm 2) (expected: True):", algorithm_2(M))
    #     print("Result (Algorithm 3) (expected: False):", algorithm_3(M))
    #     # vs = algorithm_3(M)[1][0]
    #     # V = VectorSpace(F, t)
    #     # V_sub = V.subspace([vector([F(0), F(0), F(0), F(0), F(1), F(0), F(-1), F(0)]), vector([F(0), F(0), F(0), F(0), F(0), F(1), F(0), F(-1)])])
    #     # #V_sub = V.subspace([vector([F(1), F(0), F(0), F(0), F(0), F(0), F(0), F(0)]), vector([F(0), F(1), F(1), F(1), F(0), F(0), F(0), F(0)]), vector([F(0), F(1), F(1), F(1), F(0), F(1), F(0), F(-1)])])
    #     # print("Spaces match:", vs == V_sub)
    exit()

print("--- GF(p) ---")
print("--- Testing algorithms ---")
for test_algorithm in test_algorithms:
    print(test_algorithm.__name__)
print("--- Discarding by algorithms ---")
for discard_algorithm in discard_by_algorithms:
    print(discard_algorithm.__name__)
print("--- n = " + str(n) + ", t = " + str(t) + " ---")
print("Prime number:", prime)

## Poseidon matrices
print("Spec matrices...")
num_secure = 0
num_vulnerable = 0
num_tests = 1 # This matrix is deterministically defined for (n, t) according to the paper
print("Tests:", num_tests)
for i in range(0, num_tests):
    next_test = False
    M = create_mds_p_spec(n, t)
    # Discard first (discard matrices which are vulnerable w.r.t. an algorithm in discard_by_algorithms)
    for discard_algorithm in discard_by_algorithms:
        if discard_algorithm(M)[0] == False:
            next_test = True
            break
    if next_test == True:
        num_secure += 1
        continue
    # Current test algorithms
    all_secure = True
    for test_algorithm in test_algorithms:
        if test_algorithm(M)[0] == False:
            all_secure = False
            break
    if all_secure == True:
        num_secure += 1
    else:
        num_vulnerable += 1

print("[T1] Secure:", num_secure / (float(num_tests) / 100.0))
print("[T1] Vulnerable:", num_vulnerable / (float(num_tests) / 100.0))

### Random invertible matrices
print("Random invertible matrices...")
num_secure = 0
num_vulnerable = 0
num_tests = sample_size
print("Tests:", num_tests)
time_start = time.time()
for i in range(0, num_tests):
    next_test = False
    M = MS.random_element()
    while M.is_invertible() == False:
        M = MS.random_element()
    # Discard first (discard matrices which are vulnerable w.r.t. an algorithm in discard_by_algorithms)
    for discard_algorithm in discard_by_algorithms:
        if discard_algorithm(M)[0] == False:
            next_test = True
            break
    if next_test == True:
        num_secure += 1
        continue
    # Current test algorithms
    all_secure = True
    for test_algorithm in test_algorithms:
        if test_algorithm(M)[0] == False:
            all_secure = False
            break
    if all_secure == True:
        num_secure += 1
    else:
        num_vulnerable += 1

time_end = time.time()
print("[T2] Secure:", num_secure / (float(num_tests) / 100.0))
print("[T2] Vulnerable:", num_vulnerable / (float(num_tests) / 100.0))
print("[T2] Average Execution Time:", (time_end - time_start) / float(num_tests))

### Random Cauchy matrices
print("Random Cauchy matrices...")
num_secure = 0
num_vulnerable = 0
num_tests = sample_size
print("Tests:", num_tests)
time_start = time.time()
for i in range(0, num_tests):
    next_test = False
    M = create_mds_p(n, t)
    # Discard first (discard matrices which are vulnerable w.r.t. an algorithm in discard_by_algorithms)
    for discard_algorithm in discard_by_algorithms:
        if discard_algorithm(M)[0] == False:
            next_test = True
            break
    if next_test == True:
        num_secure += 1
        continue
    # Current test algorithms
    all_secure = True
    for test_algorithm in test_algorithms:
        if test_algorithm(M)[0] == False:
            all_secure = False
            break
    if all_secure == True:
        num_secure += 1
    else:
        num_vulnerable += 1

time_end = time.time()
print("[T3] Secure:", num_secure / (float(num_tests) / 100.0))
print("[T3] Vulnerable:", num_vulnerable / (float(num_tests) / 100.0))
print("[T3] Average Execution Time:", (time_end - time_start) / float(num_tests))
