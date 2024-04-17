# This is the Poseidon2 implementation used to generate the test for the implementation in Plonky3.
# See https://eprint.iacr.org/2023/323 for an introduction to Poseidon2.

# Given a 4x4 matrix mat_4, create the size x size matrix given by
# Circ(2 mat_4, mat_4, ..., mat_4)
# Used in external layers.
def mat_epsilon(field, size, mat_4):
    block_diag_mat = block_diagonal_matrix([mat_4] * (size//4))
    repeat_mat = matrix(field, size, lambda i, j: mat_4[i%4, j%4])
    return block_diag_mat + repeat_mat

# Create the matrix 1 + Diag(vec)
# Used in internal layers.
def mat_iota(field, size, vec):
    const_mat = matrix.ones(field, size);
    diag_mat  = diagonal_matrix(field, vec);
    return const_mat + diag_mat

# Run a single external (full) round of the protocol.
def full_round(size, mat_e, round_constant, alpha, vec):
    vec += round_constant
    for i in range(size):
        vec[i] = vec[i]**alpha
    
    return mat_e*vec

# Run a single internal (partial) round of the protocol.
def partial_round(size, mat_i, round_constant, alpha, vec):
    vec[0] += round_constant
    
    vec[0] = vec[0]**alpha
    
    return mat_i*vec

# Sometimes we want to work with the internal martix being the MONTY form of the true matrix. (For delayed reduction purposes)
# This can be modelled by shifting mat_i by 2^{-32} as the MONTY constant is 2^32
# This implements Poseidon2 in the usual way if shift is false and the shifted version if shift is true.
def poseidon2(field, size, external_rounds, internal_rounds, external_constants, internal_constants, mat_4, diag_vec, alpha, shift, vec):
    assert(gcd(alpha, field.characteristic() - 1) == 1)
    mat_e = mat_epsilon(field, size, mat_4)
    if shift:
        mat_i = field(2^-32) * mat_iota(field, size, diag_vec)
    else:
        mat_i = mat_iota(field, size, diag_vec)
    
    # The initial linear layer
    output = mat_e*vec
    
    half_external = external_rounds/2
    
    # The first half of the External (full) rounds
    for i in range(half_external):
        output = full_round(size, mat_e, external_constants[i], alpha, output)
    
    # The internal (partial) rounds
    for i in range(internal_rounds):
        output = partial_round(size, mat_i, internal_constants[i], alpha, output)
    
    # The second half of the External (full) rounds
    for i in range(half_external):
        output = full_round(size, mat_e, external_constants[half_external + i], alpha, output)
    
    return output


# The optimal round numbers for a given width, alpha pair.
# Only need these for width = 16, 24 and alpha = 3, 5, 7.
def get_round_numbers(width, alpha):
#   Number of external rounds is always 8.  
    external_rounds = 8
#   Sage does not currently support the python match-case method.
    if (width, alpha) == (16, 3):
        internal_rounds = 20
    elif (width, alpha) == (16, 5):
        internal_rounds = 14
    elif (width, alpha) == (16, 7):
        internal_rounds = 13
    elif (width, alpha) == (24, 3):
        internal_rounds = 23
    elif (width, alpha) == (24, 5):
        internal_rounds = 22
    elif (width, alpha) == (24, 7):
        internal_rounds = 21
    
    return (external_rounds, internal_rounds)
    

# Handy to have a somewhat random but deterministic method to generate constants.
# Not for production but just for testing.
# Good to pick power to be large but shouldn't really matter.
def constants_from_seed(field, seed, power, width, external_rounds, internal_rounds):
    assert(gcd(power, field.characteristic() - 1) == 1)
    EXTERNAL_CONSTANTS = matrix(field, external_rounds, width, lambda i, j: (j + width*i + seed)**power)
    seed2 = seed**(power**2)
    INTERNAL_CONSTANTS = vector(matrix(field, internal_rounds, 1, lambda i, j: (i + seed2)**power))
    return (EXTERNAL_CONSTANTS, INTERNAL_CONSTANTS)

# Generate a Poseidon2 implementation pseudo-randomly from a seed and power.
def poseidon2_from_seed(field, width, alpha, mat_4, diag_vec, shift, seed, power, vec):
    assert(gcd(alpha, field.characteristic() - 1) == 1)
    
    external_rounds, internal_rounds = get_round_numbers(width, alpha)
    
    external_constants, internal_constants = constants_from_seed(field, seed, power, width, external_rounds, internal_rounds)
    
    return poseidon2(field, width, external_rounds, internal_rounds, external_constants, internal_constants, mat_4, diag_vec, alpha, shift, vec)

### Usage examples:

###
### Mersenne31 Field
###


## Start by defining the Field Constants:
# M31_FIELD = GF(2^31 - 1)
# M31_MAT_4 = matrix(M31_FIELD, [[2, 3, 1, 1], [1, 2, 3, 1], [1, 1, 2, 3], [3, 1, 1, 2]])
# M31_MAT_DIAG_16 = [-2,  1,   2,   4,   8,  16,  32,  64, 128, 256, 1024, 4096, 8192, 16384, 32768, 65536]
# M31_MAT_DIAG_24 = [-2,  1,   2,   4,   8,  16,  32,  64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304]
# M31_ALPHA = 5

# Next, generate parameters:
# P = Primes()
# M31_power = P.next(11111)
# set_random_seed(11111)
# M31_seed = M31_FIELD.random_element()

# This yields:
# M31_power = 11113
# M31_seed = 389183646

# Now generate some random inputs and run Poseidon2 width 16:
# set_random_seed(16)
# rand_input_M31_16 = vector([M31_FIELD.random_element() for t in range(16)])
# rand_ouput_16 = poseidon2_from_seed(M31_FIELD, 16, M31_ALPHA, M31_MAT_4, M31_MAT_DIAG_16, False, M31_seed, M31_power, rand_input_M31_16)

# We find:
# rand_input_M31_16 = [894848333, 1437655012, 1200606629, 1690012884, 71131202, 1749206695, 1717947831, 120589055, 19776022, 42382981, 1831865506, 724844064, 171220207, 1299207443, 227047920, 1783754913]
# rand_output_M31_16 = [687671392, 187739990, 474872297, 1025723782, 1958464721, 1004876398, 972043176, 1231017992, 1815473754, 997812498, 1891950360, 94240126, 1834774779, 2146393033, 1194588914, 1694651572]

# Can do the same for Poseidon2 width 24:
# set_random_seed(24)
# rand_input_M31_24 = vector([M31_FIELD.random_element() for t in range(24)])
# rand_ouput_24 = poseidon2_from_seed(M31_FIELD, 24, M31_ALPHA, M31_MAT_4, M31_MAT_DIAG_24, False, M31_seed, M31_power, rand_input_M31_24)

# We get:
# rand_input_M31_24 = [886409618, 1327899896, 1902407911, 591953491, 648428576, 1844789031, 1198336108, 355597330, 1799586834, 59617783, 790334801, 1968791836, 559272107, 31054313, 1042221543, 474748436, 135686258, 263665994, 1962340735, 1741539604, 2026927696, 449439011, 1131357108, 50869465]
# rand_output_M31_24 = [1647982233, 532201789, 2015259638, 680414033, 1176515591, 1163262902, 2052469886, 679834118, 2079721199, 1936663286, 66010933, 593633651, 69437652, 889870443, 1128148983, 1865068789, 649133836, 1434401472, 648402879, 1081496263, 388204549, 380976594, 1146523540, 841407217]