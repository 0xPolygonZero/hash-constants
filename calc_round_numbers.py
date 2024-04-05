###
### This file is from
### https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/b5434fd2b2785926dd1dd386efbef167da57c064/code/calc_round_numbers.py
###

from math import *
import sys

# Given p, t, R_P, alpha, M compute the minimal possible value of R_F.
def sat_inequiv_alpha(p, t, R_P, alpha, M):
    n = ceil(log(p, 2))
    if alpha > 0:
        # R_F needs to be bigger or equal to R_F_1, ... R_F_5.
        R_F_1 = 6 if M <= ((floor(log(p, 2) - ((alpha-1)/2.0))) * (t + 1)) else 10 # Statistical
        R_F_2 = 1 + ceil(log(2, alpha) * min(M, n)) + ceil(log(t, alpha)) - R_P # Interpolation
        R_F_3 = (log(2, alpha) * min(M, log(p, 2))) - R_P # Groebner 1
        R_F_4 = t - 1 + log(2, alpha) * min(M / float(t + 1), log(p, 2) / float(2)) - R_P # Groebner 2
        R_F_5 = (t - 2 + (M / float(2 * log(alpha, 2))) - R_P) / float(t - 1) # Groebner 3
        R_F_min = ceil(max(R_F_1, R_F_2, R_F_3, R_F_4, R_F_5))
    else:
        print("Invalid value for alpha!")
        exit(1)

    return R_F_min

# An additional check due to the paper: https://eprint.iacr.org/2023/537.pdf
def extra_check(t, R_F, R_P, alpha, M):
    r_temp = floor(t / 3.0)
    over = (R_F - 1) * t + R_P + r_temp + r_temp * (R_F / 2) + R_P + alpha
    under = r_temp * (R_F / 2) + R_P + alpha
    binom_log = log(comb(int(over), int(under)), 2) # comb(n, m) calculates binomial(n, m)
    if binom_log == float("inf"):
        return True
    
    cost_gb4 = ceil(2 * binom_log) # Paper uses 2.3727, we are more conservative here
    
    return cost_gb4 >= M

def get_sbox_cost(R_F, R_P, N, t):
    return int(t * R_F + R_P)

def get_size_cost(R_F, R_P, N, t):
    n = ceil(float(N) / t)
    return int((N * R_F) + (n * R_P))

def get_depth_cost(R_F, R_P, N, t):
    return int(R_F + R_P)

def find_FD_round_numbers(p, t, alpha, M, cost_function, security_margin):
    n = ceil(log(p, 2))
    N = int(n * t)

    sat_inequiv = sat_inequiv_alpha
    
    R_P = 0
    R_F = 0
    min_cost = float("inf")
    max_cost_rf = 0
    # Brute-force approach
    for R_P_t in range(1, 100):
        R_F_t = sat_inequiv_alpha(p, t, R_P_t, alpha, M)
        if R_F_t % 2 == 1:
            R_F_t += 1
        if (extra_check(t, R_F_t, R_P_t, alpha, M) == True):
            if security_margin == True:
                R_F_t += 2
                R_P_t = int(ceil(float(R_P_t) * 1.075))
            cost = cost_function(R_F_t, R_P_t, N, t)
            if (cost < min_cost) or ((cost == min_cost) and (R_F_t < max_cost_rf)):
                R_P = ceil(R_P_t)
                R_F = ceil(R_F_t)
                min_cost = cost
                max_cost_rf = R_F
    return (int(R_F), int(R_P))

def calc_final_numbers_fixed(p, t, alpha, M, security_margin):
    # [Min. S-boxes] Find best possible for t and N
    n = ceil(log(p, 2))
    N = int(n * t)
    cost_function = get_sbox_cost
    ret_list = []
    (R_F, R_P) = find_FD_round_numbers(p, t, alpha, M, cost_function, security_margin)
    min_sbox_cost = cost_function(R_F, R_P, N, t)
    ret_list.append(R_F)
    ret_list.append(R_P)
    ret_list.append(min_sbox_cost)

    # [Min. Size] Find best possible for t and N
    # Minimum number of S-boxes for fixed n results in minimum size also (round numbers are the same)!
    min_size_cost = get_size_cost(R_F, R_P, N, t)
    ret_list.append(min_size_cost)

    return ret_list # [R_F, R_P, min_sbox_cost, min_size_cost]

def print_round_numbers(primes, monomial_bound = 12, state_widths = [8, 12]):
    print('Number of Poseidon rounds for 128 bit security margin: [full rounds, partial rounds, min. sbox cost, min size cost]\n')
    for name, prime in primes:
        print('  Prime: {} = 0x{:X}'.format('<noname>' if len(name) == 0 else name, prime))
        for alpha in range(3, monomial_bound, 2):
            print(f'    Monomial degree {alpha}', end='')
            if gcd((prime - 1), alpha) != 1:
                print(f'... skipping as {alpha} and p-1 are not coprime')
                continue
            print('')
            for width in state_widths:
                print(f'      state width {width:3}:   ',
                      calc_final_numbers_fixed(prime, width, alpha, 128, True))
        print('')

if __name__ == "__main__":
    crandall_prime = ('Crandall', 2**64 - 9 * 2**28 + 1)
    goldilocks_prime = ('Goldilocks', 2**64 - 2**32 + 1)
    primes = [crandall_prime, goldilocks_prime]

    if len(sys.argv) > 1:
        p = int(sys.argv[1])
        #if not p.is_prime():
        #    print('Argument {p} is not prime')
        primes = [('', p)]

    print_round_numbers(primes)
