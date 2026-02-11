
###
# ###
# Accompanying script for the paper
# "Round-Optimal Identity-Based Blind Signature from Lattice Assumptions"
# The script provides security estimates for the id-based blind signature construction
# as well as estimate for the signature size
# ###
# Security and size estimates derived from [LNP22] ***version 1***
#(dated March 22, 2022), available at https://eprint.iacr.org/2022/
# The script is used from [AKSY 22] available at 
# https://gitlab.com/ElenaKirshanova/onemoresis_estimates/ with several necessary changes
#to accomodate our scheme
 


from MSIS_security import SIS_optimize_attack, SIS_l2_cost, SIS_linf_cost
from MLWE_security import MLWE_optimize_attack, LWE_primal_cost, LWE_dual_cost
from model_BKZ import svp_classical
from math import sqrt, floor, ceil, log, exp
from sympy import nextprime
#from estimator_ntru import combined_attack_prob
#from one_more_sis_security import combinatorial, lattice_attack

pi = 3.1415926535897932384626433832795028841971693993751058209749

###############################################################################
#
# variables that define LWE/SIS hardness differ in [LNP22] from Fig.10 to Fig.6
# below is the transition dictionary.
# We instantiate Fig.10, monitor how the parameters get updated to Fig.6
# and instantiate Thm. 4.2. with the updated parameters. The result is given
# in the HARDNESS section below.
#
#               |  Fig.10    |    Fig.8                 | Fig.6
# ______________________________________________________________________________
#   A1.nrows    |   n        |    n                     |    n
#   s1.dim      |   m1       |    m1+ve                 |  m1+ve
#   s2.dim      |   m2       |    m2                    |    m2
#   m.dim       |   ell      |    ell+(256/d+1)*{0,1,2} | Fig.8+lamb/2
# ______________________________________________________________________________
#                               HARDNESS (Thm. 4.2)
#                       here m1 and ell are those of Fig.10
#   LWE:        |
#      dim      |   m2-n-ell-{0,1,2}*(256/d+1)-lamb/2-1
#      Nsamples |   n+ell+{0,1,2}*(256/d+1)+lamb/2+1
# ______________________________________________________________________________
#   SIS:
#    hight of A |   n
#    width of A |   m1+m2+ve
#
#
#   multiplication by a set in m.dim for Fig.8 and hardness must be
#   understood as follows:  we multiply
#       by 0 if ve=vd=0,
#       by 1 if only one of {ve,vd} is at least 1,
#       by 2 if both {ve,vd} are at least 1
#
#
############################################################################
# the basic signature scheme [Sec 4.1] parameters 
# we adapt a modified version of the scheme in [BEP+ 21] 
# Available at: 
#####################

print(' ------------------ GPV (G-trapdoor) signature scheme  ------------------ ' )

lamda = 128 #192 # security parameter
q_s = 536871029 #524309 #16777333 #1073741909
n_s = 256 #512 # degree of the ring R_{q_s}
k = 30 #25#31 # bitlength of q
# This is not the 30 bit prime modulus chosen by [BEP+ 21]. 
# We chose the smallest 31 bit prime that satisfies:
# (1) q_s mod 8 = 5 -> x^128+1 splits into two nice factors (see Lemma 2.5 in LNP)
# (2) q_s > 2^{128/lamda}
# the smallness of q_s makes beta_y small, hence milder lower bound on q_zk later

# following computation by [BEP+ 21] Section A.2 (archive version)

d = 3 #2 # no. of rows of A
m = 2*d + d*k

# TrapGen -> (A,R)
eps= 2**(-lamda) 
sigma_bound = sqrt((log((2*n_s*2*d)/eps)/pi)) # lower bound on trapgen sd
sigma = 71 #135 # sd of TrapGen so that A looks like uniform
assert(sigma > sigma_bound)
print('sigma:', sigma)

# keygen security, the trapdoor must look uniform
TrapGen_LWE_security = MLWE_optimize_attack(q_s, n_s*d, n_s*d, sigma, cost_attack = LWE_primal_cost, cost_svp = svp_classical, verbose = False)
print('TrapGen security as LWE:', TrapGen_LWE_security)


# DelTrap -> R_id
alpha_bound = 10*sqrt(log(2*k*(1+eps**(-1)))/ pi) # lower bound on alpha 
alpha = 89
assert(alpha > alpha_bound)
print('alpha:',alpha)
s1_R_bound = 0.39*sigma*(sqrt(2*n_s*d)+sqrt(n_s*d*k)+4.7) # upper bound on the spectral norm s1 [BEP+ 21]
s1_R = 55
assert(s1_R < s1_R_bound)
eta = sqrt((log((2*n_s*m)/eps)/pi)) 
s_bound = sqrt(((alpha**2) +1)*(s1_R**2)+(eta**2)) # lower bound on SamplePre sd 
s = 4992 #1944 #1943 # standard deviation for SamplePre in sign
assert(s > s_bound)
print('s:', s)
t = 12 #15  #gaussian tailcut
length_bound = t*s*sqrt(m*n_s)

DelTrap_SIS_security = SIS_optimize_attack(q_s, n_s*m, n_s*d, length_bound, cost_attack = SIS_l2_cost, cost_svp = svp_classical, verbose = False)
print('signing key generation as SIS :', DelTrap_SIS_security)


# signature forgery decurity
m_id = m + d*k # no. of cols in A_id # dimension of the sign y
s1_Rid_bound = 0.39*s*(sqrt(n_s*m)+sqrt(n_s*d*k)+4.7) # upper bound on the spectral norm s1 [BEP+ 21]
s1_Rid = 18 
assert(s1_Rid < s1_Rid_bound)
eta_prime = sqrt((log((2*n_s*m_id)/eps)/pi)) 
s_prime_bound = sqrt(((alpha**2) +1)*(s1_Rid**2)+(eta_prime**2)) # lower bound on SamplePre
s_prime = 3297 #1287 #1258 # standard deviation for SamplePre in DelTrap
assert(s_prime > s_prime_bound)
print('s_prime:', s_prime)
beta_y = t*s_prime*sqrt(m_id*n_s) # l2 norm bound on signature 
print('beta_y:', beta_y)

sign_forgery = SIS_optimize_attack(q_s, n_s*m_id, n_s*d, beta_y, cost_attack = SIS_l2_cost, cost_svp = svp_classical, verbose = False)
print('signature forgery as SIS:', sign_forgery)

#############################
# message blinding

row = 3 #2 # no. of rows in A_1
sigma_x = 60999 #27999 # each coefficient of x sampled with standard deviation sigma_x 
beta_x = sigma_x*sqrt(m*n_s) # l2 norm bound on x
print('beta_x:', beta_x)

blinding_security = SIS_optimize_attack(q_s, n_s*m, n_s*d, beta_x, cost_attack = SIS_l2_cost, cost_svp = svp_classical, verbose = False)
print('blinding security as SIS:', blinding_security)


############################
# Ajtai commitment to x parameters

# Hiding & Binding
m_xr_bound = row + ((2*lamda)/log(q_s,2)) # to apply left over hash lemma
m_xr = 97 #67 # dimension of [x  r]^T
assert(m_xr > m_xr_bound)
assert(m_xr > m)
eta_r = sqrt((log((2*n_s*m_xr)/eps)/pi)) # smooting parameter # lower bound on r 
sigma_r = 195 #5.6 #6.8 # each coefficient of r sampled with standard deviation sd_r
assert(sigma_r > eta_r)
dim_r = m_xr - m #dimension of of randomness r
print("dimension of r:", dim_r)

beta_r = sigma_r*sqrt(dim_r*n_s) # l2 norm bound on r
print('beta_r:', beta_r)
beta = sqrt((beta_x**2) + (beta_r**2))

commitment_security = SIS_optimize_attack(q_s, n_s*m_xr, n_s*row, beta, cost_attack = SIS_l2_cost, cost_svp = svp_classical, verbose = False)
print('commitment security as SIS:', commitment_security)


#####################
# ZK proof    (see LNP22, Fig. 10)
# zk_secret s_zk = [y| x| r]
# E1 = [100] (for y small)
# norm bound for E1: beta_y
# E2 = [010]  (for x small)
# norm bound for E2: beta_x
# E3 = [001]  (for r small)
# norm bound for E3: beta_r
# Plus two linear equations, which don't affect proof size or security
#####################


print(' ------------------ ZK proof  area  ------------------ ' )

# in the unforgeability proof, we replace h' by y1*h+y2, where y1 and y2 have
# infinity norms tau_x, like x1 and x2.
# the extracted one-more-isis solution is (s1+x1y1 , s2+x1*y2+x2)
# below is its norm, which is hopefully not too large compared to beta
#extracted_beta = sqrt(beta**2 + tau_x**4 * n_F*n_F*2*2+ tau_x**2 * n_F*2)

n_zk = 128
ve = 3 # number of norm equations
vd = 0 # number of inf-norm equations
norm_x_zk = sqrt(ve)*sqrt(n_zk)

# norm of [s1 & x1 & x2 & s] as per Thm. 5.3 (Fig. 10), in Thm.4.2 (Fig. 6) y1||z1 -> [y1 || z1 ||x]
# divide beta**2 by two since we commit only to the first half of (y1,y2)
norm_szk = sqrt( (beta_y**2) + (beta_x)**2 + (beta_r)**2 )
#norm of the norm vector s^(e) = [y | x | r | x_zk)
alpha_e = sqrt(norm_szk**2  + norm_x_zk**2)
print('alpha_e:', alpha_e)

gamma_e = 3 # influences the repetition factor
t       = 1.64 # as in LNP22 Thm 5.3
kappa   = 2 # taken from LNP22
Be = 2*sqrt(256/26)*t*gamma_e*sqrt(337)*(alpha_e) # as per Thm. 5.3
ce = n_zk*(round(n_s*m_id/n_zk)+1+round(n_s*m/n_zk)+1+round(n_s*1/n_zk)+1) # Fig 10
lb1 = floor(Be*41*ce)+1                      # Thm. 5.3
lb2 = floor(Be**2+sqrt(ve*n_zk)*Be)+1        # Thm. 5.3
beta1 = beta_y                               # ||y||
beta2 = beta_x                               # ||x||
beta3 = beta_r                               # ||r||
lb3 = floor(Be**2+2*(max(beta1, beta2, beta3))**2)  # Thm. 5.3
lb = max(lb1, lb2, lb3)


print(2*sqrt(256/26)*t*gamma_e*sqrt(337))
print('ce:', ce, 'log(Be):', log(Be,2))
print()
print('lb1:', log(lb1,2), 'lb2:', log(lb2,2), 'lb3:',log(lb3,2))


# q_zk = q_S * q1  
# We want
# (0) q1 prime
# (1) q1 >= q_s
# (2) q_zk >= lb
# (3) (x^d+1) has two factors mod q_zk

q1 = 99327398173 #18326199821 #10103076346093 #778151346517 # 14909121133 #chosen as next prime to ceil(lb/q) such that it is congruent to 5 mod 8
q_zk = q_s*q1
assert(q_zk>=lb)
assert(q1>=q_s)
print('q1:', q1)
print('log(q_zk):', log(q_zk,2), 'q_zk:', q_zk)

l = 2 # number of factors of (x^128+1) mod q1
kappa_bound = (1./(2*sqrt(l)))*(q_s)**(1./l)
# The sigma_{-1} row of Fig 3 can be used if kappa < kappa_bound
assert(kappa < kappa_bound)
print('kappa bound:', kappa_bound)

#summary of all the params
#  We are somewhat free to choose n and m2
#  For MLWE to make sense we require m2>n+{0,1,2}*(256/d+1)+lamb/2+ell,

gamma = 2**24 # as per [AKSY 22]
D_bound = floor(log(gamma/ (kappa*n_zk), 2) + 1) 
D = D_bound - 1 # D < D_bound as per LNP22 (paragraph Dilithium compression)
assert(D < D_bound)
print('D:', D)
param_zk = {'n_zk': 128, # ring dimension
             'n': 13, # determines SIS hardness
             'm1': 569, #787 # m1 = len(y) + len(x) + len(r) + ve = 4.(m+dk) + 4.m + 4.1 + 3 = 569
             'm2': 39, # m2 > n+{0,1,2}*(256/n_zk+1)+lambda/2+ell determines LWE hardness
             'q_zk': q_zk, # ZK modulus
             'lambda': 2, # taken from LNP22
             'nu': 1,  # taken from LNP22
             'eta': 59, # taken from LNP22
             'gamma1': 10, # influences the repetition factor
             'gamma2': 1.5,  # influences the repetition factor
             'gammae': gamma_e,  # influences the repetition factor
             'alphae':(alpha_e)**2, # norm of the norm vector s^(e) = [y | x | r | x_zk)
             'norm_s1': norm_szk**2, # norm of the zk secret
             'ell': 0,
             've': 3,  # number of exact norm proofs
             'vd': 0,  # number of approxumate norm proofs
             'D': D, # cutting D low-order bits from t_A
             'gamma': gamma # cutting log_gamma bits from w
             }
scalar = 1 #scalar is the number from the set {0,1,2} (see formula above) that depends on the number ve and vd
assert(param_zk['m2']>param_zk['n']+scalar*(256/n_zk+1)+param_zk['lambda']/2+param_zk['ell']+2)


def SIS_security(paramset):
    # as per [LNP22] Thm 4.2 (see above dictionary)
    n_zk  = paramset['n_zk']
    q_zk  = paramset['q_zk']
    n  = paramset['n']
    m1 = paramset['m1']
    m2 = paramset['m2']
    norm_s1  = paramset['norm_s1']
    ve = paramset['ve']
    nu = paramset['nu']
    eta    = paramset['eta']
    gamma1 = paramset['gamma1']
    gamma2 = paramset['gamma2']
    gamma  = paramset['gamma'] #Dilithium-G - compression
    D  = paramset['D'] #Dilithium-G - compression
    sigma1 = gamma1*eta*sqrt(norm_s1+ve*n_zk)
    # sqrt(norm_s1/2+ve*n_zk)is an upper bound on the ell2-norm of ABDLOP's s1.
    sigma2 = gamma2*eta*nu*sqrt(m2*n_zk)
    # sqrt(m2*nu*n_zk) is an upper bound on the ell2-norm of ABDLOP's s2.

    # Bound on the SIS solution. From LNP22
    B1 = 2*sigma1*sqrt(2*m1*n_zk)
    B2 = 2*sigma2*sqrt(2*m2*n_zk)+2**D*eta*sqrt(n*n_zk)+gamma*sqrt(n*n_zk)
    BSIS = 4*eta*sqrt(B1**2+B2**2)

    assert(BSIS<q_zk)
    # maximal width
    max_w = (m1+m2)*n_zk
    h = n*n_zk
    m_pc, b_pc, c_pc = SIS_optimize_attack(q_zk, max_w, h, BSIS, cost_attack=SIS_l2_cost, cost_svp=svp_classical, verbose=False)
    return (m_pc, b_pc, c_pc)
    
def LWE_security(paramset, attack=LWE_primal_cost):
    # as per [LNP22] Thm 4.2 (see above dictionary)
    n_zk  = paramset['n_zk']
    q_zk  = paramset['q_zk']
    n  = paramset['n']
    m1 = paramset['m1']
    m2 = paramset['m2']
    nu = paramset['nu']
    ve = paramset['ve']
    vd = paramset['vd']
    ell = paramset['ell']
    lamb = paramset['lambda']
    assert(lamb%2 == 0)
    if vd==0 and ve==0:
        scalar = 0
    elif vd>=1 and ve>=1:
        scalar = 2
    else: scalar = 1
    
    ell_updated = ell+(round(256/n_zk)+1)*scalar+round(lamb/2)
    mLWE = (n+ell_updated+1)*n_zk
    nLWE = (m2-n-ell_updated-1)*n_zk
    m_pc_, b_pc_, c_pc_ = MLWE_optimize_attack(q_zk, nLWE, mLWE, nu, cost_attack=attack, cost_svp=svp_classical, verbose=False)
    return(m_pc_, b_pc_, c_pc_)


def proof_size(paramset):
    n_zk  = paramset['n_zk']
    q_zk  = paramset['q_zk']
    n  = paramset['n']
    m1 = paramset['m1']
    m2 = paramset['m2']
    assert(m1>=512/n_zk)
    assert(m2>=512/n_zk)      # as per Thm 4.3 [LNP 22]
    alphae  = paramset['alphae']
    norm_s1 = paramset['norm_s1']
    ell = paramset['ell']
    ve = paramset['ve']
    vd = paramset['vd']
    lamb   = paramset['lambda']
    assert(lamb%2 == 0)    # as per Section 4.4 [LNP 22]
    if vd==0 and ve==0:
        scalar = 0
    elif vd>=1 and ve>=1:
        scalar = 2
    else: scalar = 1
    D     = paramset['D']
    nu = paramset['nu']
    eta    = paramset['eta']
    gamma1 = paramset['gamma1']
    gamma2 = paramset['gamma2']
    gammae = paramset['gammae']
    sigma1 = gamma1*eta*sqrt(norm_s1+ve*n_zk) # norm_s1 is the squared norm of ABDLOP's s1 as per Thm. 5.3 (Fig. 10), hence add ve*d (norm of x)
    sigma2 = gamma2*eta*nu*sqrt(m2*n_zk)

    sigmae = gammae*sqrt(337)*(sqrt(alphae))
    # See [LNP22] fig 10; this is an upper bound on the norm of s^(e). Note that the bound in "public information" is flawed (misplaced square roots)
    print('proof-size log-sigmas:', log(sigma1, 2), log(sigma2, 2), log(sigmae, 2))
    lgq             = (log(q_zk,2))
    challnge_size   = ceil(log(2*kappa+1,2))*n_zk
    hint_size = 2.25*n*n_zk
    # p.49 of LNP22 above the paragraph "Dilithium compression" + we already have ve in m1 + we do not send t^{(d)}'s
    size_plain      = (n+ell+(256/n_zk+1)+1+lamb)*n_zk*lgq + m1*n_zk*log(4.13*sigma1,2) + m2*d*log(4.13*sigma2,2) + 256*log(4.13*sigmae,2) + challnge_size
    # p.50 (top) of LNP22
    size_cut        = n*n_zk*(lgq - D)+(ell+(256/n_zk+1)+1+lamb)*n_zk*lgq + m1*n_zk*log(4.13*sigma1,2) + (m2-n)*d*log(4.13*sigma2,2)+ 256*log(4.13*sigmae,2) + challnge_size + hint_size
    return size_plain/(8*1024), size_cut/(8*1024)
    

# security
zk_sec_SIS = SIS_security(param_zk)
print('SIS hardness:', zk_sec_SIS)
zk_sec_LWE = LWE_security(param_zk)
print('LWE hardness primal:', zk_sec_LWE)


# sizes
proof = proof_size(param_zk)
print('proof size:', proof)
com_x_size = (n_s * row * k) / (8*1024)
print(f"commitment size:{com_x_size:.2f} KB")
print(f"overall non optimized:{com_x_size+proof[0]:.2f} KB")
print(f"overall with CUT:{com_x_size+proof[1]:.2f} KB")


# number of repetitions (rejection sampling)
reps = 2*exp(14/param_zk['gamma1'] + 1/(2*param_zk['gamma1']**2) + 1/(2*param_zk['gamma2']**2) + 1/(2*param_zk['gammae']**2))
print(f"avg nbr of repeats: {reps:.2f}")


# transcript size
transcript1 = n_s * d * log(q_s,2)   # send t
transcript2 = n_s * m_id * log(4.13*s_prime, 2) # signer sends y 
transcript = (transcript1 + transcript2)/(8*1024)
print(f"transcript size: {transcript:.2f} KB")



