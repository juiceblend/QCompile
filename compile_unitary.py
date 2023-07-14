from Rings import Zomega, Ztau, Exact_U
from norm_equation import easy_factor, easy_solveable, solve_norm_equation
from exact_synthesize import exact_synthesize, exact_synthesize_weave
from Circuits import FTCircuit, Braid
from math import sqrt, ceil, pi, log, sin, cos, tan, acos
import cmath
import random
import numpy as np
from toqito.random import random_unitary

'''
Implements method compile_unitary, which outputs a braid sequence approximating an arbitary
input unitary U to precision epsilon using the trace distance
'''

# -----------------------------------------
# Constants

phi = 0.5 * (sqrt(5) + 1)
tau = 0.5 * (sqrt(5) - 1)

def Rz(angle):
    return np.array([
        [cmath.exp(-0.5j*angle), 0],
        [0, cmath.exp(0.5j*angle)]
    ])

def Pauli_X():
    return np.array([
        [0, 1],
        [1, 0]
    ])

def F():
    return np.array([
        [tau, sqrt(tau)],
        [sqrt(tau), -1*tau]
    ])



# -----------------------------------------
# approx_real returns an element of Z[tau] which approximates the input x 
# to precision tau^{n-1}(1 - tau^n)

FibArray = [1, 1]
 
def Fib(n):
    if n<0:
        print("Incorrect input")
    elif n<= len(FibArray):
        return FibArray[n-1]
    else:
        temp_fib = Fib(n-1)+Fib(n-2)
        FibArray.append(temp_fib)
        return temp_fib    
    

def approx_real(x,n):
    p, q = Fib(n), Fib(n+1)
    u, v = (-1)**(n+1)*Fib(n), (-1)**n * Fib(n-1)
    c = round(x*q)
    a = c*v + p * round(c*u/q)
    b = c*u - q * round(c*u/q)

    return Ztau(a, b)


# -----------------------------------------
# random_sample generates an element of Z[omega] approximating exp(-i theta/2) to within precision epsilon

def random_sample(theta, epsilon, r):
    # must have 0 < theta < pi/5
    # no idea what r is honestly but needs to be >= 1

    C = sqrt(phi/(4*r))
    m = ceil(log(C*epsilon*r, tau)) + 1
    N = ceil(phi**m)

    y_min = r*phi**m *(sin(theta) - epsilon*(sqrt(4-epsilon**2)*cos(theta) + epsilon*sin(theta))/2)
    y_max = r*phi**m *(sin(theta) + epsilon*(sqrt(4-epsilon**2)*cos(theta) - epsilon*sin(theta))/2)
    x_max = r*phi**m *((1 - epsilon**2/2)*cos(theta) - epsilon*sqrt(1 - epsilon**2/4)*sin(theta))
    x_c = x_max - r*epsilon**2*phi**m/(4*cos(theta))

    j = random.randint(1, N-1)
    y = y_min + j*(y_max - y_min)/N
    a_y = approx_real(y/sqrt(2-tau), m)
    x = x_c - (a_y.value()*sqrt(2-tau) - y_min)*tan(theta)
    a_x = approx_real(x,m)

    return Zomega.fromTau(a_x) + Zomega(-1,2,-1,1)*Zomega.fromTau(a_y)



# -----------------------------------------------------------------------------
# Methods for compiling Rz(theta) and Rz(theta)X

def compile_R_rotation(angle, epsilon):
    C = sqrt(phi/4)
    m = ceil(log(C*epsilon, tau)) + 1
    theta = 0
    for k in range(10):       
        adjusted_angle = (-angle/2 - pi*k/5) - ((-angle/2 - pi*k/5)//(2*pi)) * (2*pi)
        if adjusted_angle >= 0 and adjusted_angle <= pi/5:
            theta = adjusted_angle
            break
    
    not_found, u, v = True, 0, 0

    while not_found:
        u0 = random_sample(theta, epsilon, 1)
        phi_t = Ztau(1, 1)
        xi = phi_t*(phi_t**(2*m) - Ztau.fromOmega(u0*u0.star()))
        factors = easy_factor(xi)
        if easy_solveable(factors):
            not_found = False
            u = Zomega(0,1,0,0)**k * Zomega.fromTau(Ztau(0,1))**m * u0
            v = Zomega.fromTau(Ztau(0,1))**m * solve_norm_equation(xi)
    
    # Debug block
    #print("u, v = " + str(u) + ", " + str(v))
    #print("distance = " + str(sqrt(1 - 0.5 * abs((u.value()*cmath.exp(1j*angle/2) + (u.star().value()*cmath.exp(-1j*angle/2)))))))
    
    Circuit = exact_synthesize(Exact_U(u, v, 0))

    return Circuit

def compile_RX_rotation(angle, epsilon):
    r = sqrt(phi)
    C = sqrt(phi/(4*r))
    m = ceil(log(C*epsilon*r, tau)) + 1
    theta = 0
    for k in range(10):
        adjusted_angle = (angle/2 + pi/2 - pi*k/5) - ((angle/2 + pi/2 - pi*k/5)//(2*pi)) * (2*pi)
        if adjusted_angle >= 0 and adjusted_angle <= pi/5:
            theta = adjusted_angle
            break
    
    not_found, u, v = True, 0, 0

    while not_found:
        u0 = random_sample(theta, epsilon, r)
        tau_t = Ztau(0,1)
        phi_t = Ztau(1,1)
        xi = phi_t**(2*m) - tau_t * Ztau.fromOmega(u0*u0.star())
        factors = easy_factor(xi)
        if easy_solveable(factors):
            not_found = False
            v = Zomega(0,1,0,0)**k * Zomega.fromTau(Ztau(0,1))**m * u0
            u = Zomega.fromTau(Ztau(0,1))**m * solve_norm_equation(xi)
    
    v = Zomega(-1,0,0,0) * v.star()
    ex = Exact_U(u, v, 0)
    # print(ex.value())
    Circuit = exact_synthesize(ex)

    return Circuit


# -----------------------------------------------------------------------
# Methods for synthesizing an arbitrary unitary U

# decompose U into a product of rotation matrices
# Specifically outputs (alpha, beta, gamma, delta, X) such that
# U = e^i*delta Rz(alpha)*F*Rz(beta)*F*Rz(gamma)     if X is True
# U = e^i*delta Rz(alpha)*F*Rz(beta)*F*Rz(gamma) X   if X is False

def decompose_unitary(U):
    X = False
    a, b, c, d = U[0,0], U[0,1], U[1,0], U[1,1]
    if abs(a) < 0.5:
        X = True
        a, b, c, d = U[0,1], U[0,0], U[1,1], U[1,0]
    
    delta = 0.5 * cmath.phase(a*d - b*c)

    a, b, c, d = cmath.exp(-1j * delta)*a, cmath.exp(-1j * delta)*b, cmath.exp(-1j * delta)*c, cmath.exp(-1j * delta)*d, 
    
    beta = acos(1 - abs(c)**2/(2*tau**3))
    alpha = -(cmath.phase(a) - cmath.phase(cmath.exp(1j*beta) + tau)
              - cmath.phase(c) + cmath.phase(1 - cmath.exp(1j*beta)))
    gamma = -2*(cmath.phase(a) - cmath.phase(cmath.exp(1j*beta) + tau)) - alpha - beta
    
    return (alpha, beta, gamma, delta, X)


# ---------------------------------------------------------------------------------------
# Compiles an arbitrary unitary to precision epsilon

def compile_unitary(U, epsilon):
    epsilon = epsilon**6
    # break U into rotations and then compile them individually
    (alpha, beta, gamma, delta, X) = decompose_unitary(U)
    A = compile_R_rotation(alpha, epsilon)
    B = compile_R_rotation(beta, epsilon)
    if X:
        C = compile_RX_rotation(gamma, epsilon)
    else:
        C = compile_R_rotation(gamma, epsilon)


    # debug block
    # A_num, B_num, C_num = A.numerical_approximate(), B.numerical_approximate(), C.numerical_approximate()
    # estimate = cmath.exp(delta*pi*1j) * np.matmul(A_num, np.matmul(F(), np.matmul(B_num, np.matmul(F(), C_num))))
    # print(abs(estimate[0]))

    F_circ = FTCircuit()
    F_circ.join("F", 1)

    A_braid, B_braid, C_braid, F_braid = A.toBraid(), B.toBraid(), C.toBraid(), F_circ.toBraid()

    total_braid = A_braid + F_braid + B_braid + F_braid + C_braid
    total_braid.addPhase(delta)

    return total_braid