import math
import sympy # to check primality
from Rings import Ztau, Zomega

'''
Implements method solve_norm_equation(xi) which takes in some xi in Z[tau] and outputs
an x in Z[omega] which solves the norm equation |x|^2 = xi
'''

def legendre(a, p):
    return pow(a, (p - 1) // 2, p)

def modular_sqrt(n, p):
    # returns the square root of n modulo p
    assert legendre(n, p) == 1, "not a square (mod p)"
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        return pow(n, (p + 1) // 4, p)
    for z in range(2, p):
        if p - 1 == legendre(z, p):
            break
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r

def splitting_root(xi):
    # Find a square root of tau-2 modulo xi
    p = xi.N()
    b1 = pow(xi.b, -1, p)
    return modular_sqrt(-xi.a*b1 - 2, p)


def unit_dlog(u):
    # Takes a unit u in Ztau and returns (s,k) where u = s*tau^k
    s, k, a, b = 1, 0, u.a, u.b

    if a < 0:
        a, b, s = -a, -b, -s

    mu = a*b

    # repeatedly multiply by tau or tau^-1 until the remainder has |mu| <= 1
    while abs(mu) > 1:
        if mu > 1:
            a, b, k = b, a-b, k-1
        else:
            a, b, k = a+b, a, k+1
        mu = a*b

    if a == -1 and b == 0:
        s = -s
    elif a == 0 and b == 1:
        k = k+1
    elif a == 0 and b == -1:
        s, k = -s, k+1
    elif a == 1 and b == -1:
        k = k+2
    elif a == -1 and b == 1:
        s, k = -s, k+2
    elif a == 1 and b == 1:
        k = k-1
    elif a == -1 and b == -1:
        s, k = -s, k-1

    return (s, k)


def easy_factor(xi):
    # factors out the integer part as well as (2 - tau) if it exists
    c = math.gcd(xi.a, xi.b)
    a1, b1 = int(xi.a/c), int(xi.b/c)
    xi_1 = Ztau(a1, b1)

    factors = []

    if math.sqrt(c).is_integer():
        if c > 1:
            factors = [(Ztau(int(math.sqrt(c)),0),2)]
    else:
        if math.sqrt(c/5).is_integer():
            factors = [(Ztau(int(math.sqrt(c/5)), 0),2), (Ztau(5,0),1)]
        else:
            factors = [(xi, 1)]
            return factors # equation won't be solveable
    
    n = xi_1.N()
    if n%5 == 0:
        p1 = Ztau(2, -1)
        xi_2 = xi_1//p1
        factors.append((p1, 1))
        factors.append((xi_2, 1))
    else:
        factors.append((xi_1, 1))
    
    return factors


def easy_solveable(factors):
    # Checks whether a list of factors easily solves the norm equation
    for i in range(len(factors)):
        (xi, k) = factors[i]
        if k%2 == 1:
            if xi.a != 5 or xi.b != 0:
                p = xi.N()
                r = p%5
                if (not sympy.isprime(p)) or (r != 0 and r != 1):
                    return False
    
    return True


def gcd(a,b):
    if a.N() == 0:
        return b
    if b.N() == 0:
        return a
    
    if a.N() <= b.N():
        q = b//a
        r = b - q*a
        # print(b, a, q, r, r.N())
        return gcd(a,r)
    else:
        q = a//b
        r = a - q*b
        # print(a, b, q, r, r.N())
        return gcd(b,r)
    

def solve_norm_equation(xi):
    # Takes in as input xi an element of Z[tau] and
    # outputs x from Zomega which solves |x|^2 = xi

    # first check solveability
    if xi.value() < 0 or xi.dot().value() < 0:
        return Zomega(0,0,0,0)
    factors = easy_factor(xi)
    if not easy_solveable(factors):
        return Zomega(0,0,0,0)
    
    x = Zomega(1,0,0,0)

    # solve for each prime factor separately, then multiply solutions together
    for i in range(len(factors)):
        (xi_i, m) = factors[i]
        x = x * (Zomega.fromTau(xi_i)**(int(m/2)))
        if m%2 == 1:
            if xi_i == Ztau(5,0):
                x = x * (Zomega.fromTau(Ztau(1, 2)))
            elif xi_i == Ztau(2, -1):
                x = x*Zomega(-1,2,-1,1)
            else:
                M = splitting_root(xi_i)
                y = gcd(Zomega.fromTau(xi_i), Zomega(M,0,0,0) - Zomega(-1,2,-1,1))
                u = xi_i // y.Ni()
                (s,m) = unit_dlog(u)
                x = x * (Zomega.fromTau(Ztau(0,1)**int(m/2))) * y
    
    return x        

