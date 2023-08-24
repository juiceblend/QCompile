from Rings import Zomega, Exact_U, Ztau
from Circuits import FTCircuit
import numpy as np

'''
Implements methods taking in an exact unitary U and outputting an (F,T)-circuit which multiplies to U.
Also implements exact synthesis of weaves
'''


omega = Zomega(0,1,0,0)
tau = Zomega(0,0,1,-1)

T = Exact_U(Zomega(1,0,0,0), Zomega(0,0,0,0), 6)
F = Exact_U(tau, Zomega(1,0,0,0), 0)

def exact_synthesize(U):
    g = U.mu()
    V = U
    Circuit = FTCircuit()
    while g.value() >= 2: 
        min = 0
        min_value = V.j_mu(min).value()
        for j in range(0,10):
            if V.j_mu(j).value() < min_value:
                min = j
                min_value = V.j_mu(min).value()

        V = F*(T**min)*V
        g = V.mu()
        if min != 0:
            Circuit.join("T", 10 - min)
        Circuit.join("F", 1)
    
    # first factor out the phase of omega^k
    for k in range(10):
        if omega**k == V.u:
            Circuit.addPhase(omega**k)
            V = Exact_U(omega**(10-k), Zomega(0,0,0,0), 5 - 2*k) * V
            break
    
    # Am left with a matrix of the form U(1, 0, k) which we can match to T^j
    if V.k%10 != 5:
        Circuit.join("T",(V.k + 5)%10)
    
    return Circuit

def exact_synthesize_weave(U):
    g = U.mu()
    V = U
    Circuit = FTCircuit()
    Remainder = FTCircuit()
    f_parity = True

    while g.value() >= 4: # need to be careful this doesn't break from rounding error? if it does change to Gauss complexity which is always an integer
        
        min = 0
        #print("V = " + str(V))
        for j in range(0,5):
            if (F*(T**(2*j))*V).mu().value() < (F*(T**min)*V).mu().value():
                min = 2*j

        V = F*(T**min)*V
        g = V.mu()
        Circuit.join("T", 10 - min)
        Circuit.join("F", 1)
        f_parity = not f_parity
    
    # print(f_parity)
    # check if |u| = 1
    if V.u*V.u.star() == Zomega(1,0,0,0):
        # factor out phase so that u = 1
        for k in range(10):
            if omega**k == V.u:
                Remainder.addPhase(omega**k)
                V = Exact_U(omega**(10-k), Zomega(0,0,0,0), 5 - 2*k) * V
                break
        
        # Am left with a matrix of the form U(1, 0, k) which we can match to T^j, for j = k-5
        # if j is even and f_parity = True, V is a weave
        if V.k%10 != 5:
            if f_parity:
                if (V.k-5)%2 == 0: 
                    Circuit.join("T", V.k-5)
                    Circuit.addPhase(omega**k)
            else:
                Remainder.join("T", V.k-5)
            
                
    # otherwise |u| = tau
    # We write V = U[w^j t, w^i, k] and find i,j,k
    elif Ztau.fromOmega(V.u*V.u.star()) == Ztau(1, -1):

        # find v
        for i in range(10):
            if omega**i == V.v:
                break

        # find j
        for j in range(10):
            if (omega**(10-j))*V.u == Zomega(0,0,1,-1):
                break
        
        # find k
        k = V.k

        # print((i,j,k))

        # now V = w^j T^(-i-j) F T^(-i-j+k)
        if not f_parity:
            if (i-j)%2 == 0 and (-i-j+k)%2 == 0:
                Circuit.join("T",i-j)
                Circuit.join("F", 1)
                Circuit.join("T",-i-j+k)
                Circuit.addPhase(omega**j)

                return (Circuit, Remainder)
        
        if not f_parity:
            if Circuit.remove_F():
                Remainder.join("F", 1)
                
        Remainder.join("T",i-j)
        Remainder.join("F", 1)
        Remainder.join("T",-i-j+k)
        Remainder.addPhase(omega**j)
            
        
    else:
        print("meep! error: V.u not as expected")
    
    return (Circuit, Remainder)