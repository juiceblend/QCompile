from Rings import Zomega
import numpy as np
import math
import cmath

'''
Define classes for (F,T)-circuits and braids
Annoying literature convention: (F,T) circuits multiply the wrong way, e.g. (FTT) is equal to T*T*F
'''

class FTCircuit:
    def __init__(self):
        self.circuit = []
        self.phase = Zomega(1,0,0,0)
    
    def join(self, M, pow):
        self.circuit = [(M, pow%10)] + self.circuit
    
    def append(self, pair):
        self.circuit.append(pair)
    
    def remove_F(self):
        if len(self.circuit) == 0:
            return False
        
        if self.circuit[0][0] == "F":
            self.circuit = self.circuit[1:]
            return True

        return False
        
        

    def addPhase(self, phase):
        self.phase = self.phase * phase
    
    def __add__(self, other):
        result = FTCircuit()
        for pair in self.circuit:
            result.append(pair)
        for pair in other.circuit:
            result.append(pair)
        
        result.addPhase(self.phase)
        result.addPhase(other.phase)

        return result
    
    def toBraid(self):
        braid = Braid()
        braid.addPhase(cmath.phase(self.phase.value()))

        l = len(self.circuit)
        i = 0

        while i < l:
            if self.circuit[l-1-i][0] == "T":
                pow = self.circuit[l-1-i][1]
                braid.addPhase(2*pow*math.pi/5)
                braid.add(1,3*pow)

                i = i+1

            elif self.circuit[l-1-i][0] == "F":
                
                pow = self.circuit[l-1-i][1]

                # case which shouldn't happen
                if pow != 1:
                    braid.addPhase(4*pow*math.pi/5)
                    for j in range(pow):
                        braid.add(1,1)
                        braid.add(2,1)
                        braid.add(1,1)
                    
                    i = i+1
                    print("wrong power of F")
                    continue


                # check if a power of sigma_2 can be implemented
                if l-3-i >= 0:
                    if self.circuit[l-3-i][0] == "F" and self.circuit[l-2-i][0] == "T":
                        pow = self.circuit[l-2-i][1]
                        braid.addPhase(2*pow*math.pi/5)
                        braid.add(2,3*pow)

                        i = i+3
                        continue
                
                braid.addPhase(4*math.pi/5)
                braid.add(1,1)
                braid.add(2,1)
                braid.add(1,1)
                i = i+1

        return braid
    
    # implemented into toBraid
    '''
    def toWeave(self):
        braid = Braid()
        braid.addPhase(cmath.phase(self.phase.value()))

        l = len(self.circuit)

        for i in range(int(l/4)):
            pow = self.circuit[l-1-4*i][1]
            braid.addPhase(2*pow*math.pi/5)
            braid.add(1,3*pow)

            sandwich_pow = self.circuit[l-3-4*i][1]
            braid.addPhase(2*sandwich_pow*math.pi/5)
            braid.add(2,3*sandwich_pow)
        
        if l%4 == 1:
            pow = self.circuit[0][1]
            braid.addPhase(2*pow*math.pi/5)
            braid.add(1,3*pow)
    
        return braid
    '''

    
    def numerical_approximate(self):
        T = np.array([
            [1, 0],
            [0, cmath.exp(0.2j * math.pi)]
        ])

        tau = 0.5 * (math.sqrt(5) - 1)

        F = np.array([
            [tau, math.sqrt(tau)],
            [math.sqrt(tau), -tau]
        ])

        result = np.identity(2)
        l = len(self.circuit)
        for i in range(l):
            if self.circuit[l-1-i][0] == "T":
                result = np.matmul(result, np.linalg.matrix_power(T, self.circuit[l-1-i][1]))

            elif self.circuit[l-1-i][0] == "F":
                result = np.matmul(result, np.linalg.matrix_power(F, self.circuit[l-1-i][1]))
        
        return self.phase.value() * result

    
    def __str__(self):
        s = str(self.phase)
        for elt in self.circuit:
            if elt[1] == 1:
                s = s + "({0})".format(elt[0])
            else:
                s = s + "({0}^{1})".format(elt[0], elt[1])
        
        return s


class Braid:
    def __init__(self):
        self.braid = []
        self.phase = 0
    
    def add(self, sigma, pow):
        if len(self.braid) == 0:
            if pow%10 != 0:
                self.braid.append((sigma, pow%10))
            return
        
        if self.braid[-1][0] == sigma:
            prev_pow = self.braid.pop(-1)[1]
            if (prev_pow + pow)%10 != 0:
                self.braid.append((sigma, (prev_pow+pow)%10)) # since sigma1 and sigma2 both have order 10
        else:
            if pow%10 != 0:
                self.braid.append((sigma, pow%10))

    def addPhase(self, phase):
        self.phase = self.phase + phase
        self.phase = self.phase - (self.phase//(2*math.pi))*2*math.pi

    def __add__(self, other):
        result = Braid()
        for pair in self.braid:
            result.braid.append(pair)
        for pair in other.braid:
            result.braid.append(pair)
        
        result.addPhase(self.phase)
        result.addPhase(other.phase)

        return result
    
    def numerical_approximate(self):
        sigma1 = np.array([
            [cmath.exp(1.2j * math.pi), 0],
            [0, cmath.exp(0.6j * math.pi)]
        ])

        tau = 0.5 * (math.sqrt(5) - 1)

        F = np.array([
            [tau, math.sqrt(tau)],
            [math.sqrt(tau), -tau]
        ])

        sigma2 = np.matmul(F, np.matmul(sigma1, F))

        bases = [sigma1, sigma2]

        result = np.identity(2)
        for elt in self.braid:
            result = np.matmul(result, np.linalg.matrix_power(bases[elt[0]-1], elt[1]))
        
        return cmath.exp(1j * self.phase) * result

    
    def __str__(self):
        '''
        s = "exp({0}\u03C0i)".format(self.phase/math.pi) 
        for elt in self.braid:
            if elt[1] == 1:
                if elt[0] == 1:
                    s = s + "(\u03C3\u2081)"
                elif elt[0] == 2:
                    s = s + "(\u03C3\u2082)"
            else:
                if elt[0] == 1:
                    s = s + "(\u03C3\u2081^{0})".format(elt[1])
                elif elt[0] == 2:
                    s = s + "(\u03C3\u2082^{0})".format(elt[1])      
        
        return s
        '''
        s = "exp({0}\u03C0i)".format(self.phase/math.pi) 
        for elt in self.braid:
            if elt[0] == 1:
                for i in range(elt[1]):
                    s = s + "1"
            else:
                for i in range(elt[1]):
                    s = s + "2"   
        
        return s
    
    def evaluate_braid(braid):
        sigma1 = np.array([
            [cmath.exp(1.2j * math.pi), 0],
            [0, cmath.exp(0.6j * math.pi)]
        ])

        tau = 0.5 * (math.sqrt(5) - 1)

        F = np.array([
            [tau, math.sqrt(tau)],
            [math.sqrt(tau), -tau]
        ])

        sigma2 = np.matmul(F, np.matmul(sigma1, F))

        bases = [sigma1, sigma2]

        result = np.identity(2)
        for elt in braid:
            result = np.matmul(result, bases[int(elt)-1])
        
        return result