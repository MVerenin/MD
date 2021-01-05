import numpy as np
import math as m
import random as r

class Particle:
    def __init__(self, coord, velocity, mass, L):
        self.coord = np.array(coord)
        self.velocity = np.array(velocity)
        self.mass = mass
        self.L = L
    def Distance(self, other):
        return m.sqrt(sum((self.coord[i]-other.coord[i])**2 for i in range(3)))
    def VectorDistance(self, other):
        r = other.coord - self.coord
        for e in range(3):
            if abs(r[e])>self.L/2:
                r[e] -= self.L*abs(r[e])/r[e]
        return r
    def Force(self, other):
        r = self.Distance(other)
        return 8*r**(-7)*(2*r**(-6)-1)
    def VectorForce(self, other):
        return self.Force(other)/self.Distance(other)*self.VectorDistance(other)
    def Acceleration(self, other):
        return self.Force(other)/self.mass
    def VectorAcceleration(self, other):
        return self.VectorForce(other)/self.mass
    def Move(self, delta_x):
        tmp = self.coord + np.array(delta_x)
        for e in range(3):
            if tmp[e] > self.L:
                tmp[e] -= self.L
            if tmp[e] < 0:
                tmp[e] += self.L
        self.coord = tmp

L = 2
N = 2
k_B = 1.38e-23
T = 300
mass = 5
average = m.sqrt(3*k_B*T/mass)
epsilon = 1
Gas=[]
for i in range(N):
    for j in range(N):
        for k in range(N):
            Gas.append(Particle([float(i),float(j),float(k)],[float(r.gauss(1,1)),float(r.gauss(1,1)),float(r.gauss(1,1))], mass, L))

print(Gas[6].Distance(Gas[7]))            
for t in range (10):
    for PARTICLE in Gas:
        w = np.zeros(3)
        v = np.zeros(3)
        for particle in Gas:
            if PARTICLE != particle:
                w += PARTICLE.VectorAcceleration(particle)
        v += w*0.5
        PARTICLE.Move(v*0.5+0.125*w)
        #if PARTICLE == Gas[6]:
            #print(v*0.5+0.125*w)
        #print(Gas[6].coord)
        #print(PARTICLE.Distance(Gas[7]))
for i in range(len(Gas)):
    print(Gas[i].coord)

