import numpy as np
import math as m
import random as r
import matplotlib.pyplot as plt

class Particle:
    def __init__(self, coord, velocity, mass, L): #задаются координаты частицы, ее скорость, масса и ограничение на перемещение
        #self.initial_coord = np.array(coord)
        self.coord = np.array(coord)
        self.velocity = np.array(velocity)
        self.w = np.zeros(3)
        self.mass = mass
        self.L = L
    def Distance(self, other): #расстояние между двумя частицами
        return m.sqrt(sum(delta**2 for delta in self.VectorDistance(other)))
    def VectorDistance(self, other): #координаты вектора от данной частицы до другой
        delta = other.coord - self.coord
        for e in range(3):
            if abs(delta[e])>self.L/2: #частица взаимодействует только с одной частицей из всех ее копий в соседних клетках
                delta[e] -= self.L*abs(delta[e])/delta[e]
        return delta
    def Force(self, other): #сила взаимодействия между двумя частицами
        r = self.Distance(other)
        return -r**(-7)*(2*r**(-6)-1) #используем потенциал Леннарда-Джонса
    def VectorForce(self, other): #координаты вектора силы взаимодействия между двумя частицами
        return self.Force(other)/self.Distance(other)*self.VectorDistance(other)
    def Acceleration(self, other): #ускорение, приобретаемое частицей под действием такой силы
        return self.Force(other)/self.mass
    def VectorAcceleration(self, other): #координаты вектора ускорения
        return self.VectorForce(other)/self.mass
    def Move(self, delta_x): #перемещение частицы на заданный вектор
        new_coord = self.coord + np.array(delta_x)
        for e in range(3):
            if new_coord[e] > self.L:
                new_coord[e] -= self.L
            if new_coord[e] < 0:
                new_coord[e] += self.L
        self.coord = new_coord

N = 4 #число частиц N^3 = 8
L = (N-1)*3+3 #длина стороны куба, в котором находятся частицы
k_B = 1 #постоянная Больцмана
T = 5 #температура
mass = 1
av = m.sqrt(3*k_B*T/mass) #средняя скорость частиц при заданной температуре (из распределения Максвелла)
Gas=[] #создаем список из элементов класса Particle
for i in range(N):
    for j in range(N):
        for k in range(N):
            Gas.append(Particle([float(3*i+1.5),float(3*j+1.5),float(3*k+1.5)],[float(r.gauss(av,av/m.sqrt(3))*r.choice([-1,1])),float(r.gauss(av,av/m.sqrt(3))*r.choice([-1,1])),float(r.gauss(av,av/m.sqrt(3))*r.choice([-1,1]))], mass, L))

t = 0.002
plt.ion()
for step in range (350):
    plt.clf()
    for PARTICLES in Gas: #проходим по всем частицам из списка
        plt.scatter(PARTICLES.coord[0], PARTICLES.coord[1])
        plt.axis([0, L, 0, L])
        PARTICLES.Move(PARTICLES.velocity*t + PARTICLES.w*t**2/2) #двигаем частицу
        tmp = PARTICLES.w
        for particle in Gas: #проходим по всем ДРУГИМ частицам из списка
            if PARTICLES != particle:
                PARTICLES.w += PARTICLES.VectorAcceleration(particle) #добавляем ускорение от силы взаимодействия с текущей частицей
        PARTICLES.velocity += (tmp+PARTICLES.w)*t/2
    plt.draw()
    plt.pause(0.001)
plt.ioff()
plt.show()

"""
P = Particle([0,0,0],[0,0,0],1,6)
Q = Particle([1,0,0],[0,0,0],1,6)
print(P.VectorForce(Q))
"""
