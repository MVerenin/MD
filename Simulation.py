import numpy as np
import math as m
import random as r
import matplotlib.pyplot as plt

class Particle:
    def __init__(self, coord, velocity, mass, L): #задаются координаты частицы, ее скорость, масса и ограничение на перемещение
        self.initial_coord = np.array(coord)
        self.coord = np.array(coord)
        self.velocity = np.array(velocity)
        self.mass = mass
        self.L = L
    def Distance(self, other): #расстояние между двумя частицами
        return m.sqrt(sum((self.coord[i]-other.coord[i])**2 for i in range(3)))
    def VectorDistance(self, other): #координаты вектора от данной частицы до другой
        delta = other.coord - self.coord
        for e in range(3):
            if abs(delta[e])>self.L/2: #частица взаимодействует только с одной частицей из всех ее копий в соседних клетках
                delta[e] -= self.L*abs(delta[e])/delta[e]
        return delta
    def Force(self, other): #сила взаимодействия между двумя частицами
        r = self.Distance(other)
        return -24*r**(-7)*(2*r**(-6)-1) #используем потенциал Леннарда-Джонса
    def VectorForce(self, other): #координаты вектора силы взаимодействия между двумя частицами
        return self.Force(other)/self.Distance(other)*self.VectorDistance(other)
    def Acceleration(self, other): #ускорение, приобретаемое частицей под действием такой силы
        return self.Force(other)/self.mass
    def VectorAcceleration(self, other): #координаты вектора ускорения
        return self.VectorForce(other)/self.mass
    def Move(self, delta_x): #перемещение частицы на заданный вектор
        new_coord = self.coord + np.array(delta_x)
        for e in range(3):
            if new_coord[e] > 6:
                new_coord[e] -= self.L
            if new_coord[e] < 0:
                new_coord[e] += self.L
        self.coord = new_coord

d = {1.5: 'green', 4.5: 'blue'}
L = 6 #длина стороны куба, в котором находятся частицы
N = 2 #число частиц N^3 = 8
k_B = 1380 #постоянная Больцмана
T = 300 #температура
mass = 0.05 #масса частицы (масса молекулы азота)
av = 1
#av = m.sqrt(3*k_B*T/mass) #средняя скорость частиц при заданной температуре (из распределения Максвелла)
Gas=[] #создаем список из элементов класса Particle
for i in range(N):
    for j in range(N):
        for k in range(N):
            Gas.append(Particle([float(3*i+1.5),float(3*j+1.5),float(3*k+1.5)],[float(r.gauss(av,av/m.sqrt(3))*r.choice([-1,1])),float(r.gauss(av,av/m.sqrt(3))*r.choice([-1,1])),float(r.gauss(av,av/m.sqrt(3))*r.choice([-1,1]))], mass, L))

for i in range(len(Gas)): #Скорости получаются порядка 500 м/c, что примерно соответствует реальности
    print(Gas[i].velocity)

plt.ion()
for t in range (250):
    plt.clf()
    for PARTICLES in Gas: #проходим по всем частицам из списка
        plt.scatter(PARTICLES.coord[0], PARTICLES.coord[1],color = d[PARTICLES.initial_coord[2]])
        plt.axis([0, 6, 0, 6])
        w = np.zeros(3) #здесь будут складываться ускорения
        #v = PARTICLES.velocity #начальная скорость
        for particle in Gas: #проходим по всем ДРУГИМ частицам из списка
            if PARTICLES != particle:
                w += PARTICLES.VectorAcceleration(particle) #добавляем ускорение от силы взаимодействия с текущей частицей
        #v += w*0.005 #изменяем скорость
        PARTICLES.velocity += 0.01*w
        PARTICLES.Move(PARTICLES.velocity*0.005) #двигаем частицу
    plt.draw()
    plt.pause(0.001)
plt.ioff()
plt.show()
for i in range(len(Gas)):
    print(Gas[i].coord)
"""
P = Particle([0,0,0],[0,0,0],5,5)
Q = Particle([2,0,0],[0,0,0],5,5)
print(P.Acceleration(Q))
"""
