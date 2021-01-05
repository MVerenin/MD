import numpy as np
import math as m
import random as r

class Particle:
    def __init__(self, coord, velocity, mass, L): #задаются координаты частицы, ее скорость, масса и ограничение на перемещение
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
        return 8*r**(-7)*(2*r**(-6)-1) #используем потенциал Леннарда-Джонса
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

L = 2 #длина стороны куба, в котором находятся частицы
N = 2 #число частиц N^3 = 8
k_B = 1.38e-23 #постоянная Больцмана
T = 300 #температура
mass = 4.67e-26 #масса частицы (масса молекулы азота)
av = m.sqrt(3*k_B*T/mass) #средняя скорость частиц при заданной температуре (из распределения Максвелла)
Gas=[] #создаем список из элементов класса Particle
for i in range(N):
    for j in range(N):
        for k in range(N):
            Gas.append(Particle([float(i),float(j),float(k)],[float(r.gauss(av,av/m.sqrt(3))),float(r.gauss(av,av/m.sqrt(3))),float(r.gauss(av,av/m.sqrt(3)))], mass, L))
"""
for i in range(len(Gas)): Скорости получаются порядка 500 м/c, что примерно соответствует реальности
    print(Gas[i].velocity)
"""        
for t in range (1):
    for PARTICLE in Gas: #проходим по всем частицам из списка
        w = np.zeros(3) #здесь будут складываться ускорения
        v = PARTICLE.velocity #начальная скорость
        for particle in Gas: #проходим по всем ДРУГИМ частицам из списка
            if PARTICLE != particle:
                w += PARTICLE.VectorAcceleration(particle) #добавляем ускорение от силы взаимодействия с текущей частицей
        v += w*0.05 #изменяем скорость
        PARTICLE.Move(v*0.05) #двигаем частицу
for i in range(len(Gas)):
    print(Gas[i].coord)

