import sys
import numpy as np
import random
from numpy import save
import time
import multiprocessing
import glob
import os


class space:
    def __init__(self, x_min, x_max, y_min, y_max, z_min, z_max, N):
        #inicializace prostoru, mezní podmínky
        print('Space init')
        self.x_min = x_min 
        self.x_max = x_max 
        self.y_min = y_min 
        self.y_max = y_max 
        self.z_min = z_min 
        self.z_max = z_max
        self.N = N
        self.particle_posicion_generator = setuper(x_min = self.x_min, x_max = self.x_max, y_min = self.y_min, y_max = self.y_max, z_min = self.z_min, z_max = self.z_max)
        self.getcells((self.x_max, self.y_max, self.z_max))
        self.populatecells(self.N)
        self.find_centroids()

    def getcells(self, D_cells):
        #rozdělení prostoru do buněk
        print('getcells', end = '\t')
        start = time.perf_counter()
        self.cells = []
        for x in range(D_cells[0]):
            cellsx = []
            for y in range(D_cells[1]):
                cellsy = []
                for z in range(D_cells[2]):
                    cellsy.append(cell(x = x, y = y, z = z))
                cellsx.append(cellsy)
            self.cells.append(cellsx)
        end = time.perf_counter()
        print(f'done in {end - start}')

    def populatecells(self, N):
        #Vygenerování N částic a jejich přidělení do patřičné buňky
        print('populatecells', end = '\t')
        start = time.perf_counter()
        points = np.array([])
        for _ in range(N):
            c = self.particle_posicion_generator.coordinates()
            self.cells[int(c[0])][int(c[1])][int(c[2])].particles.append(particle(coordinates = np.array([c[0], c[1], c[2]]), F = np.array([0.0, 0.0, 0.0]), v = np.array([0.0, 0.0, 0.0]), m = 18.0))
            points = np.append(points, c, axis=0)
        points = np.resize(points, (len(points), 3))
        
        #uložení počátečního stavu
        save('Out/points_000', points)
        end = time.perf_counter()
        print(f'done in {end - start}')

    def find_centroids(self):
        #nalezení všech těžišť
        print('find_centroids', end = '\t')
        start = time.perf_counter()
        for cellx in self.cells:
            for celly in cellx:
                for cellz in celly:
                    cellz.find_centroid()
        end = time.perf_counter()
        print(f'done in {end - start}')

    def update(self, dt, t):
        #krok simulace
        print(f'Space update, time = {t}, ', end = '\t')
        start = time.perf_counter()

        #zjisteni sil z potencialu
        for cellx in self.cells:
            for celly in cellx:
                #process_list = []
                #for x in range(0, len(celly), 12):
                #    for n in range(12):
                #        v = n + x
                #        t = multiprocessing.Process(target=celly[v].get_forces_in_cell, args=[self])
                #        t.start()
                #        process_list.append(t)
                #    for process in process_list:
                #        process.join()
                    


                for cellz in celly:
                    cellz.get_forces_in_cell(self)

        #aplikace sil
        for cellx in self.cells:
            for celly in cellx:
                for cellz in celly:
                    cellz.aplly_forces(dt)

        #kontarola prislustnosti castic
        for cellx in self.cells:
            for celly in cellx:
                for cellz in celly:
                    cellz.check_grid(self)

        #aktualizace tezist
        for cellx in self.cells:
            for celly in cellx:
                for cellz in celly:
                    cellz.find_centroid()

        end = time.perf_counter()
        print(f'done in {end - start}')

    def Save(self):
        #ulozeni pozic
        n = len(glob.glob('Out/points_*'))
        print(f'Saving\t{n}')
        if n < 10:
            n = '00' + str(n)
        elif n < 100:
            n = '0' + str(n)
        else:
            n = str(n)

        points = np.array([])
        for cellx in self.cells:
            for celly in cellx:
                for cellz in celly:
                    for particle in cellz.particles:
                        points = np.append(points, particle.coordinates, axis=0)
        points = np.resize(points, (len(points), 3))
        save('Out/points_' + n, points)

class cell:
    def __init__(self, x, y, z):
        #deklarace prostoru pro částice, polohy buňky a jejího středu
        self.particles = []
        self.x = x
        self.y = y
        self.z = z
        self.c = np.array([x + 0.5, y + 0.5, z + 0.5])
        self.epsilon = 0.0103
        self.sigma = 3.4

    def find_centroid(self):
        #nalezení těžiště buňky a celkové hmotnosti
        self.N_particles = len(self.particles)
        if self.N_particles > 0:
            self.m = 0
            self.centroid = np.array([0.0, 0.0, 0.0])        
            for particle in self.particles:
                self.centroid += particle.coordinates*particle.m
                self.m += particle.m
            self.centroid = self.centroid/self.m
        else:
            self.centroid = np.array([0.0, 0.0, 0.0])
            self.m = 0

    def get_forces_in_cell(self, space):
        #nalezeni vnitrich a okolnich sil
        for particle in self.particles:
            #interakce castice s ostatnimi v bunce
            for other_particle in (other_particle for other_particle in self.particles if particle is not other_particle):
                r = np.linalg.norm(particle.coordinates - other_particle.coordinates)
                if r > 0.25:
                    di = -(particle.coordinates - other_particle.coordinates) / r
                    particle.F = particle.F + self.LJ_force(r = r)*di
                else:
                    particle.v = particle.v*0.9

            #interakce castice s ostatnimi bunkami
            for cellx in space.cells:
                for celly in cellx:
                    for cellz in (cellz for cellz in celly if cellz is not self):
                        r = np.linalg.norm(particle.coordinates - cellz.centroid)
                        di = -(particle.coordinates - cellz.centroid) / r
                        particle.F = particle.F + self.LJ_force(r = r)*di*cellz.m

    def aplly_forces(self, dt):
        #aplikace sil 
        for particle in self.particles:
            particle.update(dt)

    def check_grid(self, space):
        #kontrola prislusnosti castic
        for particle in self.particles:
            if not self.x == int(particle.coordinates[0]) or not self.y == int(particle.coordinates[1]) or not self.z == int(particle.coordinates[2]):
                if int(particle.coordinates[0]) >= space.x_max or int(particle.coordinates[0]) < space.x_min or int(particle.coordinates[1]) >= space.y_max or int(particle.coordinates[1]) < space.y_min or int(particle.coordinates[2]) >= space.z_max or int(particle.coordinates[2]) < space.z_min:
                    self.particles.remove(particle)
                else:
                    space.cells[int(particle.coordinates[0])][int(particle.coordinates[1])][int(particle.coordinates[2])].particles.append(particle)
                    self.particles.remove(particle)

    def LJ_force(self, r):
        #
        return 48 * self.epsilon * np.power(self.sigma, 12) / np.power(r, 13) - 24 * self.epsilon * np.power(self.sigma, 6) / np.power(r, 7)


class particle:
    def __init__(self, coordinates, F, v, m):
        self.coordinates = coordinates
        self.F = F
        self.v = v
        self.m = m

    def update(self, dt):
        dp = dt*self.F
        self.v += + dp/self.m
        self.coordinates += self.v*dt
        self.F = np.array([0.0, 0.0, 0.0])

class setuper:
    def __init__(self, x_min, x_max, y_min, y_max, z_min, z_max):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.z_min = z_min
        self.z_max = z_max

    def x(self):
        return self.x_min + (self.x_max-self.x_min)*random.random()

    def y(self):
        return self.y_min + (self.y_max-self.y_min)*random.random()

    def z(self):
        return self.z_min + (self.z_max-self.z_min)*random.random()

    def si(self):
        return 1 if random.random() < 0.5 else -1

    def coordinates(self): #valec
        X = self.x()*self.si()/2
        Y = self.y_min + (np.sqrt(((self.x_max-self.x_min)/2)**2-X**2)-self.y_min)*random.random()*self.si()
        return np.array([X+2.5, Y+2.5, self.z()])

def main():
    #
    start = time.perf_counter()
    Space = space(x_min = 0, x_max = 5, y_min = 0, y_max = 5, z_min = 0, z_max = 15, N = 2250)
    end = time.perf_counter()
    print('\n', f'Done in {end - start}')

    T = 10000; dt = 0.00000001; t = 0
    while(T >= t):
        Space.update(dt = dt, t = t*dt)
        t = t + 1
        if t % 50 == 0:
            Space.Save()

def Check():
    for f in glob.glob('Out/po*'):
        os.remove(f)

def Help():
    #
    print('[-u, -l]')

def Options(opt):
    #
    if len(opt) > 0:
        if opt[0] == '-h':
            Help()
        elif opt[0] == '-u':
            main()
        elif opt[0] == '-d':
            pass
            #clear()
        sys.exit()
    else:
        Help()
        sys.exit()

if __name__ == "__main__":
    #
    Check()
    Options(sys.argv[1:])