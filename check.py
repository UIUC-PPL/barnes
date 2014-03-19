#!/usr/bin/env python
import math
import sys

#calculates the nearest k neighbor ball radius for each particle read from the input file

k=10
if len(sys.argv) > 1: k = int(sys.argv[1])
print 'Nearest k: ', k
print 'Particle \t Ball radius'

#read the particle file and fill the particles array
particles = []
with open("100_particles", "r") as particles_file:
	for particle in particles_file:
		sp = particle.split(" ")
		p = []
		for s in sp:
			p.append(float(s))
		particles.append(p)

#function to calculate distance^2 between two particles
def distance(p1, p2):
	return ((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)

for i in range(0,len(particles)):
	distances = []
	for po in particles:
		#calculate all of the distances from a particle and store it in distances
		distances.append(distance(particles[i], po))
	distances.sort()
	print str(i), " \t ", distances[k] 
	#print str(i),"\t",particles[i][0],"\t",particles[i][1],"\t",particles[i][2],"\t",distances[k] 

