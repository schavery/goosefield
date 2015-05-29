# -*- coding: utf-8 -*-

import numpy as np
import cmath as cm
import math
from sets import Set
import warnings

import sys

debug = True

def main():
	# p = [1,0,1,1,0,1,1,1,0,1]
	# q = [1,0,0,1,1,0]
	p = [7,3,-2,4]
	q = [2]

	poly_mult_fft(p, q)


def poly_mult_fft(p, q):
	# steps.
	# 1. make sure the polynomials are the same length
	# 2. figure out how many points we need to generate (2d+2)
	# 3. get the nth roots of unity
	# 4. evaluate both p and q at all of the points
	# 4a. this is actually the FFT
	# 4b. split the polynomial into even and odd parts
	# 4c. recurse
	# 5. multiply the results together
	# 5a. now we have a vector of values of pq
	# 5b. Mc = VP
	# 5c. c = M^-1 VP
	# 6. make the M matrix using the roots of unity
	# 7. take the conjugate, divide by number of points
	# 8. matrix mult M^-1 and VP
	# 9. return result

	# global debug

	# step 1 ################
	pn2 = 2**next_power_of_2(len(p))
	qn2 = 2**next_power_of_2(len(q))

	# we need n-1 degree in both p and q,
	# extend with 0s
	n = max(pn2,qn2)

	# make q's length a power of 2
	while len(q) < n:
		q.append(0)

	while len(p) < n:
		p.append(0)

	lenp = len(p)
	lenq = len(q)

	# if debug:
	# 	print 'p = ' + str(p)
	# 	print 'q = ' + str(q)

	# step 2 ################

	# is this wrong?
	# points = 2*lenq
	points = lenq

	# step 3 ################

	# fancy class for roots
	# we can share a single object for both
	# polynomials here, because we have ensured that
	# they are the same length.
	roots = Nthroots(points)

	# if debug:
	# 	print 'Primitive root of unity = ' + str(roots)

	# step 4 ################

	# step 4 is recursive, so its a function
	VP = fft(p, roots)
	VQ = fft(q, roots)

	# if debug:
	# 	print VP
	# 	print VQ

	# sys.exit()

	values = [VP[x]*VQ[x] for x in xrange(lenp)]

	# if debug:
	# 	print "pq values:"
	# 	print values

	# sys.exit()

	# step 5 ################
	# get the inverse matrix of omegas
	M = [[0 for x in xrange(lenp)] for x in xrange(lenp)]

	sqroot = roots.sq()
	for x in xrange(lenp):
		for xx in xrange(lenp):
			M[x][xx] = sqroot.kth(x*xx)
		M[x] = np.conjugate(M[x])

	# if debug:
	# 	print M

	# sys.exit()

	c = np.dot(M, values)
	warnings.simplefilter("ignore")
	print [round(math.copysign(1, c[x]) * abs(c[x]/lenp), 2) for x in xrange(len(c))]

def next_power_of_2(int):
	n = 1
	while 2**n < int:
		n += 1

	return n


def fft(p, omega):
	if len(p) == 1:
		return p
	else:
		global debug

		# first, split p into evens and odds
		pe, po = [], []
		for idx, val in enumerate(p):
			if idx % 2 == 0:
				pe.append(val)
			else:
				po.append(val)

		# if debug:
		# 	print 'Evens: ' + str(pe)
		# 	print 'Odds: ' + str(po)

		# the next iteration needs w^2
		omega_squared = omega.sq()

		# if debug:
		# 	print 'Next root of unity: ' + str(omega_squared)

		U = fft(pe, omega_squared)
		W = fft(po, omega_squared)

		V = [0 for x in xrange(len(p))]

		for x in xrange(len(p) / 2):
			V[x] = U[x] + omega_squared.kth(x)*W[x]
			V[x + len(p) / 2] = U[x] - omega_squared.kth(x)*W[x]

		# if debug:
		# 	print V

		return V


class Nthroots(object):
	"""store and work on the values of the nth roots of unity"""

	def __init__(self, n):
		self.n = n

	def w(self):
		return math.pi / self.n

	def kth(self, k):
		theta = k * math.pi / self.n
		return cm.cos(theta) + cm.sin(theta) * 1j

	def sq(self):
		# doubles the angle omega
		# returns new instance
		new_n = self.n / 2
		return self.__class__(new_n)

	def __str__(self):
		return 'π/' + str(self.n)

	def __repr__(self):
		return repr(self.__str__())


if __name__ == "__main__":
	main()
	# n = Nthroots(16)
	# print n.kth(0)
