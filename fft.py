import numpy as np
import cmath as cm
import math
from sets import Set
from fractions import Fraction

# debug = True

def main():
	p = [1,0,1,1,0,1,1,1,0,1]
	q = [1,0,0,1,1,0]

	poly_mult(p, q)
	

def poly_mult(p, q):
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


	# make q's length a power of 2
	while len(q) < 16:
		q.append(0)

	while len(p) < 16:
		p.append(0)

	lenp = len(p)
	lenq = len(q)

	# if debug:
	# 	print 'p = ' + str(p)
	# 	print 'q = ' + str(q)

	# step 2 ################
	points = 2*lenq

	# step 3 ################
	
	# fancy class for roots
	roots = Nthroots(points)

	# if debug:
	# 	print 'roots of unity = ' + str(roots.w())

	# step 4 ################

	# step 4 is recursive, so its a function
	VP = fft(lenp, p, roots)
	VQ = fft(lenq, q, roots)

	values = [VP[x]*VQ[x] for x in xrange(lenp)]

	# step 5 ################
	# get the inverse matrix of omegas
	M = [[0 for x in xrange(lenp)] for x in xrange(lenp)]

	for x in xrange(lenp):
		for xx in xrange(lenp):
			M[x][xx] = roots.kth(x*xx)
		M[x] = np.conjugate(M[x])

	# print M

	c = np.dot(M, values)
	print [c[x]/lenp for x in xrange(len(c))]


def fft(lenp, p, omega):
	if lenp == 1:
		return [p[0]]
	else:
		# global debug

		# first, split p into evens and odds
		pe, po = [], []

		for idx, val in enumerate(p):
			if idx % 2 == 0:
				pe.append(val)
			else:
				po.append(val)

		omega_squared = omega.sq()
		halflen = lenp/2

		U = fft(halflen, pe, omega_squared)
		W = fft(halflen, po, omega_squared)

		V = {}
		for x in xrange(halflen):
			V[x] = U[x] + omega.kth(x)*W[x]
			V[x + halflen] = U[x] - omega.kth(x)*W[x]

		return V


class Nthroots(object):
	"""store and work on the values of the nth roots of unity"""

	def __init__(self, n):
		self.n = n
		self.frac = Fraction(1, n)

	def w(self):
		return self.frac

	def kth(self, k):
		theta = k * math.pi * self.frac
		return cm.cos(theta) + cm.sin(theta) * 1j

	def sq(self):
		# doubles the angle omega, stored in frac
		new_n = self.n / 2

		return self.__class__(new_n)


if __name__ == "__main__":
	main()
	# n = Nthroots(16)
	# print n.kth(0)
