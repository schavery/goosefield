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
	q = [1,1]

	# poly_mult_fft(p, list(reversed(q)))
	poly_mult_fft(p, q)
	print
	poly_mult_naive(p, q)

def poly_mult_naive(p, q):
	v = [0 for x in xrange(len(p)+len(q))]

	for idx,x in enumerate(p):
		for idy,y in enumerate(q):
			v[idx+idy] += x*y

	print v



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

	# Figure out the degree of the resulting polynomial
	dp = poly_degree(p)
	dq = poly_degree(q)

	result_degree = dp+dq

	# so the number of points we need to evaluate
	points = result_degree + 1

	# to make the recursion simple, make the lengths
	# go to a power of 2

	next2 = 2**next_power_of_2(points)

	while len(q) < next2:
		q.append(0)

	while len(p) < next2:
		p.append(0)

	lenp = len(p)

	roots = Nthroots(next2)

	VP = fft(p, roots)
	VQ = fft(q, roots)


	values = [VP[x]*VQ[x] for x in xrange(lenp)]


	# get the inverse matrix of omegas
	M = [[0 for x in xrange(lenp)] for x in xrange(lenp)]

	for x in xrange(lenp):
		for xx in xrange(lenp):
			M[x][xx] = roots.kth(x * xx)
		M[x] = np.conjugate(M[x])


	c = np.dot(M, values)


	# there is a warning about the copysign call dropping the
	# imaginary part of the complex number, so I hide it.
	warnings.simplefilter("ignore")
	print [round(math.copysign(1, c[x]) * abs(c[x]/lenp), 2) for x in xrange(len(c))]

def next_power_of_2(int):
	n = 1
	while 2**n < int:
		n += 1

	return n


def poly_degree(p_array):
	for x in xrange(len(p_array), 0, -1):
		if p_array[x-1] != 0:
			return x-1


def fft(p, omega):
	if len(p) == 1:
		return p
	else:
		# first, split p into evens and odds
		pe, po = [], []
		for idx, val in enumerate(p):
			if idx % 2 == 0:
				pe.append(val)
			else:
				po.append(val)

		# the next iteration needs w^2
		omega_squared = omega.sq()

		U = fft(pe, omega_squared)
		W = fft(po, omega_squared)

		V = [0 for x in xrange(len(p))]

		for x in xrange(len(p) / 2):
			V[x] = U[x] + omega.kth(x)*W[x]
			V[x + len(p) / 2] = U[x] - omega.kth(x)*W[x]

		return V


class Nthroots(object):
	"""store and work on the values of the nth roots of unity"""

	def __init__(self, n):
		self.n = n

	def w(self):
		return 2 * math.pi / self.n

	def kth(self, k):
		theta = 2 * k * math.pi / self.n
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
