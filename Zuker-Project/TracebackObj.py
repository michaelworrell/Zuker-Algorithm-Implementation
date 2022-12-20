"""
Class that stores traceback information, including the optimal pairing case,
which alongside other class variables can be used to redirect
to another traceback object. 
"""
import sys


class Traceback:
	"""
	Initializes the Traceback object.
	NOTE: zukers.py noes not access object values until functions
	independently setting those values have already been called.
	"""
	def __init__(self):
		self.iprime = 0
		self.jprime = 0
		self.k = 0
		self.case = 0

	"""
	Sets the optimal case for this object
	(a positive integer between 1 and 4 inclusive).
	"""
	def setCase(self, case):
		self.case = case
		return

	"""
	Sets the iPrime and jPrime for this object
	(a positive integer for both iprime and jprime).
	"""
	def setijprime(self, iprime, jprime):
		self.iprime, self.jprime = iprime, jprime
		return

	"""
	Sets k for this traceback object.
	"""
	def setK(self, k):
		self.k = k
		return
