###JC69.py

from math import *
def JC69(m, t):
	a = 0
	b = 1
	matrixq = [[a,b,b,b],[b,a,b,b],[b,b,a,b],[b,b,b,a]]
	x = 1/4.0 + 3/4.0 * exp(-4*m*t)
	y = 1/4.0 - 1/4.0 * exp(-4*m*t)
	a = x
	b = y
	matrixq = [[a,b,b,b],[b,a,b,b],[b,b,a,b],[b,b,b,a]]
	print(matrixq)
	
JC69(1, 1)
m = 1
t = 1
x = 1/4.0 + 3/4.0 * exp(-4*m*t)

print(x)