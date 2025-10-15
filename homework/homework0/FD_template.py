"""
Numerical Mathematics for Engineering II WS 25/26
Homework 00 Exercise 0.4
Finite difference quotients
"""

import matplotlib.pyplot as plt
import numpy as np

def forwardDiff(f,x,h):
	return  ## question (b): complete this line

def centralDiff(f,x,h):
	return  ## question (b): complete this line

# Sanity checks
def f_lin(x):
	return x

def f_lin_prim(x):
	return 1

def f_quad(x):
	return x**2

def f_quad_prim(x):
	return 2*x
x = 1
h = 0.1
err_FD_lin =   ## question (b): complete this line
err_CD_quad =  ## question (b): complete this line
print(f"Error of (FD) for a linear function: {err_FD_lin:.3e}")
print(f"Error of (CD) for a quadratic function: {err_CD_quad:.3e}")



# Convergence curves 
def f(x):
	return np.exp(x)

def fprim(x):
	return np.exp(x)

# Compute derivatives at x = 1
x = "todo" ## Question (c): complete this line 
fprim_exact = "todo" ## Question (c): complete this line 

h = 10**np.linspace(-1,-5,5)	  ## question (d): comment this line
#h = 10**np.linspace(-1,-12,12)   ## question (d): decomment this line

fprimF = forwardDiff(f,x,h)
fprimC = "todo" ## Question (c): complete this line 

errF = np.abs((fprimF - fprim_exact)/fprim_exact)
errC = np.abs((fprimC - fprim_exact)/fprim_exact)

plt.loglog(h,errF,label='Forward difference')
plt.loglog(h,errC,label='Central difference')
plt.loglog(h,h**1,label='$ch$', linestyle="--", color='tab:blue')
plt.loglog(h,h**2,label='$ch^2$', linestyle="--" , color='tab:orange')
plt.xlabel('Step size')
plt.ylabel('Relative error')
plt.grid()
plt.legend()
#plt.savefig('DQ_coarse_h.png')
#plt.savefig('DQ_fine_h.png')
plt.show()
