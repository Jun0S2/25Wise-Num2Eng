"""
Numerical Mathematics for Engineering II WS 25/26
Homework 00 Exercise 0.3
Revision on convergence curves and tables 
"""

# Import necessary libraries
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

def f(x):
    return 2 + np.cos(x)

a = 0
b = np.pi/2

def approximate_integral(n,rule):
    h = (b-a)/n
    x = np.linspace(a, b, n+1)

    if(rule == "trapezoid"):
        integral = ## add
    elif(rule == "simpson"):
        integral = scipy.integrate.simpson(f(x), x)

    relative_error = ## add
    return relative_error 




if __name__ == '__main__':
    N = 5
    error_t = np.zeros(N)
    error_s = np.zeros(N)
    for k in range(N):
        error_t[k] = approximate_integral(2**(k+1), "trapezoid")
        error_s[k] = approximate_integral(2**(k+1), "simpson")

    h = (b-a) / (2**(np.arange(N)+1))

    plt.plot(h, ## add , label="trapezoid") 
    plt.plot(h, ## add , label="Simpson")
    plt.xlabel("$h$")
    plt.ylabel("relative error")
    plt.legend()
    plt.show()

    plt.loglog(h, error_t, "x-", label="trapezoid")
    plt.loglog(h, error_s, "x-", label="Simpson")
    plt.loglog(## add label="$h^2$")  
    plt.loglog(## add label="$h^4$")  
    plt.xlabel("$h$")
    plt.ylabel("relative error")
    plt.legend()
    plt.show()
    
    ### Convergence tables
    rate_t = ## add
    rate_s = ## add
    print("")
    print(f"   h     |   E_T    |     p_T    |   E_S    |    p_S    ")
    print(f"{h[0]:.2e} | {error_t[0]:.2e} |     --     | {error_s[0]:.2e} | -- ")
    for i in range(1,N):
        print(f"{h[i]:.2e} | {error_t[i]:.2e} | {rate_t[i-1]:.4e} | {error_s[i]:.2e} | {rate_s[i-1]:.4e}")
    print("")
