# Exercise 0.1

# Numpy

```cpp
for (auto i = 0 ; i < N ; i++) {
    w[k] = u[k] + v[k];
}
```

```python
w = np.dot(u,v)
```

Numpy is much simpler and faster to use. Thus, whenever possible, use numpy rather than developing your own.

```python
import numpy as np
a = np.array([1.0.0])
...
# creating arrays
a = np.zeros((3,4))
a.flatten()
a.reshape(6,2)
a.T
**2
```

## Dot Product vs Elementwise Product

```python
A = np.array([[1.1][0.1]])
B = np.array([[2,0].[3,4]])
A * B # elementwise
A @ B # matrix product
A.dot(B) # matrix product
```

## Indexing

```python
# indexed, sliced, iterated ... all possible
b[2,3]
b[0:5.1]
b[:,1]
```

# Matplotlib Package

Visualize results

```
import matplotlib.pyplot as plt
x = np.linspace(0,2*np.pi)
plt.plot(x, np.sin(x).pi)
plt.xlabel("$x$") # latex ok
plt.ylabel("$\\sin(x)$")
plt.show()
```

# Exercise 0.2

## Big O Notation

- x-> infinitive : 5
- x-> 0 : 3

# Convergence_template.py

```
   h     |   E_T    |     p_T    |   E_S    |    p_S
7.85e-01 | 1.25e-02 |     --     | 5.50e-04 | --
3.93e-01 | 3.11e-03 | 2.0113e+00 | 3.25e-05 | 4.0824e+00
1.96e-01 | 7.76e-04 | 2.0028e+00 | 2.00e-06 | 4.0200e+00
9.82e-02 | 1.94e-04 | 2.0007e+00 | 1.25e-07 | 4.0050e+00
4.91e-02 | 4.85e-05 | 2.0002e+00 | 7.79e-09 | 4.0012e+00

```

## 수렴률

- Trapezoid Rule : p_T == 2.0 -> O(h^2)
- Simpsons Rule : p_S == 4.0 => O(h^4)

## Error

As h /2

- Trapezoid -> 1/4 (h^2)
- Simpson -> 1/16 (h^4)

-> Matches our expectation

### Which rule is more accurate? Can you easily read the convergence order on the curve? Can you easily read the convergence order on the curve?
