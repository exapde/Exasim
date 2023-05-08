import numpy as np
import matplotlib.pyplot as plt

n1=62
n2=52
elem_idx1 = np.arange(n1)
elem_h1 = np.zeros(n1)
elem_idx2 = np.arange(n2)
elem_h2 = np.zeros(n2)
c1=1e-6
b1=1.05
b2 = 1.08

elem_h1 = c1*b1**elem_idx1
c2 = elem_h1[-1]
elem_h2 = c2*b2**elem_idx2

elem_size=np.hstack((elem_h1, elem_h2[1:]))

# Compare to exponential series (can't use integral because not compounded continuously):
# https://en.wikipedia.org/wiki/Geometric_series
elem_h_sum1 = c1*(1-b1**(elem_idx1+1))/(1-b1)   # Have to do elems+1 beacuse the formula expects the n-1 th term
elem_h_sum2 = c2*(1-b2**(elem_idx2+1))/(1-b2)   # Have to do elems+1 beacuse the formula expects the n-1 th term
elem_h_sum2 -= elem_h_sum2[0]   # Subtracting offset
elem_h_sum2 += elem_h_sum1[-1]  # Shifting over to start at the end of the first set

elem_h_sum = np.hstack((elem_h_sum1, elem_h_sum2[1:]))

slope1 = (elem_h1[-1]-elem_h1[0])/(elem_h_sum1[-1]-elem_h_sum1[0])
x_growth1 = np.linspace(0,elem_h_sum1[-1])+1e-6
y_int1 = 1e-6*(1-slope1)
y_growth1 = slope1*x_growth1+y_int1

slope2 = (elem_h2[-1]-elem_h2[0])/(elem_h_sum2[-1]-elem_h_sum2[0])
x_growth2 = np.linspace(elem_h_sum1[-1],elem_h_sum2[-1])
y_int2 = elem_h2[0] -slope2*elem_h_sum2[0]
y_growth2 = y_int2 + slope2*x_growth2      # Point-slope form

x_combined = np.linspace(0,14e-3,1000)

alpha = 0.5*(np.tanh(100000*(x_combined-elem_h_sum2[0]))+1)
y_combined = (1-alpha)*(slope1*x_combined+y_int1) + alpha*(slope2*x_combined+y_int2)

print('elem_h_sum2[0]', elem_h_sum2[0])
print('slope1', slope1)
print('y_int1', y_int1)
print('slope2', slope2)
print('y_int2', y_int2)

plt.scatter(elem_h_sum1*1e3, elem_h1*1e3)
plt.scatter(elem_h_sum2*1e3, elem_h2*1e3)
plt.plot(x_growth1*1e3, (slope1*x_growth1+y_int1)*1e3)
plt.plot(x_growth2*1e3, (y_growth2)*1e3)
plt.plot(x_combined*1e3, y_combined*1e3, 'red', linewidth=3, linestyle='--')

plt.axhline(y=0, c='black')
plt.axvline(x=0, c='black')
plt.title('Mesh size progression from needle tip, Chen 2017')
plt.xlabel('Distance from wall [mm]')
plt.ylabel('Element size h [mm]')
plt.show()
