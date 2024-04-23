from cProfile import label
import numpy as np
import matplotlib.pyplot as plt


t = []
point = []
with open('data_rho.txt', 'r') as f:
    data = f.readlines()
    
title, *x = data[0].split()
x = [float(el) for el in x]
for row in data[1:]:
    T, *U = row.split()
    U = [float(el) for el in U]
    t.append(float(T))
    point.append(U)

point = np.array(point)

plt.plot(x, point[len(point) // 2])
plt.grid()

print(t[len(point) // 2])
l_x=x[0]
r_x=x[-1]
l_t=t[0]
r_t=t[-1]

figure, axes = plt.subplots()

c = axes.pcolormesh(x, t, point)
axes.set_title('rho')
axes.axis([l_x, r_x, l_t, r_t])
figure.colorbar(c)

plt.show()
