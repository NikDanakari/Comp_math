import numpy as np
import matplotlib.pyplot as plt


t = []
point = []
with open('data-LW-1.01.txt', 'r') as f:
    data = f.readlines()
    
T, step = data[0].split()
T_max = float(T)
step = float(step)

for row in data[1:]:
    T, *U = row.split()
    U = [float(el) for el in U]
    t.append(float(T))
    point.append(U)

point = np.array(point)

l_x = 0
N = len(point[0])
x = np.linspace(0, step * (N - 1), N)

figure, axes = plt.subplots()

c = axes.pcolormesh(x, t, point)
axes.set_title('LW-1.01')
axes.axis([x[0], x[-1], 0, T_max])
figure.colorbar(c)

plt.show()


