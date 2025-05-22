import numpy as np
import matplotlib.pyplot as plt

a = 2

def f1(x1, x2):
    return a*x1**2 - x1 + x2**2 - 1

def f2(x1, x2):
    return x2 - np.tan(x1)

x1 = np.linspace(-1.5, 1.5, 400)
x2 = np.linspace(-2, 2.5, 400)
X1, X2 = np.meshgrid(x1, x2)

F1 = f1(X1, X2)
F2 = f2(X1, X2)

plt.figure(figsize=(10, 8))

plt.contour(X1, X2, F1, levels=[0], colors='blue', linewidths=2)
plt.contour(X1, X2, F2, levels=[0], colors='red', linewidths=2)
plt.xlabel('x1')
plt.ylabel('x2')
plt.axvline(0, color ='black', linewidth=0.5)
plt.axhline(0, color ='black', linewidth=0.5)
plt.title('Система уравнений:\n(2x1 - x1 + x2² - 1 = 0) и (x2 - tan(x1) = 0)')
plt.grid()
plt.savefig('system_plot.png', dpi=100)
print("График сохранен в system_plot.png")