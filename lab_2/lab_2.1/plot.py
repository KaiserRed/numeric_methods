import numpy as np
import matplotlib.pyplot as plt

def f1(x):
    return np.log(x + 2)

def f2(x):
    return x**4 - 0.5

x = np.linspace(0, 1.5, 500)
y1 = f1(x)
y2 = f2(x)

plt.figure(figsize=(10, 6))
plt.plot(x, y1, label='ln(x+2)', color='blue')
plt.plot(x, y2, label='x⁴ - 0.5', color='green')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Пересечение функций ln(x+2) и x⁴ - 0.5')
plt.grid()
plt.legend()



plt.axvline(x=0.5, color='gray', linestyle='--', alpha=0.5)
plt.axvline(x=1.2, color='gray', linestyle='--', alpha=0.5)

plt.savefig('intersection_plot.png', dpi=100)
print("График сохранен в intersection_plot.png")
