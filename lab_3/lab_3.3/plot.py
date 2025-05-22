import matplotlib.pyplot as plt
import numpy as np

# Данные из таблицы
x = np.array([0.1, 0.5, 0.9, 1.3, 1.7, 2.1])
y = np.array([10.0, 2.0, 1.1111, 0.76923, 0.58824, 0.47618])

# Аппроксимирующие многочлены (коэффициенты из Go-кода)
# P1(x) = 6.5919 - 3.7283 * x
# P2(x) = 10.0988 - 14.1074 * x + 4.7178 * x^2

def P1(x): return 6.5919 + (-3.7283) * x
def P2(x): return 10.0988 + (-14.1074) * x + 4.7178 * x**2

# Создаём значения x для построения графика
x_vals = np.linspace(0.05, 2.15, 500)
y_p1 = P1(x_vals)
y_p2 = P2(x_vals)

# Построение графика
plt.figure(figsize=(8, 6))
plt.plot(x, y, 'ro', label='Табличные точки')
plt.plot(x_vals, y_p1, 'b-', label='Линейная аппроксимация P1(x)')
plt.plot(x_vals, y_p2, 'g--', label='Квадратичная аппроксимация P2(x)')

plt.title("Аппроксимация табличной функции")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.legend()
plt.savefig("approximation_plot.png")  # сохранение файла
plt.close()
