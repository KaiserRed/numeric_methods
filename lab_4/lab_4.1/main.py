import matplotlib.pyplot as plt

def exact_solution(x):
    return x**2 + x + 1/x

def f_system(x, state):
    """Преобразование ОДУ 2-го порядка в систему 1-го порядка"""
    y, dy = state
    ddy = (3*x**2 + y - x*dy) / (x**2)
    return [dy, ddy]

def euler_system(F, x0, state0, h, n):
    """Векторный метод Эйлера"""
    x = [x0]
    states = [state0.copy()]
    
    for _ in range(n):
        current_state = states[-1]
        derivatives = F(x[-1], current_state)
        new_state = [current_state[i] + h * derivatives[i] for i in range(len(current_state))]
        x.append(x[-1] + h)
        states.append(new_state)
    
    return x, states

def rk4_system(F, x0, state0, h, n):
    """Векторный метод Рунге-Кутты 4-го порядка"""
    x = [x0]
    states = [state0.copy()]
    
    for _ in range(n):
        xi = x[-1]
        yi = states[-1]
        
        k1 = F(xi, yi)
        k2 = F(xi + h/2, [yi[j] + h/2 * k1[j] for j in range(len(yi))])
        k3 = F(xi + h/2, [yi[j] + h/2 * k2[j] for j in range(len(yi))])
        k4 = F(xi + h, [yi[j] + h * k3[j] for j in range(len(yi))])
        
        new_state = [
            yi[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])
            for j in range(len(yi))
        ]
        
        x.append(xi + h)
        states.append(new_state)
    
    return x, states

def adams4_system(F, x0, state0, h, n):
    """Векторный метод Адамса 4-го порядка (предиктор-корректор)"""
    if n < 4:
        return rk4_system(F, x0, state0, h, n)
    
    x, states = rk4_system(F, x0, state0, h, 3)
    
    derivatives_history = [F(xi, si) for xi, si in zip(x, states)]
    # Предиктор
    for i in range(3, n):
        predictor = [
            states[i][j] + h/24*(
                55*derivatives_history[i][j] 
                - 59*derivatives_history[i-1][j] 
                + 37*derivatives_history[i-2][j] 
                - 9*derivatives_history[i-3][j]
            )
            for j in range(len(state0))
        ]
        
        # Корректор (Адамс-Моултон)
        f_pred = F(x[i] + h, predictor)
        corrector = [
            states[i][j] + h/24*(
                9*f_pred[j] 
                + 19*derivatives_history[i][j] 
                - 5*derivatives_history[i-1][j] 
                + derivatives_history[i-2][j]
            )
            for j in range(len(state0))
        ]
        
        x.append(x[i] + h)
        states.append(corrector)
        derivatives_history.append(F(x[-1], states[-1]))
    
    return x, states

def runge_romberg(Fh, Fkh, k, p):
    error = abs(Fkh - Fh) / (k**p - 1)
    return error

def print_results_table(x_euler, y_euler, x_rk, y_rk, x_adams, y_adams, exact):
    headers = ["i", "Метод", "x", "y", "f(x)", "|y - f(x)|"]
    data = []
    
    for i in range(len(x_euler)):
        # Эйлер
        data.append([
            i,
            "Эйлер",
            f"{x_euler[i]:.6f}",
            f"{y_euler[i]:.6f}",
            f"{exact[i]:.6f}",
            f"{abs(y_euler[i] - exact[i]):.2e}"
        ])
        
        # Рунге-Кутта
        if i < len(x_rk):
            data.append([
                i,
                "Рунге-Кутта",
                f"{x_rk[i]:.6f}",
                f"{y_rk[i]:.6f}",
                f"{exact[i]:.6f}",
                f"{abs(y_rk[i] - exact[i]):.2e}"
            ])
            
        if i < len(x_adams):
            data.append([
                i,
                "Адамс",
                f"{x_adams[i]:.6f}",
                f"{y_adams[i]:.6f}",
                f"{exact[i]:.6f}",
                f"{abs(y_adams[i] - exact[i]):.2e}"
            ])

    print("\nРезультаты вычислений:")
    print(f"{'':-^90}")
    print(f"{headers[0]:^12} | {headers[1]:<12} | {headers[2]:^5} | {headers[3]:^9} | {headers[4]:^9} | {headers[5]:^12}")
    print(f"{'':-^90}")
    
    for row in data:
        print(f"{row[0]:^3} | {row[1]:<12} | {row[2]:^5} | {row[3]:^9} | {row[4]:^9} | {row[5]:^12}")
    
    print(f"{'':-^90}")


print("Дифференциальное уравнение: x²y'' + xy' - y - 3x² = 0")
print("Начальные условия: \n y(1) = 3 \n y'(1) = 2")
print("Интервал интегрирования:vx ∈ [1, 2]")
print("Шаг сетки: h = 0.1")
print("Точное решение: y(x) = x² + x + 1/x ")

x0 = 1.0
y0 = 3.0
dy0 = 2.0
h = 0.01
n = 15

x_euler, states_euler = euler_system(f_system, x0, [y0, dy0], h, n)
x_rk, states_rk = rk4_system(f_system, x0, [y0, dy0], h, n)
x_adams, states_adams = adams4_system(f_system, x0, [y0, dy0], h, n)

y_euler = [state[0] for state in states_euler]
y_rk = [state[0] for state in states_rk]
y_adams = [state[0] for state in states_adams]

x_exact = [x0 + i*h for i in range(n+1)]
y_exact = [exact_solution(x) for x in x_exact]

print_results_table(x_euler, y_euler, x_rk, y_rk, x_adams, y_adams, y_exact)

h2 = h / 2
n2 = int((2.0 - 1.0) / h2)

_, states_euler_h2 = euler_system(f_system, x0, [y0, dy0], h2, n2)
_, states_rk_h2 = rk4_system(f_system, x0, [y0, dy0], h2, n2)
_, states_adams_h2 = adams4_system(f_system, x0, [y0, dy0], h2, n2)

y_euler_h2 = [s[0] for s in states_euler_h2]
y_rk_h2 = [s[0] for s in states_rk_h2]
y_adams_h2 = [s[0] for s in states_adams_h2]

exact_end = exact_solution(x0 + n*h)

rr_euler = runge_romberg(y_euler[-1], y_euler_h2[-1], 2, 1)
rr_rk = runge_romberg(y_rk[-1], y_rk_h2[-1], 2, 4)
rr_adams = runge_romberg(y_adams[-1], y_adams_h2[-1], 2, 4)

# Вывод таблицы
print("\nОценка погрешности методом Рунге-Ромберга:")
print("+------------+----------------+----------------+----------------+")
print("| Метод      | Значение (h)   | Значение (h/2) | Погрешность RR |")
print("+------------+----------------+----------------+----------------+")
print(f"| Эйлер      | {y_euler[-1]:>12.6f} | {y_euler_h2[-1]:>12.6f} | {rr_euler:>12.6e} |")
print(f"| Рунге-Кутты| {y_rk[-1]:>12.6f} | {y_rk_h2[-1]:>12.6f} | {rr_rk:>12.6e} |")
print(f"| Адамс      | {y_adams[-1]:>12.6f} | {y_adams_h2[-1]:>12.6f} | {rr_adams:>12.6e} |")
print("+------------+----------------+----------------+----------------+")
print(f"| Точное     | {exact_end:>12.6f} | {'':>12} | {'':>12} |")
print("+------------+----------------+----------------+----------------+")

plt.figure(figsize=(12, 6))
plt.plot(x_exact, y_exact, 'k-', label='Точное решение')
plt.plot(x_euler, y_euler, 'r--', label='Эйлер')
plt.plot(x_rk, y_rk, 'g-.', label='Рунге-Кутта 4')
plt.plot(x_adams, y_adams, 'b:', label='Адамс 4')
plt.xlabel('x')
plt.ylabel('y(x)')
plt.legend()
plt.grid(True)
plt.title('Сравнение методов (векторный подход)')
plt.savefig('vector_solution_comparison.png', dpi=300)

