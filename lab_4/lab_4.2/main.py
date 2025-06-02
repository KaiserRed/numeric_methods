import math
import matplotlib.pyplot as plt

EPS = 1e-10

def y_exact(x):
    return x**3 + x**2 + 2

def y_exact_prime(x):
    return 3*x**2 + 2*x

def p(x):
    return -4*(x**2 + 3)/(x*(x**2 + 6)) if x != 0 else 0

def q(x):
    return 6/(x**2 + 6) if x**2 + 6 != 0 else 0

def shooting_method(a, b, n, alpha, y2b):
    h = (b - a)/n
    x = [a + i*h for i in range(n+1)]
    x[0] = EPS
    
    y0 = y_exact(EPS)
    y0_prime = y_exact_prime(EPS)
    
    s1 = y0_prime
    s2 = y0_prime + 1.0
    
    y1 = integrate(x, n, y0, s1)
    y2 = integrate(x, n, y0, s2)
    
    F1 = y1[-1] - y_prime(b, n, h, y1) - y2b
    F2 = y2[-1] - y_prime(b, n, h, y2) - y2b
    
    s = s1 - F1*(s2 - s1)/(F2 - F1 + EPS)
    count = 0
    for _ in range(10):
        count+=1
        y = integrate(x, n, y0, s)
        F = y[-1] - y_prime(b, n, h, y) - y2b
        print("==", count)
        if abs(F) < 1e-10:
            break
        s1, F1 = s2, F2
        s2, F2 = s, F
        s = s1 - F1*(s2 - s1)/(F2 - F1 + EPS)
    return x, y

def integrate(x, n, y0, y0_prime):
    y = [0.0]*(n+1)
    z = [0.0]*(n+1)
    y[0] = y0
    z[0] = y0_prime
    
    for i in range(n):
        h_step = x[i+1] - x[i]
        
        k1y = h_step * z[i]
        k1z = h_step * (-p(x[i])*z[i] - q(x[i])*y[i])
        
        k2y = h_step * (z[i] + 0.5*k1z)
        k2z = h_step * (-p(x[i]+0.5*h_step)*(z[i]+0.5*k1z) - q(x[i]+0.5*h_step)*(y[i]+0.5*k1y))
        
        k3y = h_step * (z[i] + 0.5*k2z)
        k3z = h_step * (-p(x[i]+0.5*h_step)*(z[i]+0.5*k2z) - q(x[i]+0.5*h_step)*(y[i]+0.5*k2y))
        
        k4y = h_step * (z[i] + k3z)
        k4z = h_step * (-p(x[i]+h_step)*(z[i]+k3z) - q(x[i]+h_step)*(y[i]+k3y))
        
        y[i+1] = y[i] + (k1y + 2*k2y + 2*k3y + k4y)/6
        z[i+1] = z[i] + (k1z + 2*k2z + 2*k3z + k4z)/6
    
    return y

def y_prime(b, n, h, y):
    return (3*y[-1] - 4*y[-2] + y[-3])/(2*h)

def finite_difference_dirichlet(a, b, n, alpha, beta):
    h = (b - a)/n
    x = [a + i*h for i in range(n+1)]
    x[0] = EPS
    
    A = [0.0]*(n+1)
    B = [0.0]*(n+1)
    C = [0.0]*(n+1)
    D = [0.0]*(n+1)
    
    B[0] = 1.0
    C[0] = 0.0
    D[0] = alpha
    
    for i in range(1, n):
        A[i] = 1/h**2 - p(x[i])/(2*h)
        B[i] = -2/h**2 + q(x[i])
        C[i] = 1/h**2 + p(x[i])/(2*h)
        D[i] = 0.0
    
    A[n] = 0.0
    B[n] = 1.0
    D[n] = beta
    
    return x, thomas_algorithm(A, B, C, D)

def thomas_algorithm(a, b, c, d):
    n = len(d)
    cp = [0.0]*n
    dp = [0.0]*n
    
    cp[0] = c[0]/b[0]
    dp[0] = d[0]/b[0]
    
    for i in range(1, n):
        denom = b[i] - a[i]*cp[i-1]
        cp[i] = c[i]/denom
        dp[i] = (d[i] - a[i]*dp[i-1])/denom
    
    y = [0.0]*n
    y[-1] = dp[-1]
    
    for i in range(n-2, -1, -1):
        y[i] = dp[i] - cp[i]*y[i+1]
    
    return y

def runge_romberg(yh, y2h, p):
    n = len(y2h)
    rr = [0.0]*n
    for i in range(n):
        if i % 2 == 0 and (i//2) < len(yh):
            rr[i] = yh[i//2] + (yh[i//2] - y2h[i])/(2**p - 1)
    return rr

def max_abs_error(y1, y2):
    return max(abs(a - b) for a, b in zip(y1, y2))

if __name__ == "__main__":
    a, b = 0.0, 4.0
    alpha = 2.0
    beta = 26.0
    n = 10

    print("Исходные условия:")
    print("x(x²+6)y'' -4(x²+3)y' +6xy = 0")
    print("y'(0) = 0, y(4) - y'(4) = 26")

    x_shoot, y_shoot = shooting_method(a, b, n, alpha, beta)
    y_exact_vals = [y_exact(xi) for xi in x_shoot]
    print("\nМетод стрельбы:")
    for i in range(len(x_shoot)):
        print(f"i = {i:3d}  x={x_shoot[i]:.4f}  y={y_shoot[i]:.6f}  y_exact={y_exact_vals[i]:.6f}  |err|={abs(y_shoot[i]-y_exact_vals[i]):.2e}")
    max_err_shoot = max_abs_error(y_shoot, y_exact_vals)
    print(f"Максимальная ошибка метода стрельбы: {max_err_shoot:.2e}")

    alpha_fd = y_exact(0.0)
    beta_fd = y_exact(4.0)
    x_fd1, y_fd1 = finite_difference_dirichlet(a, b, n, alpha_fd, beta_fd)
    x_fd2, y_fd2 = finite_difference_dirichlet(a, b, 2*n, alpha_fd, beta_fd)
    y_exact_fd = [y_exact(xi) for xi in x_fd1]
    
    print("\nКонечно-разностный метод:")
    for i in range(len(x_fd1)):
        print(f"i = {i:3d}  x={x_fd1[i]:.4f}  y={y_fd1[i]:.6f}  y_exact={y_exact_fd[i]:.6f}  |err|={abs(y_fd1[i]-y_exact_fd[i]):.2e}")
    max_err_fd = max_abs_error(y_fd1, y_exact_fd)
    print(f"Максимальная ошибка КР метода: {max_err_fd:.2e}")

    y_rr = runge_romberg(y_fd1, y_fd2, 2)
    y_exact_rr = [y_exact(xi) for xi in x_fd2]
    max_err_rr = max_abs_error(y_rr, y_exact_rr)
    print(f"\nМаксимальная ошибка Рунге-Ромберга: {max_err_rr:.2e}")


    plt.figure(figsize=(12, 6))
    
    plt.plot(x_shoot, y_exact_vals, 'k-', linewidth=2, label='Точное решение')
    
    plt.plot(x_shoot, y_shoot, 'r--', linewidth=1.5, label='Метод стрельбы')
    
    plt.plot(x_fd1, y_fd1, 'b:', linewidth=1.5, label='Конечно-разностный метод')
    
    plt.xlabel('x', fontsize=12)
    plt.ylabel('y(x)', fontsize=12)
    plt.title('Сравнение методов решения ОДУ', fontsize=14)
    plt.legend(fontsize=12, loc='upper left')
    plt.grid(True, linestyle='--', alpha=0.5)
    
    plt.savefig('methods_comparison.png', dpi=300, bbox_inches='tight')
