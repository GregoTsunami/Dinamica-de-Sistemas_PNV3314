import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Parâmetros do sistema (ajustável)
M = 200.0   # (kg)
L = 1.5     # (m)
k_T = 5000.0 # (N·m/rad)
k = 10000.0  # (N/m)
c = 300.0    # (N·s/m)
F2 = 500.0  # (N)
m = 75.0    # (kg)
F1 = 1000.0  # (N)

###############################
#SIMULAÇÃO DO SISTEMA ACOPLADO#
###############################

def sistema(t, y):
    theta, theta_ponto, x, x_ponto = y

    # Equações de movimento
    dtheta = theta_ponto
    dtheta_ponto = -(k_T * theta + k * theta * L + c * theta_ponto * L - F2 * L / 2) / (M*L**2 / 12)
    dx = x_ponto
    dx_ponto = (F1 - 2 * k * x - 2 * c * x_ponto) / m

    return [dtheta, dtheta_ponto, dx, dx_ponto]

# Condições iniciais (ajustável) 
y0 = [0.1, 0.0, 0.1, 0.0]  # pequenas perturbações iniciais
#[theta, theta_ponto, x, x_ponto]

# Tempo da Simulação
t_span = (0, 10)
t_eval = np.linspace(t_span[0], t_span[1], 1000)

# Solução do Sistema
sol = solve_ivp(sistema, t_span, y0, t_eval=t_eval)

# Resultado Gráfico do Deslocamento do Sistema
plt.figure(figsize=(12, 8))


plt.subplot(2, 1, 1)
plt.plot(sol.t, sol.y[0], label='theta (rad)')
plt.plot(sol.t, sol.y[1], label='theta_ponto (rad/s)')
plt.xlabel('Tempo (s)')
plt.ylabel('theta e theta_ponto')
plt.legend()


plt.subplot(2, 1, 2)
plt.plot(sol.t, sol.y[2], label='x (m)')
plt.plot(sol.t, sol.y[3], label='x_ponto (m/s)')
plt.xlabel('Tempo (s)')
plt.ylabel('x e x_ponto')
plt.legend()

plt.tight_layout()
plt.show()



############################################################################
def transferencia(t, y):
    theta_ponto = y

    # Equações de movimento
    dtheta = theta_ponto    
    Xm = 1/(-(dtheta**2)*m + 2*c*dtheta+2*k)
    XM = 1/((-100/3)*(dtheta**2)+100*dtheta+2200)
    return [Xm, XM]

# Solução da Função
t = sol.t
dtheta = sol.y[0]
Xm = sol.y[0]
XM = sol.y[1]

# Resultado Gráfico da Função de Transferência
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

ax1.plot(dtheta, Xm)
ax1.set_title('Xm x Velocidade Angular')
ax1.set_ylabel('Xm')
ax1.grid(True)

ax2.plot(dtheta, XM)
ax2.set_title('XM vs. Velocidade Angular')
ax2.set_xlabel('Velocidade Angular (rad/s)')
ax2.set_ylabel('XM')
ax2.grid(True)

plt.tight_layout()
plt.show()