import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

#Parâmetros do Sistema (ajustável)
M = 100.0               #kg
m = 50.0                #kg
ks = 2000.0             #N/m
cs = 100.0              #Ns/m
k1 = 3000.0             #N/m
c1 = 150.0              #Ns/m
k2 = 2500.0             #N/m
c2 = 120.0              #Ns/m
N = M + m               #kg
kT = 1000.00            #N*m/rad
R = 2.0                 #m

###############################
#SIMULAÇÃO DO SISTEMA ACOPLADO#
###############################

#Função Característica do Movimento do Sistema
def sistema(t, y, M, m, ks, cs, k1, c1, k2, c2, kT, R, N):
    ym, ym_o, yM, yM_o, xN, xN_o = y

    # Equações de movimento
    dym = ym_o
    dvm = -(c1 + c2)*(ym_o - yM_o)/m - (k1 + k2)*(ym - yM)/m
    dyM = yM_o
    dvM = ((c1 + c2)*(ym_o - yM_o) + (k1 + k2)*(ym - yM) - 2*cs*yM_o - 2*ks*yM) / M
    dxN = xN_o
    dvxN = -kT * xN / (N * R**2)

    return [dym, dvm, dyM, dvM, dxN, dvxN]

# Condições iniciais (ajustável) 
y0 = [0.1, 0.0, 0.1, 0.0, 0.1, 0.0]  #excitação inicial
     #[ym, ym_o, yM, yM_o, xN, xM_o]

# Tempo da Simulação
t_span = (0, 10)
t_eval = np.linspace(t_span[0], t_span[1], 500)

# Solução do Sistema
sol = solve_ivp(sistema, t_span, y0, args=(M, m, ks, cs, k1, c1, k2, c2, kT, R, N), t_eval=t_eval)

# Resultado Gráfico do Deslocamento do Sistema
plt.figure(figsize=(12, 8))
plt.plot(sol.t, sol.y[0], label='Deslocamento de y_m')
plt.plot(sol.t, sol.y[2], label='Deslocamento de y_M')
plt.xlabel('Tempo (s)')
plt.ylabel('Deslocamento (m)')
plt.legend()
plt.grid()
plt.show()


plt.figure(figsize=(12, 8))
plt.plot(sol.t, sol.y[4], label='Deslocamento de x_N')
plt.plot(sol.t, sol.y[5], label='Velocidade de x_N')
plt.xlabel('Tempo (s)')
plt.ylabel('Deslocamento e Velocidade em x (m, m/s)')
plt.legend()
plt.grid()
plt.show()