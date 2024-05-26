import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


#Função Característica do Movimento do Sistema
def sistema(t, y, M, m, ks, cs, k1, c1, k2, c2):
    yM, yM_o, ym, ym_o = y
    yM_oo = (-2*cs*yM_o  - 2*ks*yM + (k1+k2)*(ym-yM) + (c1+c2)*(ym_o-yM_o)) / M
    ym_oo = (- (c1+c2)*(ym_o-yM_o) - (k1+k2)*(ym-yM)) / m
    return [yM_o, yM_oo, ym_o, ym_oo]

#Parâmetros do Sistema
M = 100.0               #kg
m = 50.0                #kg
ks = 2000.0             #N/m
cs = 100.0              #Ns/m
k1 = 3000.0             #N/m
c1 = 150.0              #Ns/m
k2 = 2500.0             #N/m
c2 = 120.0              #Ns/m

#Condições iniciais [yM(0), yM'(0), ym(0), ym'(0)]
y0 = [0.0, 0.0, 1.0, 0.0] #ym(0) é para deixar o sistema excitado

#Tempo da simulação
t_span = (0, 10) #intervalo de tempo
t_eval = np.linspace(0, 10, 500) #500 pontos de tempo para plotar

#Solução do Sistema
sol = solve_ivp(sistema, t_span, y0, args=(M, m, ks, cs, k1, c1, k2, c2), t_eval=t_eval)

#Resultado Gráfico do Deslocamento do Sistema
plt.plot(sol.t, sol.y[0], label='y_M(t)')
plt.plot(sol.t, sol.y[2], label='y_m(t)')
plt.xlabel('Tempo')
plt.ylabel('Deslocamento')
plt.legend()
plt.grid()
plt.show()