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
def sistema(t, y, M, m, ks, cs, k1, c1, k2, c2):
    yM, yM_o, ym, ym_o = y
    yM_oo = (-2*cs*yM_o  - 2*ks*yM + (k1+k2)*(ym-yM) + (c1+c2)*(ym_o-yM_o)) / M
    ym_oo = (- (c1+c2)*(ym_o-yM_o) - (k1+k2)*(ym-yM)) / m
    return [yM_o, yM_oo, ym_o, ym_oo]

#Condições iniciais (ajustável) [yM(0), yM'(0), ym(0), ym'(0)]
y0 = [0.0, 0.0, 1.0, 0.0] #ym(0) = 1 é para deixar o sistema excitado

#Tempo da simulação (ajustável)
t_span = (0, 10) #Intervalo de tempo
t_eval = np.linspace(0, 10, 500) #500 pontos de tempo para plotar

#Solução do Sistema
sol = solve_ivp(sistema, t_span, y0, args=(M, m, ks, cs, k1, c1, k2, c2), t_eval=t_eval)

#Resultado Gráfico do Deslocamento do Sistema
plt.plot(sol.t, sol.y[0], label='y_M(t)')
plt.plot(sol.t, sol.y[2], label='y_m(t)')
plt.xlabel('Tempo[s]')
plt.ylabel('Deslocamento[m]')
plt.legend()
plt.title('Sistema Dinâmico Acoplado')
plt.grid()
plt.show()



##################################
#SIMULAÇÃO DO SISTEMA DESACOPLADO#
##################################

#Função para M
def sistema_M(t, y, M, ks, cs):
    yM, yM_o = y
    yM_oo = (-2*cs*yM_o - 2*ks*yM) / M
    return [yM_o, yM_oo]

#Função para m
def sistema_m(t, y, m, k1, c1, k2, c2):
    ym, ym_o = y
    ym_oo = (- (c1+c2)*ym_o - (k1+k2)*ym) / m
    return [ym_o, ym_oo]

#Condições Iniciais para M (ajustável) [yM(0), yM'(0)]
y0_M = [-1.0, 0.0] #yM(0) = -1 é para deixar o sistema excitado

#Condições Iniciais para m (ajustável) [ym(0), ym'(0)]
y0_m = [1.0, 0.0] #ym(0) = 1 é para deixar o sistema excitado

#Tempo da simulação (ajustável)
t_span2 = (0, 10) #Intervalo de tempo
t_eval2 = np.linspace(0, 10, 500) #500 pontos de tempo para plotar

#Solução para M
sol_M = solve_ivp(sistema_M, t_span2, y0_M, args=(M, ks, cs), t_eval=t_eval2)

#Solução para m
sol_m = solve_ivp(sistema_m, t_span2, y0_m, args=(m, k1, c1, k2, c2), t_eval=t_eval2)

#Resultados Gráficos dos Deslocamentos
plt.plot(sol_M.t, sol_M.y[0], label='y_M(t)')
plt.plot(sol_m.t, sol_m.y[0], label='y_m(t)')
plt.xlabel('Tempo[s]')
plt.ylabel('Deslocamento[m]')
plt.legend()
plt.title('Sistema Dinâmico Desacoplado')
plt.grid()
plt.show()



#######################################
#SIMULAÇÃO DA SEGUNDA PARTE DO SISTEMA#
#######################################

#Função do movimento da segunda parte do sistema
def sistema2(N, kT, R, x0, v0, tmax, dt):
    #Vetor intervalo de tempo
    t = np.arange(0, tmax, dt)

    # Vetor para armazenar resultados
    x = np.zeros_like(t)
    v = np.zeros_like(t)

    # Condições iniciais
    x[0] = x0
    v[0] = v0

    for i in range(1, len(t)):
        # Atualiza a aceleração
        a = - (kT / (N * R)) * x[i-1]
        
        # Atualiza a velocidade
        v[i] = v[i-1] + a * dt
        
        # Atualiza a posição
        x[i] = x[i-1] + v[i-1] * dt
    return t, x, v

#Condições Iniciais
x0 = 1.0            #m
v0 = 0.0            #m/s   
tmax = 10.0         #s
dt = 0.01           #s

#Resultado da Função
t, x, v = sistema2(N, kT, R, x0, v0, tmax, dt)

#Resultados Gráficos da Segunda Parte do Sistema
plt.figure(figsize=(10, 5))
plt.plot(t, x, label='Posição(x)')
plt.xlabel('Tempo[s]')
plt.ylabel('Deslocamento[m]')
plt.legend()
plt.title('Sistema Dinâmico Pt2')
plt.grid()
plt.show()