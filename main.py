#===================================================================#
#     ANÁLISE DE VIBRAÇÕES LIVRES EM UMA VIGA ENGASTADA E LIVRE     #
#===================================================================#

import numpy as np
from math import *
from quadratura_gaussiana import quadraturaGaussiana

# Entrada
b = float(input('bw: '))
h = float(input('h: '))
L = float(input('L: '))
E = float(input('E: '))
m = float(input('Massa específica: '))

# Cálculo de propriedades
I = b*h**3/12
m_linear = m*b*h

# Domínio discreto
x = np.arange(0, L+L/100, L/100)

# Funções de forma e suas segundas derivadas no domínio discreto
fi_cosseno = 1 - np.cos(x*pi/(2*L))
fi_cosseno_xx = (pi/(2*L))**2*np.cos(x*pi/(2*L))

fi_polinomio = -x**3/(2*L**3)+3*x**2/(2*L**2)
fi_polinomio_xx = -6*x/(2*L**3)+6/(2*L**2)

# Integração
m_modal_cosseno = quadraturaGaussiana(10, 0, L, x, m_linear*fi_cosseno**2)
k_modal_cosseno = quadraturaGaussiana(10, 0, L, x, E*I*fi_cosseno_xx**2)

m_modal_polinomio = quadraturaGaussiana(5, 0, L, x, m_linear*fi_polinomio**2)
k_modal_polinomio = quadraturaGaussiana(5, 0, L, x, E*I*fi_polinomio_xx**2)

# Saída
print(50*'=')
print(f'Função de forma triginométrica:')
print(f'Massa modal: {m_modal_cosseno} kg')
print(f'Rigidez modal: {k_modal_cosseno} N.m')
print(f'Frequência natural: {sqrt(k_modal_cosseno/m_modal_cosseno)/(2*pi)} Hz')
print(40*'-')
print(f'Função de forma polinomial:')
print(f'Massa modal: {m_modal_polinomio} kg')
print(f'Rigidez modal: {k_modal_polinomio} N.m')
print(f'Frequência natural: {sqrt(k_modal_polinomio/m_modal_polinomio)/(2*pi)} Hz')
print(50*'=')
