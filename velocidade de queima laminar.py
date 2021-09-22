# Prova: Simulação Numérica 03
# Autor: Michel Bernardino Machado
# Data: 16/03/2021

# Pacotes utilizados
import matplotlib.pyplot as plt
import numpy as np
import cantera as ct

# Parâmetros da simulação 
Pin = 709275 
Tin = np.linspace(400, 650, 6)
Phi = np.arange(0.1, 3, 0.05) # para combustão do H2, alterar Phi[0] para um valor maior que 0.1 
width = 0.03
loglevel = 1

# Componentes da mistura
# fuel_species = 'H2' 
fuel_species = 'C12H23'
oxidizer_species = 'O2:1.0, N2:3.76'

# Mecânismos utilizados 
# mixture = ct.Solution('gri30.yaml')
mixture = ct.Solution('C12H23.cti')

# Matriz para armazenar a velocidade de chama
flame_velocity_1 = np.zeros((len(Tin),len(Phi)))

# Função para plotagem dos gráficos
def flame_velocity_plot(flame_velocity_1):
    
    for k in range(len(Tin)):

        plt.plot(Phi, flame_velocity_1[k, :], label = '$T = $'+ str(Tin[k]) + ' $K$')
        plt.xlabel('Razão de equivalência')
        plt.ylabel('Velocidade de queima [m/s]')
        plt.legend()
        plt.title('VELOCIDADE DE QUEIMA X RAZÃO DE EQUIVALÊNCIA')

    plt.show()

    for k in range(len(Tin)):

        plt.plot(Phi, 39.3701*flame_velocity_1[k, :], label = '$T = $'+ str(Tin[k]) + ' $K$')
        plt.xlabel('Razão de equivalência')
        plt.ylabel('Velocidade de queima [in/s]')
        plt.legend()
        plt.title('VELOCIDADE DE QUEIMA X RAZÃO DE EQUIVALÊNCIA')

    plt.show()

for i in range(len(Tin)):

    for j in range(len(Phi)):

        # Define as condições iniciais da mistura a ser queimada
        mixture.TP = Tin[i], Pin
        mixture.set_equivalence_ratio(Phi[j], fuel_species, oxidizer_species)

        # Cria um objeto para representar a chama
        f = ct.FreeFlame(mixture, width = width)
        f.set_refine_criteria(ratio = 3, slope = 0.06, curve = 0.12)
        f.show_solution()

        # Resolve utilizando o modelo de tranporte de mistura média
        f.transport_model = 'Mix'
        f.solve(loglevel = loglevel, auto = True)
        f.show_solution()
        flame_velocity_1[i, j] = f.velocity[0]

flame_velocity_plot(flame_velocity_1)


