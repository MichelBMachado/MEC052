# Prova: Simulação Numérica 01
# Autor: Michel Bernardino Machado
# Data: 26/01/2021

# Instruções de uso:
# Para usar o código é necessário mover o arquivo: "gri30etanol.cti" para a pasta data.

# Bibliotecas utilizadas
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Temperatura dos reagente variando de 600 k  à 1300 K.
T = np.linspace(600, 1300, 5)

# Pressão dos reagentes, é mantida constante ao longo do processo.
P = 1.8*10**6

# Fases
gas = ct.Solution('gri30etanol.cti')
carbon = ct.Solution('graphite.xml')
mix_phases = [(gas, 1.0), (carbon, 0.0)]

# Espécies químicas.
fuel_species = 'CH4'
#fuel_species = 'C2H5OH'
oxidizer_species = 'O2:1.0, N2:3.76'

# Razão de equivalência.
npoints = 50
phi = np.linspace(0.1, 4.0, npoints)

# Mistura dos reagentes.
mix = ct.Mixture(mix_phases)

# Matrizes para armazenar os dados.
tad = np.zeros(npoints)
xeq = np.zeros((mix.n_species,npoints))

CO = np.zeros((npoints, 5))
CO2 = np.zeros((npoints, 5))
NO = np.zeros((npoints, 5))
NO2 = np.zeros((npoints, 5))

for i in range(5):

    for j in range(npoints):

        # Define o estado do gás.
        gas.set_equivalence_ratio(phi[j], fuel_species, oxidizer_species)

        # Cria uma mistura com 1 mol de gás e 0 mol de carbono.
        mix = ct.Mixture(mix_phases)
        mix.T = T[i]
        mix.P = P

        # Equilibra a mistura adiabaticamente à pressão constante.
        mix.equilibrate('HP', solver='gibbs', max_steps=1000)

        tad[j] = mix.T
        xeq[:,j] = mix.species_moles
        
        NO[j,i] = (xeq[mix.species_index(gas, 'NO'),j] / mix.phase_moles(gas))
        NO2[j,i] = (xeq[mix.species_index(gas, 'NO2'),j] / mix.phase_moles(gas))
        CO[j,i] = (xeq[mix.species_index(gas, 'CO'),j] / mix.phase_moles(gas))
        CO2[j,i] = (xeq[mix.species_index(gas, 'CO2'),j] / mix.phase_moles(gas))

    plt.plot(phi, tad, label = '$T = $'+ str(T[i]) + ' $K$')
    plt.xlabel('Razão de equivalência')
    plt.ylabel('Temperatura de chama adiabática [K]')
    plt.legend()
    plt.text(1.0, 1000,fuel_species, fontsize = 12)

plt.show()

# Função para plotagem dos gráficos de fração molar por razão de equivalência
def molar_fraction_plot(species_name, species_molar_fraction):
    for k in range(5):
        plt.plot(phi, species_molar_fraction[:,k], label = '$T = $'+  str(T[k]) + ' $K$')
        plt.xlabel('Razão de equivalência')
        plt.ylabel('Fração molar de ' + species_name)
        plt.legend()
        plt.text(3.5, 0.01,fuel_species, fontsize = 12)
    plt.show()

species_name = 'CO'
species_molar_fraction = CO
molar_fraction_plot(species_name, species_molar_fraction)

species_name = 'CO2'
species_molar_fraction = CO2
molar_fraction_plot(species_name, species_molar_fraction)

species_name = 'NO'
species_molar_fraction = NO
molar_fraction_plot(species_name, species_molar_fraction)

species_name = 'NO2'
species_molar_fraction = NO2
molar_fraction_plot(species_name, species_molar_fraction)





