# Prova Simulação Numérica 02
# Autores:
#   Michel Bernardino Machado
#   Núbia Andrade Bergamini
# Data: 14/02/2021

import matplotlib.pyplot as plt
import numpy as np
import cantera as ct

# Temperatura e pressão
T_combustor = np.arange(300,800,100)
P_combustor = 607950

# Razão de equivalência
phi_1 = 0.5
phi_2 = 0.8
phi = 2*phi_2 - phi_1

# Volume do combustor [m^3]
vol_1 = 7.2*10**(-5) 
vol_2 = 8.6*10**(-5)

# Tempo de residência
initial_residence_time = 1
extinction_residence_time = 1.7*10**(-5)

# Componentes da mistura
fuel_species = 'CH4:1.0'
oxidizer_species = 'O2:1.0, N2:3.76'
mixture_species = fuel_species + ', ' + oxidizer_species

for i in range(6):

    # Mecanismo utilizado
    mixture_1 = ct.Solution('gri30etanol.cti')
    mixture_1.TP = T_combustor[i], P_combustor
    mixture_1.set_equivalence_ratio(phi_1, fuel_species, oxidizer_species)

    mixture_2 = ct.Solution('gri30etanol.cti')
    mixture_2.TP = T_combustor[i], P_combustor
    mixture_2.set_equivalence_ratio(phi, fuel_species, oxidizer_species)

    # Reservatório para armazenar os reagentes 
    mixture_1_inlet = ct.Reservoir(mixture_1)
    mixture_1.equilibrate('HP')
    mixture_2_inlet = ct.Reservoir(mixture_2)
    mixture_2.equilibrate('HP')

    # Criação do combustor com duas zonas
    reactor_1 = ct.Reactor(contents = mixture_1, volume = vol_1, name = 'Primary zone')
    reactor_2 = ct.Reactor(contents = mixture_2, volume = vol_2, name = 'Secondary zone')

    # Reservatório para armazenar os gases de exaustão
    exhaust = ct.Reservoir(mixture_1)

    def mdot_1(t):
        return reactor_1.mass + reactor_2.mass / residence_time

    def mdot_2(t):
        return reactor_2.mass / residence_time 

    # Controladores de fluxo de massa das misturas de ar e combustível
    inlet_1_mfc = ct.MassFlowController(upstream = mixture_1_inlet, downstream = reactor_1, mdot = mdot_1) 
    inlet_2_mfc = ct.MassFlowController(upstream = mixture_2_inlet, downstream = reactor_2, mdot = mdot_2)

    # Fluxo de mistura do reator 1 para o reator 2
    inter_reactors_mfc = ct.Valve(upstream = reactor_1, downstream = reactor_2, K = 1.0) 

    # Fluxo de mistura do reator 2 para o reservatório
    outlet_mfc = ct.Valve(upstream = reactor_2, downstream = exhaust, K = 1.0)

    sim = ct.ReactorNet([reactor_1, reactor_2])

    states = ct.SolutionArray(mixture_1, extra = ['tres'])
    states_2 = ct.SolutionArray(mixture_2, extra = ['tres'])

    residence_time = initial_residence_time
    while residence_time > extinction_residence_time:
        sim.set_initial_time(0.0)  
        sim.advance_to_steady_state(max_steps = 100000)
        print('tres = {:.2e}; T1 = {:.1f}; T2 = {:.1f}; Phi1 = {:.1f}; Phi2 = {:.2f}'.format(residence_time, reactor_1.T, reactor_2.T, mixture_1.equivalence_ratio(basis = 'mole'), mixture_2.equivalence_ratio(basis = 'mole')))
        states.append(reactor_1.thermo.state, tres = residence_time)
        states_2.append(reactor_2.thermo.state, tres = residence_time)
        residence_time *= 0.9  # decrease the residence time for the next iteration

    # Gráficos dos resultados
    plt.plot(states.tres[:-1], states.T[:-1], '.-', color='C1')
    plt.xlabel('Tempo de residência [s]')
    plt.ylabel('Temperatura de chama adiabática [K]')
    plt.show()

    plt.plot(states_2.tres[:-1], states_2.T[:-1], '.-', color='C1')
    plt.xlabel('Tempo de residência [s]')
    plt.ylabel('Temperatura de chama adiabática [K]')
    plt.show()

