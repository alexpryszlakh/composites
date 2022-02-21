# Effective Lamina Properties
# Author: Alex Hernandez-Pryszlak, Daniel Donate-Perez
# Date:   2/15/2022

# import pandas
import math
import numpy as np

# Elastic Constant
# Inputs

units = int(input('PSI (1) or Pa (2) (enter 1 or 2): '))
Vf = float(input('Fiber Volume Fraction, Vf: '))  # Volume Fraction

E1f = float(input('Elastic Modulus in Fiber Direction, E1f [Pa]: '))
E2f = float(input('Elastic Modulus in Transverse Direction, E2f [Pa]: '))
G12f = float(input('Shear Modulus in 1-2 Plane, G12f [Pa]: '))
v12f = float(input('Poisson Ratio in 1-2 Plane, v12f: '))
Em = float(input('Elastic Modulus of Matrix, Em [Pa]: '))
vm = float(input('Poisson Ratio of Matrix, vm: '))

Gm = Em / (2 * (1 + vm))
E1 = E1f * Vf + Em * (1 - Vf)  # Rule of Mixtures
v12 = v12f * Vf + vm * (1 - Vf)  # Rule of Mixtures
E2 = 1 / ((1 - math.sqrt(Vf)) / Em + math.sqrt(Vf) / (
        E2f * math.sqrt(Vf) + Em * (1 - math.sqrt(Vf))))  # Improved Mechanics of Materials model
print(E2)
G12 = Gm * (((Gm + G12f) - Vf * (Gm - G12f)) / ((Gm + G12f) + Vf * (Gm - G12f)))  # CCM
print(G12)

# Outputs
# lamina_properties = [(E1f, E2f, G12f, v12f), (Em, Em, Gm, vm), (E1, E2, G12, v12)]
lamina_properties = [E1f, E2f, G12f, v12f, Em, Em, Gm, vm, E1, E2, G12, v12] # creating the list of properties 

# final_prop = np.array(lamina_properties) * 1e-9

# if else statement to then obtain proper results
if units == 1:
    final_prop = np.array(lamina_properties) * 1e-6
    print(final_prop)
else:
    final_prop = np.array(lamina_properties) * 1e-9
    print(final_prop)

'''def lamina():
    units = int(input('PSI (1) or Pa (2) (enter 1 or 2): '))
    Vf = float(input('Fiber Volume Fraction, Vf: '))  # Volume Fraction

    E1f = float(input('Elastic Modulus in Fiber Direction, E1f [Pa]: '))
    E2f = float(input('Elastic Modulus in Transverse Direction, E2f [Pa]: '))
    G12f = float(input('Shear Modulus in 1-2 Plane, G12f [Pa]: '))
    v12f = float(input('Poisson Ratio in 1-2 Plane, v12f: '))
    Em = float(input('Elastic Modulus of Matrix, Em [Pa]: '))
    vm = float(input('Poisson Ratio of Matrix, vm: '))

    Gm = Em / (2 * (1 + vm))
    E1 = E1f * Vf + Em * (1 - Vf)  # Rule of Mixtures
    v12 = v12f * Vf + vm * (1 - Vf)  # Rule of Mixtures
    E2 = 1 / ((1 - math.sqrt(Vf)) / Em + math.sqrt(Vf) / (
            E2f * math.sqrt(Vf) + Em * (1 - math.sqrt(Vf))))  # Improved Mechanics of Materials model
    G12 = Gm * (((Gm + G12f) - Vf * (Gm - G12f)) / ((Gm + G12f) + Vf * (Gm - G12f)))  # CCM
    return Gm, E1, v12, E2, G12'''
