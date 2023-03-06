#!/usr/bin/env python3

import numpy as np
import random
import matplotlib.pyplot as plt


def random_sign():
    if random.random() < 0.5:
        return -1
    else:
        return 1


class Photon:
    def __init__(self, position, direction, energy):
        """
        A photon with a position (x, y) and an energy.
        :param position:    2D coordinates in cm
        :type position:     numpy.ndarray
        :param direction:   direction of the photon in [0, 2pi)
        :type direction:    float
        :param energy:      photon energy in MeV
        :type energy:       float
        """
        self.position = position
        self.x = self.position[0]
        self.y = self.position[1]
        self.direction = direction
        self.energy = energy

    def set_position(self, position):
        self.position = position
        self.x = position[0]
        self.y = position[1]

    def set_direction(self, direction):
        self.direction = direction

    def set_energy(self, energy):
        self.energy = energy

    def in_rectangle(self, rectangle):
        return (rectangle[0, 0] <= self.x <= rectangle[1, 0]) and (rectangle[0, 1] <= self.y <= rectangle[1, 1])

    def do_rayleigh(self):
        zeta = random.random()
        xi = random.random()

        # p_1 = I_1(π) / ( I_1(π) + I_2(π) )
        # with I_1(π) = integral sin(θ) dθ from 0 to π
        #      I_2(π) = integral cos²(θ) sin(θ) dθ from 0 to π
        p_1 = 3/4

        if zeta <= p_1:
            polar_angle = np.arccos(1 - 2*zeta)
        elif xi <= 0.5:
            polar_angle = np.arccos((1 - 2*xi) ** (1/3))
        else:
            polar_angle = np.arccos(- (2*xi - 1) ** (1/3))

        self.set_direction(self.direction + random_sign() * polar_angle)

    def do_photoelectric(self):
        self.set_energy(0)

    def do_compton(self):
        unsuccessful = True

        while unsuccessful:
            eta_1 = random.random()
            eta_2 = random.random()
            zeta = random.random()

            E_prime_min = 1 / (1/self.energy + 2/0.511)
            E_prime_max = self.energy

            A_1 = 2 / (E_prime_max**2 - E_prime_min)
            A_2 = 1 / (np.log(E_prime_max/0.511) - np.log(E_prime_min/0.511))
            # I'm unsure about the
            # energy units here...

            k_1 = 1/(A_1 * self.energy)
            k_2 = self.energy/A_2

            p_1 = k_1 / (k_1 + k_2)

            if eta_1 < p_1:
                xi = np.sqrt(np.abs(E_prime_min**2 + 2*eta_2/A_1))
            else:
                xi = E_prime_min * np.exp(eta_2/A_2)

            polar_angle = np.arcsin(np.sqrt(1 - (1 - 0.511/xi + 0.511/self.energy)**2))
            # I'm unsure about the
            # energy units here...

            fcg_quotient = 1 - np.sin(polar_angle) ** 2 / (xi/self.energy + self.energy/xi)

            if zeta <= fcg_quotient:
                self.set_energy(xi)
                self.set_direction(self.direction + random_sign() * polar_angle)
                unsuccessful = False

    def do_pair_production(self):
        self.set_energy(0)


def find_nearest_attenuation_coefficient(list, energy, density):
    """
    Returns the best fitting attenuation coefficient from a list
    :param list:    array of energies and attenuation coefficients
    :type list:     numpy.ndarray
    :param energy:  energy to be used (in MeV!)
    :type energy:   float
    :return:        float
    """
    idx = (np.abs(list[:, 0] - energy)).argmin()
    return list[idx, 1] * density


def find_nearest_energy_absorption_coefficient(list, energy, density):
    """
    Returns the best fitting energy absorption coefficient from a list
    :param list:    array of energies and energy absorption coefficients
    :type list:     numpy.ndarray
    :param energy:  energy to be used (in MeV!)
    :type energy:   float
    :return:        float
    """
    idx = (np.abs(list[:, 0] - energy)).argmin()
    return list[idx, 2] * density


def find_nearest_cross_sections(list, energy):
    """
    Returns the best fitting cross-sections from a list
    :param list:    array of energies and attenuation coefficients
    :type list:     numpy.ndarray
    :param energy:  energy to be used (in MeV!)
    :type energy:   float
    :return:        numpy.ndarray
    """
    idx = (np.abs(list[:, 0] - energy)).argmin()
    return list[idx, 1:]


def cross_sections_to_probabilities(cross_sections):
    """
    Converts a list of cross-sections to probabilities where the last list entry is the total cross-section.
    They should be in this order: Rayleigh, Compton, photoelectric effect, pair production (n), pair production (e).
    :param cross_sections:  list of cross-sections
    :type cross_sections:   numpy.ndarray
    :return:                numpy.ndarray
    """
    probabilities = cross_sections[:5]/cross_sections[5]
    return probabilities


# water region (depth: x, width: y)
DEPTH = 10  # cm
WIDTH = 5   # cm

# coordinates of beam source
BEAM_SOURCE = np.array([0, 0])
water_region = np.array([[BEAM_SOURCE[0], BEAM_SOURCE[1]-WIDTH/2],
                         [BEAM_SOURCE[0]+DEPTH, BEAM_SOURCE[1]+WIDTH/2]])

# initial energy
INIT_ENERGY = 0.1  # MeV

# initial direction (pencil beam: 0)
INIT_DIRECTION = 0

# number and mass density of water in @ 20 °C and 1 bar
NUMBER_DENSITY = 1E8  # 1/cm²
MASS_DENSITY = 998.23E-3  # g/cm³

# attenuation coefficients of water
# Energy (MeV) | mu/rho (cm²/g) | mu_en/rho (cm²/g)
ATTEN_COEFF = np.load('attenuation_coefficients.npy')

# cross-sections of photons in water
# Energy | Rayleigh | Compton | photoelectric effect | pair production (nuclear) | pair production (electron) | total
CROSS_SECTIONS = np.load('cross_sections.npy')

# number of events
N_EVENTS = int(1E4)

if __name__ == '__main__':

    deposited_energy = 0
    for i in range(N_EVENTS):
        # create photon
        print('Creating new photon...')
        photon = Photon(BEAM_SOURCE, INIT_DIRECTION, INIT_ENERGY)

        positions = np.array(BEAM_SOURCE)

        consider_photon = True
        while consider_photon:

            # find attenuation coefficient for the closest energy in table
            attenuation_coefficient = find_nearest_attenuation_coefficient(ATTEN_COEFF, photon.energy, MASS_DENSITY)

            # sample distance to next interaction
            distance = - 1/attenuation_coefficient * np.log(random.random())

            # assign new position to photon
            photon.set_position(photon.position + distance * np.array([np.cos(photon.direction), np.sin(photon.direction)]))
            positions = np.append(positions, photon.position)

            # check whether photon is still inside water region
            if photon.in_rectangle(water_region):
                cross_sections = find_nearest_cross_sections(CROSS_SECTIONS, photon.energy)
                probabilities = cross_sections_to_probabilities(cross_sections)
                interaction_sampling = random.random()

                if interaction_sampling <= probabilities[0]:
                    # Rayleigh scattering
                    photon.do_rayleigh()

                elif interaction_sampling <= probabilities[0] + probabilities[2]:
                    # photoelectric absorption
                    photon.do_photoelectric()
                    print('Photoelectric effect.')

                elif interaction_sampling <= np.sum(probabilities[:3]):
                    # Compton scattering
                    photon.do_compton()

                else:
                    # pair production
                    photon.do_pair_production()
                    print('Pair production.')

                if photon.energy == 0:
                    consider_photon = False

            else:
                print('The photon is gone.')
                consider_photon = False

        deposited_energy += (INIT_ENERGY - photon.energy)

        # positions = positions.reshape((int(len(positions)/2), 2))
        # print(positions)
#
        # plt.plot(positions[:, 0], positions[:, 1])
        # plt.xlim(water_region[0, 0], water_region[1, 0])
        # plt.ylim(water_region[0, 1], water_region[1, 1])
        # plt.show()

    print('Simulated energy deposition: ' + str(deposited_energy/N_EVENTS*NUMBER_DENSITY) + ' MeV')

    expected_energy = INIT_ENERGY * NUMBER_DENSITY * find_nearest_energy_absorption_coefficient(ATTEN_COEFF, INIT_ENERGY, MASS_DENSITY) * DEPTH

    print('Expected energy:             ' + str(expected_energy) + ' MeV')
