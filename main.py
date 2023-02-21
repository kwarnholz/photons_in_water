#!/usr/bin/env python3

import numpy as np
import random
import matplotlib.pyplot as plt

class Photon:
    def __init__(self, position, direction, energy):
        """
        A photon with a position (x, y) and an energy.
        :param position:    2D coordinates in cm
        :type position:     numpy.ndarray
        :param direction:   direction of the photon in [0, 2pi)
        :type direction:    float
        :param energy:      photon energy in eV
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

    def in_rectangle(self, rectangle):
        return (rectangle[0, 0] <= self.x <= rectangle[1, 0]) and (rectangle[0, 1] <= self.y <= rectangle[1, 1])


def find_nearest_attenuation_coefficient(list, energy):
    """
    Returns the best fitting attenuation coefficient from a list
    :param list:    array of energies and attenuation coefficients
    :type list:     numpy.ndarray
    :param energy:  energy to be used (in MeV!)
    :type energy:   float
    :return:        float
    """
    idx = (np.abs(list[:, 0] - energy)).argmin()
    return list[idx, 1]


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


# water region (depth: x, width: y)
DEPTH = 10  # cm
WIDTH = 5   # cm

# coordinates of beam source
BEAM_SOURCE = np.array([0, 0])
water_region = np.array([[BEAM_SOURCE[0], BEAM_SOURCE[1]-WIDTH/2],
                         [BEAM_SOURCE[0]+DEPTH, BEAM_SOURCE[1]+WIDTH/2]])

# initial energy
INIT_ENERGY = 100E3  # eV

# initial direction (pencil beam: 0)
INIT_DIRECTION = 0

# number and mass density of water in 1/cm³ @ 20 °C and 1 bar
NUMBER_DENSITY = 33.3679E21
MASS_DENSITY = 998.23E-3  # g/cm³

# attenuation coefficients of water
# Energy (MeV) | mu/rho (cm²/g) | mu_en/rho (cm²/g)
ATTEN_COEFF = np.load('attenuation_coefficients.npy')

# cross-sections of photons in water
# Energy | Rayleigh | Compton | photoelectric effect | pair production (nuclear) | pair production (electron) | total
CROSS_SECTIONS = np.load('cross_sections.npy')

# number of events
# N_EVENTS = 1E4
N_EVENTS = 1

if __name__ == '__main__':
    for i in range(N_EVENTS):
        # create photon
        photon = Photon(BEAM_SOURCE, INIT_DIRECTION, INIT_ENERGY)

        # find attenuation coefficient for the closest energy in table
        attenuation_coefficient = find_nearest_attenuation_coefficient(ATTEN_COEFF, photon.energy/1E6) * MASS_DENSITY
        print(attenuation_coefficient)

        # sample distance to next interaction
        distance = - 1/attenuation_coefficient * np.log(random.random())
        print(distance)

        # assign new position to photon
        photon.set_position(photon.position + distance * np.array([np.cos(photon.direction), np.sin(photon.direction)]))
        print(photon.position)

        # check whether photon is still inside water region
        if photon.in_rectangle(water_region):
            print('Do stuff!')
        else:
            break
