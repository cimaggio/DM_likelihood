#!/home/software/anaconda/bin/python

import config

def pass_energy_cuts(energy):

    if (energy < config.energy_cut):
        return False
    if (energy > config.energy_upper_cut):
        return False

    return True
        
