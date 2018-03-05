#!/home/software/anaconda/bin/python

import sys
sys.path.append('/media/san/astro/soft/root/root_v5.34.36/root/lib')
import ROOT

from ROOT import TFile, TObject, TH1, TAxis, TArrayD
from pathlib import Path
import config
import numpy as np
from scipy.interpolate import interp1d

inter_mass = []

masses = [ 100., 110., 120., 130., 140., 150., 160., 180.,
           200., 220., 240., 260., 280.,
           300., 330., 360.,
           400., 450.,
           500., 550.,
           600., 650.,
           700., 750.,
           800., 900.,
           1000., 1100., 1200.,1300., 1500.,1700., 
           2000., 2500., 
           3000., 4000., 5000., 6000., 7000., 8000., 9000.,
           10000.,12000., 15000., 20000., 30000., 50000., 100000. ]
           

def get_mass_histogram(directory,mass,channel):

    file = directory + "/dNdESignal_" + channel + "_" + str(int(mass)) + ".0mass.root"
    print ("Looking for file: ",file)

    if (not Path(file).is_file()):
        err_message = 'File \"'+ file + '\" does not exist!'
        print (err_message)
        raise ValueError(err_message)

    DoHisto = TFile.Open(file, "READ")
    print(DoHisto)  
    
    histogram = DoHisto.Get("hdNdE")
    print (histogram)

    bins = histogram.GetNbinsX()

    model_E  = np.empty(bins)
    model_N  = np.empty(bins)
    
    for ebin in range(0, bins):
        model_E[ebin] = histogram.GetXaxis().GetBinCenter(ebin+1)
        model_N[ebin] = histogram.GetBinContent(ebin+1)

    return interp1d(model_E,model_N, kind=config.ModelInterpolation)


def get_mass_interpolations(directory, masses, channel):
    
    global inter_mass

    for mass in masses:
        inter_mass.append(get_mass_histogram(directory,mass,channel))

def get_masses(mass_dir,channel,mmin,mmax):

    global inter_mass

    nmasses = len(masses)
    print (nmasses)
    real_masses = np.empty(nmasses)
    print ('nmasses: ',nmasses, ' type: ', type(nmasses), ' shape: ',np.shape(real_masses))    

    reali = 0
    for mass in masses:
        if (mass < mmin):
            continue
        if (mass > mmax):
            continue
        print ('i: ',reali,' mass: ',mass, 'shape: ',np.shape(real_masses))
        real_masses[reali] = mass
        reali += 1 

    real_masses = np.resize(real_masses,reali)

    get_mass_interpolations(mass_dir, real_masses, channel)

    return inter_mass
