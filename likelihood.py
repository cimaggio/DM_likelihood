#!/home/software/anaconda/bin/python

import numpy as np
import scipy as sp
from tqdm import tqdm 
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import math
import config
from scipy.optimize import minimize

from array import array

from PdfOff import off_energy_distr, create_off_interpolations, off_normalization, normoff, numtotoff, minim, get_OffPdf
from PdfOn  import OnInit, integrand_init, PdfOnFunc, integrand, normon, numtoton, ene_on

from models import get_masses, E_lowlims, E_uplims

# input files
EffAreaPath= config.GlobalDir+"extractaeff/effective_areaweighted_0.0_flute_full.root"
OffDir = []
OffDir.append(config.GlobalDir+"st0306/flute_Idx0/Samples_simultBG.root")
OffDir.append(config.GlobalDir+"st0307/flute_Idx0/Samples_simultBG.root")

OnDir  = []
OnDir.append(config.GlobalDir+"st0306/flute_Idx0/Samples.root")
OnDir.append(config.GlobalDir+"st0307/flute_Idx0/Samples_simultBG.root") #this file does not exist!!! 

EffOnTime = []
EffOnTime.append(305812.)
EffOnTime.append(185309.)

Tau      = 1./3

Channel   = 'tautau'    # possible decays: bb, mumu, tautau, WW

# main program 

# initialize model
(mass_interp, masses) = get_masses(config.ModelsDir,Channel,200.,500.)
#get_histogram(ModelsDir,500.,config.Channel)
#print (np.shape(E_lowlims))

#exit(0)

do_plot = False
create_off_interpolations(OffDir, do_plot)
off_normalization(do_plot)

print ('energy: 71.3', get_OffPdf(71.3,0))

OnInit(EffAreaPath,OnDir)
#pdfon = PdfOnFunc()
#print ('pdf =',pdfon)

for m in range (0, len(masses)):

    integrand_init(mass_interp[m],masses[m],E_lowlims[m],E_uplims[m])

    b = normoff  # should be later a nuisance parameter
    sig = 1.     # signal strength, later to be minimized
    pdf  = 0.
    for k in range (0, config.nsamples):   
        
        # create the combined PDF of On and OFF
        pbar = tqdm(total=numtoton[k])
        
        for i in range (0, numtoton[k]):
            pbar.update(1)
            #           res = quad(get_integral_etrue,35., 2000.,args=(k,fi,T,i))
            #           pdfon += math.log(quad(get_integral_etrue,35., 2000.,args=(k,fi,T,i))[0])    #it should be a np array of integrals!
            #            print ('energy: ',ene_on[k][i])
            #            print ('sum: ',np.sum(integrand[k][i]))
            pdf +=  math.log(( Tau * b[k] * get_OffPdf(ene_on[k][i],k)/normoff[k] + sig * EffOnTime[k] * np.sum(integrand[k][i]))
                            / ( Tau * b[k] + normon[k] * EffOnTime[k] ) )
        pbar.close()

        # Add the Poissonian fluctuations of the total number of events
        pdf +=  numtoton[k] * math.log(normon[k]*EffOnTime[k])  - normon[k]*EffOnTime[k]
        pdf +=  numtotoff[k]* math.log(b[k])       - b[k]

    pdf = -2. * pdf

    print ('pdf =',pdf)

#print (log_like(150000.,0))
#(b_est,l,h) = minim(0)
#print ("b_est:", b_est," low lim: ",l," high lim: ", h)
