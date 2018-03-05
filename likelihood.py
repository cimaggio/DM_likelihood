#!/home/software/anaconda/bin/python

import sys
sys.path.append('/media/san/astro/soft/root/root_v5.34.36/root/lib')
import ROOT
import numpy as np
import scipy as sp
from tqdm import tqdm 
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import math
import config
from scipy.optimize import minimize

from ROOT import TFile, TObject, TH1, TMath, TNtuple, TCanvas 
from array import array

from PdfOff import off_energy_distr, create_off_interpolations, off_normalization, normoff, log_like, minim, inter_func
from PdfOn  import OnInit, PdfOnFunc, integrand, normon, numtoton, ene_on

# input files
EffAreaPath="/media/san2/astro/M15/extractaeff/effective_areaweighted_0.0_flute_full.root"
OffDir = []
OffDir.append("/media/san2/astro/M15/st0306/flute_Idx0/Samples_simultBG.root")
OffDir.append("/media/san2/astro/M15/st0307/flute_Idx0/Samples_simultBG.root")

OnDir  = []
OnDir.append("/media/san2/astro/M15/st0306/flute_Idx0/Samples.root")
OnDir.append("/media/san2/astro/M15/st0307/flute_Idx0/Samples_simultBG.root") #this file does not exist!!! 

EffOnTime = []
EffOnTime.append(305812.)
EffOnTime.append(185309.)

Tau      = 1./3

# main program 
OnInit(EffAreaPath,OnDir,EffOnTime)
#pdfon = PdfOnFunc()
print ('pdf =',pdfon)

b = normoff  # should be later a nuisance parameter

for k in range (0, config.numsamples):   

    # create the combined PDF of On and OFF
    
    pbar = tqdm(total=numtoton[k])
    
    for i in range (0, numtoton[k]):
           pbar.update(1)
           #           res = quad(get_integral_etrue,35., 2000.,args=(k,fi,T,i))
           #           pdfon += math.log(quad(get_integral_etrue,35., 2000.,args=(k,fi,T,i))[0])    #it should be a np array of integrals!
           pdfon +=  math.log(( Tau * b[k] * inter_func[k](ene_on[k][i])/normoff[k] 
                                + np.sum(integrand[k][i]))   /  
                              ( Tau * b[k] + normon[k] ) )
           
    pbar.close()




do_plot = False
#create_off_interpolations(OffDir, do_plot)
#off_normalization(do_plot)
#print (log_like(150000.,0))
#(b_est,l,h) = minim(0)
#print ("b_est:", b_est," low lim: ",l," high lim: ", h)
