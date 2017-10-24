#!/home/software/anaconda/bin/python

import sys
sys.path.append('/media/san/astro/soft/root/root_v5.34.36/root/lib')
import ROOT
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import math
from scipy.optimize import minimize

from ROOT import TString, TFile, TObject, TH1, TMath, TNtuple, TCanvas 
from array import array
#from ROOT import *   (this command line does not work in a script)

# all global parameters here
numtot = 0
content = 0
ene = np.empty(shape=(540000),dtype='double')
enest = []
logenergy = []
pEnBiasVsEnEst = []
pEnBiasVsEnEst2 = []
hEnResVsEnEst = []
hEnResVsEnEst2 = []
hEnPartVsEnEst = []
hEnMinVsEnEst = []
hEnMaxVsEnEst = []
nbinsenest = []
est_en_edges = np.empty(shape=(2, 17), dtype='double')

sEffArea  = []
sEffAreaErec = []
AeffEnMax = []

# all cuts here
energy_cut = 50    # GeV
energy_upper_cut = 8000    # GeV

PrelimPath="/media/san2/astro/M15/extractaeff/effective_areaweighted_0.0_flute_full.root"
Dir06="/media/san2/astro/M15/st0306/flute_Idx0/Samples.root"
Dir07="/media/san2/astro/M15/st0307/flute_Idx0/Output_glike.root" #this file does not exist!!! 

#-----------------------------------------------------------------------------------------------------------------------------------------

def Init():

#    NumBckgEnBins = 40
    print("\n\n--- INITIALIZING VALUES from:", PrelimPath, "\n")

    # TFile *LivParamFile  = TFile.Open(PrelimPath, "READ")
    LivParamFile2 = TFile.Open(PrelimPath, "READ")

    SampleName = ["0306", "0307"]

    global pEnBiasVsEnEst
    global pEnBiasVsEnEst2
    global hEnResVsEnEst
    global hEnResVsEnEst2
    global hEnPartVsEnEst
    global hEnMinVsEnEst
    global hEnMaxVsEnEst
    global nbinsenest
    
    global sEffArea 
    global sEffAreaErec
    global AeffEnMax
    

    for k in range (0, 2):    # k= number of samples

        pEnBiasVsEnEst.append(LivParamFile2.FindObjectAny(ROOT.Form("pEnBiasVsEnEst__%s_smoothed" % SampleName[k])))
        pEnBiasVsEnEst[k].SetDirectory(0)        

        pEnBiasVsEnEst2.append(LivParamFile2.FindObjectAny(ROOT.Form("pEnBiasVsEnEst2__%s_smoothed" % SampleName[k])))
        pEnBiasVsEnEst2[k].SetDirectory(0)        

        hEnResVsEnEst.append(LivParamFile2.FindObjectAny(ROOT.Form("hEnResVsEnEst__%s_smoothed" % SampleName[k])))
        hEnResVsEnEst[k].SetDirectory(0)        

        hEnResVsEnEst2.append(LivParamFile2.FindObjectAny(ROOT.Form("hEnResVsEnEst2__%s_smoothed" % SampleName[k])))
        hEnResVsEnEst2[k].SetDirectory(0)        

        hEnPartVsEnEst.append(LivParamFile2.FindObjectAny(ROOT.Form("hEnPartVsEnEst__%s_smoothed" % SampleName[k])))
        hEnPartVsEnEst[k].SetDirectory(0)       
 
        nbinsenest.append(pEnBiasVsEnEst[k].GetNbinsX())

        print("Est Energy Edges for period st%s" % SampleName[k], "\n")

        print ("nbinsenest =", nbinsenest[k], "\n")

        for i in range (0, nbinsenest[k]):
            est_en_edges[k][i] = pEnBiasVsEnEst[k].GetBinLowEdge(i+1) 
            print( pEnBiasVsEnEst[k].GetBinLowEdge(i+1), flush = True)
        
        print ("\n")

        sEffArea.append(LivParamFile2.FindObjectAny(ROOT.Form("sEffArea__%s" % SampleName[k])))
        sEffAreaErec.append(LivParamFile2.FindObjectAny(ROOT.Form("sEffAreaErec__%s" % SampleName[k])))
        

        AeffEnMax.append(sEffArea[k].GetXmax())

    LivParamFile2.Close()
    LivParamFile2.Delete()

    print("--- INITIALIZATION FINISHED!\n")
    return;

#-----------------------------------------------------------------------------------------------------------------------------------------
def gaussian(x, mu, sigma):
          return 1./sigma/math.sqrt(2*math.pi)*np.exp(-0.5*np.power((x - mu)/sigma, 2.))

def get_area(entrue,k=0):
    area = sEffArea[k].Eval(math.log10(entrue))
    print ('log energy: ',math.log10(entrue), 'area: ', area)
    return math.pow(10.,area)

def get_enmigmat(entrue,k=0):
    DoHisto06 = TFile.Open(Dir06, "READ")
    non = DoHisto06.Get("fOnSample")
    numtoton = non.GetEntries()
    global enest 

    for i in range (0, numtoton):
       non.GetEntry(i)
       enest = non.GetArgs()[0]
    
    penbias = np.empty(shape=(nbinsenest[k]))
    penbias2 = np.empty(shape=(nbinsenest[k]))
    henres = np.empty(shape=(nbinsenest[k]))
    henres2 = np.empty(shape=(nbinsenest[k]))
    henpart = np.empty(shape=(nbinsenest[k]))
 
    for i in range (0, nbinsenest[k]):
        penbias[i] = pEnBiasVsEnEst[k].GetBinContent(i+1)
        penbias2[i] = pEnBiasVsEnEst2[k].GetBinContent(i+1)
        henres[i] = hEnResVsEnEst[k].GetBinContent(i+1)
        henres2[i] = hEnResVsEnEst2[k].GetBinContent(i+1)
        henpart[i] = hEnPartVsEnEst[k].GetBinContent(i+1) 

    inter_mean1 = interp1d(est_en_edges[k], penbias, kind='linear')
    inter_mean2 = interp1d(est_en_edges[k], penbias2, kind='linear')    
    inter_sigma1 = interp1d(est_en_edges[k], henres, kind='linear')
    inter_sigma2 = interp1d(est_en_edges[k], henres2, kind='linear')
    inter_part = interp1d(est_en_edges[k], henpart, kind='linear')

    mean1 = enest - inter_mean1(enest)*enest   #  pEnBiasVsEnEst = ( E' - E ) / E'
    mean2 = enest - inter_mean1(enest)*enest  #  pEnBiasVsEnEst = ( E' - E ) / E'
    sigma1 = inter_sigma1(enest)*enest
    sigma2 = inter_sigma2(enest)*enest
    part   = inter_part(enest)

    print ('mean1=', mean1, ' mean2= ', mean2, ' sigma1= ', sigma1, ' sigma2= ', sigma2, ' part=',part)

    result = part*gaussian(entrue,mean1,sigma1)+(1.-part)*gaussian(entrue,mean2,sigma2)

    return result

#-----------------------------------------------------------------------------------------------------------------------------------------

def DGauss(func, EtMin, EtMax, Eps):
    """Integral implemented by M.Martinez in a Fortran code. The inputs parameters are: the function that has to be integrated (func), the extremes of the integration range (EtMin, EtMax) and the precision required for the integral (Eps)."""
    
    W = np.array([0.1012285362903762591525313543,
                  0.2223810344533744705443559944,
                  0.3137066458778872873379622020,
                  0.3626837833783619829651504493,
                  0.02715245941175409485178057246,
                  0.06225352393864789286284383699,
                  0.09515851168249278480992510760,
                  0.1246289712555338720524762822,
                  0.1495959888165767320815017305,
                  0.1691565193950025381893120790,
                  0.1826034150449235888667636680,
                  0.1894506104550684962853967232])

    X = np.array([0.9602898564975362316835608686,
                  0.7966664774136267395915539365,
                  0.5255324099163289858177390492,
                  0.1834346424956498049394761424,
                  0.9894009349916499325961541735,
                  0.9445750230732325760779884155,
                  0.8656312023878317438804678977,
                  0.7554044083550030338951011948,
                  0.6178762444026437484466717640,
                  0.4580167776572273863424194430,
                  0.2816035507792589132304605015,
                  0.09501250983763744018531933543])

    dgauss = 0.0
    
    if EtMax == EtMin: 
       return dgauss 
    else:
       const = 0.005/(EtMax-EtMin)
       BB = EtMin  
       control = 0        
        
       while 1:
          if (control == 0):
             AA = BB
             BB = EtMax
        
          C1 = 0.05*(BB+AA)
          C2 = 0.05*(BB-AA)
          S8 = 0.0
          S16 = 0.0

          for i in range (0, 12):
           
             if (i<4):
                U = C2*X[i]
                S8 = S8+W[i]*(func(C1+U)+func(C1-U))
             else:
                U = C2*X[i]
                S16 = S16+W[i]*(func(C1+U)+func(C1-U))
          
          S8 = C2*S8
          S16 = C2*S16
              
          if (math.fabs(S16-S8)<=Eps*(1.+math.fabs(S16))):
             dgauss = dgauss+S16
             if BB != EtMax:
                control = 0
                continue
             else:
                print('\n')
                print('The integral is: ', dgauss)
                break
          else:
             BB = C1
             if (1.+math.fabs(const*C2) != 1.):
                control = 1
                continue  #C1
             else:
                dgauss = 0.0
                print('\n')
                print('The integral is: ', dgauss)
                break

      
#-----------------------------------------------------------------------------------------------------------------------------------------
     
def PdfOnFunc():
   fi = 47
   T = 5400   
   funct = lambda E: fi*(E**-5)*T*get_area(E)*get_enmigmat(E) 
   pdfon = DGauss(funct, 35, 2000, 0.5)
   return pdfon

#-----------------------------------------------------------------------------------------------------------------------------------------

Init()
PdfOnFunc()
