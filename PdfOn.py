#!/home/software/anaconda/bin/python

import sys
sys.path.append('/media/san/astro/soft/root/root_v5.34.36/root/lib')
import ROOT
import numpy as np
#import scipy as sp
#from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import math
from tqdm import tqdm 
import config
from scipy.optimize import minimize
from cuts import pass_energy_cuts
#from num_integration import DGauss
#from mpmath import dps

from ROOT import TString, TFile, TObject, TH1, TMath, TNtuple, TCanvas 
from array import array
#from ROOT import *   (this command line does not work in a script)

# all global parameters here
content = 0
ene_on = []
integrand = []
integrand_norm = []
logenergy = []
pEnBiasVsEnEst = []
pEnBiasVsEnEst2 = []
hEnResVsEnEst = []
hEnResVsEnEst2 = []
hEnPartVsEnEst = []
hEnMinVsEnEst = []
hEnMaxVsEnEst = []
nbinsenest = []
numsamples   = config.nsamples
numenergies  = config.energy_numbins
est_en_edges = np.empty(shape=(numsamples, numenergies), dtype='double')

penbias  = np.empty(shape=(numsamples, numenergies), dtype='double')
penbias2 = np.empty(shape=(numsamples, numenergies), dtype='double')
henres   = np.empty(shape=(numsamples, numenergies), dtype='double')
henres2  = np.empty(shape=(numsamples, numenergies), dtype='double')
henpart  = np.empty(shape=(numsamples, numenergies), dtype='double')

normon   = np.empty(numsamples, dtype='double') # content after applying the ON-boundaries in energy


inter_mean1  = []
inter_mean2  = []
inter_sigma1 = []
inter_sigma2 = []
inter_part   = []
enmig        = []

mean1 = []
mean2 = []
sigma1= []
sigma2= []
part  = []

numtoton = np.empty(numsamples, dtype='int')    # total number of events in ON sample (one for each period)

sEffArea  = []
sEffAreaErec = []
AeffEnMax = []

#-----------------------------------------------------------------------------------------------------------------------------------------

def OnInit(PrelimPath,OnDir,EffOnTime):

#    NumBckgEnBins = 40
    print("\n\n--- INITIALIZING VALUES from:", PrelimPath, " eff. on time: ", EffOnTime[0], "\n")

    # TFile *LivParamFile  = TFile.Open(PrelimPath, "READ")
    LivParamFile2 = TFile.Open(PrelimPath, "READ")

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

    global penbias
    global penbias2
    global henres
    global henres2
    global henpart

    global inter_mean1
    global inter_mean2
    global inter_sigma1
    global inter_sigma2
    global inter_part

    global ene_on
    global integrand
    global normon

    for k in range (0, config.nsamples):    # k= number of samples

        pEnBiasVsEnEst.append(LivParamFile2.FindObjectAny(ROOT.Form("pEnBiasVsEnEst__%s_smoothed" % config.SampleName[k])))
        pEnBiasVsEnEst[k].SetDirectory(0)        

        pEnBiasVsEnEst2.append(LivParamFile2.FindObjectAny(ROOT.Form("pEnBiasVsEnEst2__%s_smoothed" % config.SampleName[k])))
        pEnBiasVsEnEst2[k].SetDirectory(0)        

        hEnResVsEnEst.append(LivParamFile2.FindObjectAny(ROOT.Form("hEnResVsEnEst__%s_smoothed" % config.SampleName[k])))
        hEnResVsEnEst[k].SetDirectory(0)        

        hEnResVsEnEst2.append(LivParamFile2.FindObjectAny(ROOT.Form("hEnResVsEnEst2__%s_smoothed" % config.SampleName[k])))
        hEnResVsEnEst2[k].SetDirectory(0)        

        hEnPartVsEnEst.append(LivParamFile2.FindObjectAny(ROOT.Form("hEnPartVsEnEst__%s_smoothed" % config.SampleName[k])))
        hEnPartVsEnEst[k].SetDirectory(0)       
 
        nbinsenest.append(pEnBiasVsEnEst[k].GetNbinsX())

        print("Est Energy Edges for period st%s" % config.SampleName[k], "\n")

        print ("nbinsenest =", nbinsenest[k], "\n")

        for i in range (0, nbinsenest[k]):
            est_en_edges[k][i] = pEnBiasVsEnEst[k].GetBinLowEdge(i+1) 
            print( pEnBiasVsEnEst[k].GetBinLowEdge(i+1), flush = True)
        
        print ("\n")

        sEffArea.append(LivParamFile2.FindObjectAny(ROOT.Form("sEffArea__%s" % config.SampleName[k])))
        sEffAreaErec.append(LivParamFile2.FindObjectAny(ROOT.Form("sEffAreaErec__%s" % config.SampleName[k])))
        

        AeffEnMax.append(sEffArea[k].GetXmax())

    LivParamFile2.Close()
    LivParamFile2.Delete()


    for k in range (0, config.nsamples):    # k= number of samples    

        print (OnDir[k])

        file        = TFile.Open(OnDir[k], "READ")
        ntuple      = file.Get("fOnSample")
        numtoton[k] = ntuple.GetEntries()

        ene_on.append(np.empty(shape=(config.EvtMaxOn),dtype=config.ArrayDataType))

        num = 0
        for i in range (0, numtoton[k]):
            ntuple.GetEntry(i)
            energy = ntuple.GetArgs()[0]

            if (not pass_energy_cuts(energy)):
                continue               

            ene_on[k][num] = energy
            num += 1
                
        numtoton[k] = num
        ene_on[k] = np.resize(ene_on[k],num)

        print ('The number of on entries (after energy cuts) is:', num)
        print (ene_on[k])
        print ('The sum of ', np.sum(ene_on[k]))

        for i in range (0, nbinsenest[k]):

            penbias[k][i]  = pEnBiasVsEnEst[k].GetBinContent(i+1)
            penbias2[k][i] = pEnBiasVsEnEst2[k].GetBinContent(i+1)
            henres[k][i]   = hEnResVsEnEst[k].GetBinContent(i+1)
            henres2[k][i]  = hEnResVsEnEst2[k].GetBinContent(i+1)
            henpart[k][i]  = hEnPartVsEnEst[k].GetBinContent(i+1) 


        inter_mean1.append(interp1d(est_en_edges[k], penbias[k], kind='linear'))
        inter_mean2.append(interp1d(est_en_edges[k], penbias2[k], kind='linear'))
        inter_sigma1.append(interp1d(est_en_edges[k], henres[k], kind='linear'))
        inter_sigma2.append(interp1d(est_en_edges[k], henres2[k], kind='linear'))
        inter_part.append(interp1d(est_en_edges[k], henpart[k], kind='linear'))
            
#        mean1.append(np.empty(shape=(num),dtype='double'))
#        mean2.append(np.empty(shape=(num),dtype='double'))
#        sigma1.append(np.empty(shape=(num),dtype='double'))
#        sigma2.append(np.empty(shape=(num),dtype='double'))
#        part.append(np.empty(shape=(num),dtype='double'))

#        mean1[k] = ene_on[k] - inter_mean1[k](ene_on[k])*ene_on[k]  #  pEnBiasVsEnEst = ( E' - E ) / E'
#        mean2[k] = ene_on[k] - inter_mean1[k](ene_on[k])*ene_on[k]  #  pEnBiasVsEnEst = ( E' - E ) / E'
#        sigma1[k] = inter_sigma1[k](ene_on[k])*ene_on[k]
#        sigma2[k] = inter_sigma2[k](ene_on[k])*ene_on[k]
#        part[k]   = inter_part[k](ene_on[k])

        integrand.append(np.empty(shape=(num,config.BinsIntegrand),dtype=config.ArrayDataType))

        pbar = tqdm(total=num*config.BinsIntegrand)

        for i in range (0, num):

            energy = ene_on[k][i]

            mean1 = energy - inter_mean1[k](energy)*energy   #  pEnBiasVsEnEst = ( E' - E ) / E'
            mean2 =  energy - inter_mean1[k](energy)*energy  #  pEnBiasVsEnEst = ( E' - E ) / E'
            sigma1 = inter_sigma1[k](energy)*energy
            sigma2 = inter_sigma2[k](energy)*energy
            part = inter_part[k](energy)

            # define the integration limits for true energy
            lowl  = mean2 - config.SigmaLimits * sigma2
            highl = mean2 + config.SigmaLimits * sigma2
            if (lowl < config.energy_absmin):
                lowl = config.energy_absmin

            # define the bins in true energy
            binwidth = (highl - lowl) / config.BinsIntegrand

            for j in range(0,config.BinsIntegrand):
                
                entrue = lowl + (j+0.5)*binwidth

                pbar.update(1)

                area = get_area(entrue,k)
                flux = get_phi(entrue)
                enmig = part*gaussian(entrue,mean1,sigma1)+(1.-part)*gaussian(entrue,mean2,sigma2)
            
                integrand[k][i][j] = flux * area * enmig * EffOnTime[k] * binwidth

        pbar.close()

        print ("Starting Normalization\n")

        integrand_norm.append(np.empty(shape=(config.BinsIntegrandNorm,config.BinsIntegrand),dtype=config.ArrayDataType))

        on_en_edges = [] 
        binsize = TMath.Log(config.energy_upper_cut/config.energy_cut)/config.BinsIntegrandNorm

        for ebin in range (0, config.BinsIntegrandNorm+1):
            on_en_edges.append(config.energy_cut * TMath.Exp(binsize * ebin))

        for ebin in range (0, config.BinsIntegrandNorm):

            energy = math.sqrt(on_en_edges[ebin+1]*on_en_edges[ebin])
            
            mean1  =  energy - inter_mean1[k](energy)*energy   #  pEnBiasVsEnEst = ( E' - E ) / E'
            mean2  =  energy - inter_mean1[k](energy)*energy  #  pEnBiasVsEnEst = ( E' - E ) / E'
            sigma1 = inter_sigma1[k](energy)*energy
            sigma2 = inter_sigma2[k](energy)*energy
            part   = inter_part[k](energy)

            # define the integration limits for true energy
            lowl  = mean2 - config.SigmaLimits * sigma2
            highl = mean2 + config.SigmaLimits * sigma2
            if (lowl < config.energy_absmin):
                lowl = config.energy_absmin

            # define the bins in true energy
            binwidth = (highl - lowl) / config.BinsIntegrand

            for j in range(0,config.BinsIntegrand):
                
                entrue = lowl + (j+0.5)*binwidth

                pbar.update(1)

                area = get_area(entrue,k)
                flux = get_phi(entrue)
                enmig = part*gaussian(entrue,mean1,sigma1)+(1.-part)*gaussian(entrue,mean2,sigma2)
            
                integrand_norm[k][ebin][j] = flux * area * enmig * (on_en_edges[ebin+1]-on_en_edges[ebin]) * EffOnTime[k] 

        normon[k] = np.sum(integrand_norm[k])
        
        print ("Normalization=", normon[k])

#        enmig.append(part*gaussian(entrue,mean1,sigma1)+(1.-part)*gaussian(entrue,mean2,sigma2))

#        print ('mean1=', mean1, ' mean2= ', mean2, ' sigma1= ', sigma1, ' sigma2= ', sigma2, ' part=',part)



    print("--- INITIALIZATION FINISHED!\n")
    return;

#-----------------------------------------------------------------------------------------------------------------------------------------
def gaussian(x, mu, sigma):
          return 1./sigma/math.sqrt(2*math.pi)*np.exp(-0.5*np.power((x - mu)/sigma, 2.))


def get_phi(entrue):

    return config.Phi * (entrue**-5)
 
def get_area(entrue,k):

    global sEffArea

    area = sEffArea[k].Eval(math.log10(entrue))
#    print ('log energy: ',math.log10(entrue), 'area: ', area)
    return math.pow(10.,area)

def get_enmigmat(entrue,k,i):

    global enest    
    global penbias
    global penbias2
    global henres
    global henres2
    global henpart

    mean10 = ene_on[k][i] - inter_mean1[k](ene_on[k][i])*ene_on[k][i]   #  pEnBiasVsEnEst = ( E' - E ) / E'
    mean20= ene_on[k][i] - inter_mean1[k](ene_on[k][i])*ene_on[k][i]  #  pEnBiasVsEnEst = ( E' - E ) / E'
    sigma10 = inter_sigma1[k](ene_on[k][i])*ene_on[k][i]
    sigma20 = inter_sigma2[k](ene_on[k][i])*ene_on[k][i]
    part0 = inter_part[k](ene_on[k][i])

#    print ('entrue=', entrue)
#    print ('mean1=', mean1, ' mean2= ', mean2, ' sigma1= ', sigma1, ' sigma2= ', sigma2, ' part=',part)
 
    result = part0*gaussian(entrue,mean10,sigma10)+(1.-part0)*gaussian(entrue,mean20,sigma20)

#    print ('result=',result)

    return result


#-----------------------------------------------------------------------------------------------------------------------------------------

#   f(x,k,fi,T):

#   x = np.empty(len(part[k])) 
#   return fi * T * get_area(x,k) * (part[k]*gaussian(x,mean1[k],sigma1[k])+(1.-part[k])*gaussian(x,mean2[k],sigma2[k]))

    

def get_integral_etrue(entrue, k, fi, T, i):

    area = get_area(entrue,k)
    flux = fi * (entrue**-5)

    result = fi * T * get_area(entrue,k) #* get_enmigmat(entrue,k,i)

    #get_enmigmat(entrue, k)
#    print ("rsul=",result)
    return result
 
      
#-----------------------------------------------------------------------------------------------------------------------------------------
     
def PdfOnFunc():
   fi = 47.
   T = 5400.
   pdfon = 0
   print ("HERE\n")

#   f = lambda x,kk,ffi,TT: ffi * TT * get_enmigmat(x,kk)

#   print ("TEST: ",np.sum(get_integral_etrue(1000.,0,fi,T)))
#   print ("TEST1:", get_integral_etrue(1000.,0,fi,T)[0], "\n", get_integral_etrue(1000.,0,fi,T)[1], "\n", get_integral_etrue(1000.,0,fi,T)[2], "\n", get_integral_etrue(1000.,0,fi,T)[3], "\n", get_integral_etrue(1000.,0,fi,T)[4], "\n", )
#   print ("TEST2: ",[quad(f,35.,2000.,args=(i,0,fi,T,)) for i in range(0,numtoton[0])])
#   print ("TEST2: ",quad(f,35.,2000.,args=(0,0,fi,T)))
#   print ("TEST2: ",np.sum(get_integral_etrue(0,fi,T)))
#   print ("TEST2: ",quad(f,35.,2000.,args=(0,fi,T))[0])

   for k in range (0, 2):   
#      funct = lambda Etrue: fi*(Etrue**-5)*T*get_area(Etrue,k)*get_enmigmat(Etrue, k) 
#      pdfon = DGauss(funct, 35, 2000, 0.5)
#      mp.dps = 15
       
       pbar = tqdm(total=numtoton[k])
       res  = []        

       for i in range (0, numtoton[k]):
           pbar.update(1)
#           res = quad(get_integral_etrue,35., 2000.,args=(k,fi,T,i))
#           pdfon += math.log(quad(get_integral_etrue,35., 2000.,args=(k,fi,T,i))[0])    #it should be a np array of integrals!
           pdfon +=  math.log(np.sum(integrand[k][i]))
           

       pbar.close()
   return pdfon

#-----------------------------------------------------------------------------------------------------------------------------------------

#Init()
#PdfOnFunc()
