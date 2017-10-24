#!/home/software/anaconda/bin/python

import sys
sys.path.append('/media/san/astro/soft/root/root_v5.34.36/root/lib')
import ROOT
import numpy as np
#print (ROOT.TTimeStamp().AsString())
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
OffDir06="/media/san2/astro/M15/st0306/flute_Idx0/Samples.root"
OffDir07="/media/san2/astro/M15/st0307/flute_Idx0/Output_glike.root" #this file does not exist!!! 

#-----------------------------------------------------------------------------------------------------------------------------------------

def isto():
    """ 
    this does that:
    """

    print("\n\n--- INITIALIZING VALUES from:", OffDir06, "\n")

    DoHisto06 = TFile.Open(OffDir06, "READ")
    DoHisto07 = TFile.Open(OffDir07, "READ")
    noff = np.array([DoHisto06.Get("fOffSample"), DoHisto07.Get("fOffSample")])
    global content
    global numtot
    Samplename = np.array([6, 7])   
    for k in range (0,2):
       numtot = noff[k].GetEntries()
       print("The total number of entries in the normalization histogram is:", numtot[k])    

       NumBckgEnBins = 40
       EnCutHighBG = 8000  #55435.5   #to change with the real value #GeV
       EnCutLowBG = 100  #5.54355 o 300
       off_en_edges = [] 
        
       binsize = TMath.Log(EnCutHighBG/EnCutLowBG)/NumBckgEnBins    

       for ebin in range (0, NumBckgEnBins+1):
          off_en_edges.append(EnCutLowBG * TMath.Exp(binsize * ebin))

       offbins = array('d',off_en_edges)    

       isto = ROOT.TH1F("integral", "integral", NumBckgEnBins, offbins)  

       for i in range (0, numtot[k]):
          noff[k].GetEntry(i)
          energy = noff[k].GetArgs()[0]
          isto.Fill(energy[k])

       f = TFile("normalization_histogram_st030 %d .root", Samplename[k], "RECREATE")
       gramma = TCanvas("EOff_Histogram", "EOFF_Histogram")
       gramma.SetCanvasSize(800, 800)
       gramma.SetWindowSize(1000, 1000)
       gramma.SetLogx()
       gramma.SetLogy()
       isto[k].SetXTitle("Log Eest")
       isto[k].SetYTitle("Log events")
       isto[k].Draw()
       gramma.Write()
       f.Close()

       ll = isto[k].GetXaxis().FindBin(EnCutLowBG)
       hh = isto[k].GetXaxis().FindBin(EnCutHighBG)
   
       content = isto[k].Integral(ll, hh)
       print ("The integral of the normalization histogram is:", content[k]) 

#    for ebin in range (0, NumBckgEnBins):
#       off_en_center [ebin] = TMath.Sqrt(off_en_edges [ebin+1]*off_en_edges [ebin])
##      off_en_center[ebin] = (off_en_edges[ebin+1]+off_en_edges[ebin])/2.
#       log_off_en_center [ebin] = TMath.Log( off_en_center [ebin])
        
#       for k in range (0, 2):    #k in range 0,2 ???
##            if k < NumSubSamples:
#               hOffHeight[k] = ROOT.TH1F(ROOT.Form("hOffHeight_%s" % SampleName[k]), ROOT.Form("Background Level_%s" % SampleName[k]), NumBckgEnBins, off_en_edges)
##            else:
##              hOffHeight[k] = ROOT.TH1F("hOffHeight", "Background Level", NumBckgEnBins, off_en_edges)

##            if k < NumSubSamples:
##               hOnHeight[k] = ROOT.TH1F(ROOT.Form("hOnHeight_%s" % SampleName[k]), ROOT.Form("Signal Level_%s" % SampleName[k]), NumBckgEnBins, off_en_edges)
##            else:
##               hOnHeight[k] = ROOT.TH1F("hOnHeight", "Signal Level", NumBckgEnBins, off_en_edges)

    return;

    
#-----------------------------------------------------------------------------------------------------------------------------------------

def interpolate(plot=False):

    DoHisto06 = TFile.Open(OffDir06, "READ")
    noff = DoHisto06.Get("fOffSample")

    global ene
    global numtot
    global inter_func

    NumBckgEnBins = 40
    EnCutHighBG = 1.15*energy_upper_cut  
    EnCutLowBG =  0.85*energy_cut
    off_en_edges = [] 
    off_en_center = []
    log_off_en_center = []
    bin_off_en_center = []
    y_log_off_en_center = []    
    binsize = TMath.Log(EnCutHighBG/EnCutLowBG)/NumBckgEnBins    

    for ebin in range (0, NumBckgEnBins + 1):
       off_en_edges.append(EnCutLowBG * TMath.Exp(binsize * ebin))

    offbins = array('d',off_en_edges)   
    isto = ROOT.TH1F("interpolation", "interpolation", NumBckgEnBins, offbins)

    num = 0

    for i in range (0, numtot):
       noff.GetEntry(i)
       energy = noff.GetArgs()[0]
       isto.Fill(energy)
       if (energy < energy_cut):
           continue
       if (energy > energy_upper_cut):
           continue
           
       ene[num] = math.log(energy)
       num += 1

    numtot = num
    ene = np.resize(ene,num)

    print ('The number of off entries (after energy cuts) is:', num)
    print (ene)
    print ('The sum of ', np.sum(ene))

    ebin_cnt = 0
    for ebin in range (0, NumBckgEnBins):
       off_en_center.append(TMath.Sqrt(off_en_edges[ebin+1]*off_en_edges[ebin]))
       bin = isto.GetXaxis().FindBin(off_en_center[ebin])
       bin_off_en_center.append(bin)

       bin_content = isto.GetBinContent(bin)
       if (bin_content < 1.):
           continue

       log_off_en_center.append(TMath.Log(off_en_center[ebin]))
       y_log_off_en_center.append(TMath.Log(bin_content))
       print ('ebin: ',ebin,' log_off_center: ',log_off_en_center[ebin_cnt], ' content: ',y_log_off_en_center[ebin_cnt], 'bin= ',bin_content)
       ebin_cnt += 1

    inter_func = interp1d(log_off_en_center, y_log_off_en_center, kind='linear')    #inter_func = interp1d(x, y, kind='linear') 

#    print (np.sum(inter_func(ene)))

    summ = 0.
    for i in range (0,num):
        if (np.isinf(summ)):
            print ('i=',i-1,'e=',ene[i-1],'v=',inter_func(ene[i-1]))
            break
        else:
            summ += inter_func(ene[i])

    print ('found sum: ', summ)
    
    if plot:
        xnew = np.linspace(TMath.Log(80), TMath.Log(7000), num=200, endpoint=True)
        plt.plot(log_off_en_center, y_log_off_en_center, 'o', xnew, inter_func(xnew), '-')
        plt.xlabel('Log Eest')
        plt.ylabel('Log events')
        plt.legend(['data', 'linear'], loc='best')
        plt.show()
        plt.savefig('interpolation.png')
    
        
#       for k in range (0, 2):    #k in range 0,2 ???
##            if k < NumSubSamples:
#               hOffHeight[k] = ROOT.TH1F(ROOT.Form("hOffHeight_%s" % SampleName[k]), ROOT.Form("Background Level_%s" % SampleName[k]), NumBckgEnBins, off_en_edges)
##            else:
##              hOffHeight[k] = ROOT.TH1F("hOffHeight", "Background Level", NumBckgEnBins, off_en_edges)

##            if k < NumSubSamples:
##               hOnHeight[k] = ROOT.TH1F(ROOT.Form("hOnHeight_%s" % SampleName[k]), ROOT.Form("Signal Level_%s" % SampleName[k]), NumBckgEnBins, off_en_edges)
##            else:
##               hOnHeight[k] = ROOT.TH1F("hOnHeight", "Signal Level", NumBckgEnBins, off_en_edges)
    return;  

#-----------------------------------------------------------------------------------------------------------------------------------------

def minim():
    mini = sp.optimize.minimize(log_like, 10000, method='Powell') 

    mu_est = mini.x
    y_0 = log_like(mu_est)
 
    print ('mu_0: ',mu_est,'y_0: ', y_0)

    root_fn = lambda mu: log_like(mu) - y_0 - 4.0

    mu_low  = sp.optimize.brentq(root_fn, 1., mu_est)
    mu_high = sp.optimize.brentq(root_fn, mu_est, 1e10)

    return mini,mu_low,mu_high;   

#-----------------------------------------------------------------------------------------------------------------------------------------

def log_like(mu):
#    summa = np.float128(0.)

    if (mu < 0.):
        return 999999999.


    summa = np.sum(inter_func(ene))
    y = -2.* (summa + numtot*math.log(mu) - mu)
    print ('mu: ',mu, 'y: ',y)
    return y; 

#-----------------------------------------------------------------------------------------------------------------------------------------

isto()
interpolate()
y = minim()
print (y)
