#!/home/software/anaconda/bin/python

import sys
sys.path.append('/media/san/astro/soft/root/root_v5.34.36/root/lib')
import ROOT
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d
from scipy.integrate import quad
import matplotlib.pyplot as plt
import math
import config
import time
from cuts import pass_energy_cuts
from scipy.optimize import minimize

from ROOT import TFile, TObject, TH1, TMath, TNtuple, TCanvas, TLine, TColor, TPad, TLatex
from array import array
#from ROOT import *   (this command line does not work in a script)

# all global parameters here
numsamples  = config.nsamples
numtotoff   = np.empty(numsamples, dtype='int')    # total number of events in OFF sample (one for each period)
normoff     = np.empty(numsamples, dtype='double') # content after applying the ON-boundaries in energy
log_ene_off = []
inter_func = []

#-----------------------------------------------------------------------------------------------------------------------------------------

#
# this function is deprecated
#
def off_energy_distr(OffDir, NumBckgEnBins=40, EnCutHighBG=8000, EnCutLowBG=100):
    """ 
    this function does:

    - read the OFF sample input files
    - 


    """

    global normoff
    global numtotoff
    global numsamples

    for k in range (0,numsamples):

        print("\n\n--- INITIALIZING VALUES from:", OffDir[k], "\n")
        DoHisto = TFile.Open(OffDir[k], "READ")
        print(DoHisto)

        noff = np.array([ DoHisto.Get("fOffSample") ])

        numtotoff[k] = noff[0].GetEntries()
        print("The total number of entries in the normalization histogram is:", numtotoff[k])    
        
        off_en_edges = [] 
        binsize = TMath.Log(EnCutHighBG/EnCutLowBG)/config.NumBckgEnBins    

        for ebin in range (0, config.NumBckgEnBins+1):
            off_en_edges.append(EnCutLowBG * TMath.Exp(binsize * ebin))

        # create here offbins with the (larger) background limits
        offbins = array('d',off_en_edges)    

        isto = ROOT.TH1F("integral", "integral", config.NumBckgEnBins, offbins)  

        for i in range (0, numtotoff[k]):
            noff[0].GetEntry(i)
            energy = noff[0].GetArgs()[0]
            isto.Fill(energy)

        f = TFile('normalization_histogram_sample_{}.root'.format(k), "RECREATE")
        gramma = TCanvas("EOff_Histogram", 'EOFF_Histogram sample {}'.format(k))
        gramma.SetCanvasSize(800, 800)
        gramma.SetWindowSize(1000, 1000)
        gramma.SetLogx()
        gramma.SetLogy()
        isto.SetXTitle("Log Eest")
        isto.SetYTitle("Log events")
        isto.Draw()

        ll = isto.GetXaxis().FindBin(EnCutLowBG)
        hh = isto.GetXaxis().FindBin(EnCutHighBG)

        line = TLine()
        line.SetLineColor(TColor.kRed)
        line.DrawLine(ll,isto.GetMinimum(),ll,isto.GetMaximum())
        line.DrawLine(hh,isto.GetMinimum(),hh,isto.GetMaximum())

#        gPad.Update()
   
        normoff[k] = isto.Integral(ll, hh)
        print ("The integral from ",ll," to ",hh," of the normalization histogram is:", normoff[k]) 

        latex = TLatex()
        latex.DrawLatex(ll,isto.GetMaximum()/2,"The integral of the normalization histogram is: {}".format(normoff[k]))

        gramma.Write()
        f.Close()
    return;
    
#-----------------------------------------------------------------------------------------------------------------------------------------


def get_OffPdf(energy,nsample):
    return math.exp(inter_func[nsample](math.log(energy)))

def off_normalization(plot=False):

    global normoff

    for k in range (0,numsamples):

        normoff[k], errk = quad(get_OffPdf, config.energy_cut,config.energy_upper_cut, args=(k,)) 
        print ('Background Norm sample=',k,' is: ',normoff[k], ' +- ',errk)


#-----------------------------------------------------------------------------------------------------------------------------------------

def test_sum(k):

    global numtotoff
    global log_ene_off
    global inter_func

    summ  = 0.
    for i in range (0,numtotoff[k]):
        if (np.isinf(summ)):
            print ('i=',i-1,'e=',log_ene_off[k][i-1],'v=',inter_func(log_ene_off[k][i-1]))
            break
        else:
            summ += inter_func[k](log_ene_off[k][i])
            
    return summ

def create_off_interpolations(OffDir,plot=False):

    global log_ene_off
    global numtotoff
    global inter_func

    off_en_edges = [] 
    off_en_center = []
    log_off_en_center   = np.empty(shape=(numsamples,config.NumBckgEnBins),dtype='double')
    y_log_off_en_center = np.empty(shape=(numsamples,config.NumBckgEnBins),dtype='double')

    binsize = math.log(config.EnCutHighBG/config.EnCutLowBG)/config.NumBckgEnBins    

    for ebin in range (0, config.NumBckgEnBins + 1):
       off_en_edges.append(config.EnCutLowBG * math.exp(binsize * ebin))

    for ebin in range (0, config.NumBckgEnBins):
        off_en_center.append(math.sqrt(off_en_edges[ebin+1]*off_en_edges[ebin]))

    # create here offbins with the (larger) background limits
    offbins = array('d',off_en_edges)   

    for k in range (0,numsamples):

        print("\n\n--- INITIALIZING VALUES from:", OffDir[k], "\n")
        DoHisto = TFile.Open(OffDir[k], "READ")
        print(DoHisto)

        log_ene_off.append(np.empty(shape=(config.EvtMaxOff),dtype='double'))  # array of log(energies) from all OFF events
        noff = np.array([ DoHisto.Get("fOffSample") ])
        numtotoff[k] = noff[0].GetEntries()
        print("The total number of entries in the normalization histogram is:", numtotoff[k])    
        
        isto = ROOT.TH1F("interpolation", "interpolation", config.NumBckgEnBins, offbins)

        num = 0
        norm = 0.

        for i in range (0, numtotoff[k]):
            noff[0].GetEntry(i)
            energy = noff[0].GetArgs()[0]
            isto.Fill(energy)

            if (not pass_energy_cuts(energy)):
                continue               

            log_ene_off[k][num] = math.log(energy)
            num += 1

        numtotoff[k] = num
        log_ene_off[k] = np.resize(log_ene_off[k],num)

        print ('The number of off entries (after energy cuts) is:', num)
        print (log_ene_off[k])
        print ('The sum of ', np.sum(log_ene_off[k]))

        ebin_cnt = 0
        for ebin in range (0, config.NumBckgEnBins):
            bin = isto.GetXaxis().FindBin(off_en_center[ebin])
            bin_content = isto.GetBinContent(bin)
            if (bin_content < 1.):
                print ('ebin: ',ebin,' bin: ', bin, 'skipped because of zero entries')
                continue

            # calculate here the logarithm of dN/dE (which normalized will provide dP/dE)
            log_off_en_center[k][ebin_cnt]   = math.log(off_en_center[ebin])
            y_log_off_en_center[k][ebin_cnt] = math.log(bin_content / (off_en_edges[ebin+1]-off_en_edges[ebin]))
            print ('ebin: ',ebin,' ebin_cnt=', ebin_cnt, ' log_off_center: ',log_off_en_center[k][ebin_cnt], ' content: ',y_log_off_en_center[k][ebin_cnt], 'bin= ',bin, ' content=', bin_content)
            ebin_cnt += 1

        inter_func.append(interp1d(log_off_en_center[k], y_log_off_en_center[k], kind='linear'))    #inter_func = interp1d(x, y, kind='linear') 
        print(math.exp(inter_func[k](math.log(90))))

        print ('found sum: ',np.sum(inter_func[k](log_ene_off[k])))
    
        if plot:
            xnew = np.linspace(TMath.Log(config.energy_cut), TMath.Log(config.energy_upper_cut), num=200, endpoint=True)
            plt.plot(log_off_en_center[k], y_log_off_en_center[k], 'o', xnew, inter_func[k](xnew), '-')
            plt.xlabel('Log Eest')
            plt.ylabel('Log events')
            plt.legend(['data', 'linear'], loc='best')
            plt.show()
            plt.savefig('interpolation_{}.png'.format(k))
    
        
#       for k in range (0, 2):    #k in range 0,2 ???
##            if k < NumSubSamples:
#               hOffHeight[k] = ROOT.TH1F(ROOT.Form("hOffHeight_%s" % SampleName[k]), ROOT.Form("Background Level_%s" % SampleName[k]), config.NumBckgEnBins, off_en_edges)
##            else:
##              hOffHeight[k] = ROOT.TH1F("hOffHeight", "Background Level", config.NumBckgEnBins, off_en_edges)

##            if k < NumSubSamples:
##               hOnHeight[k] = ROOT.TH1F(ROOT.Form("hOnHeight_%s" % SampleName[k]), ROOT.Form("Signal Level_%s" % SampleName[k]), config.NumBckgEnBins, off_en_edges)
##            else:
##               hOnHeight[k] = ROOT.TH1F("hOnHeight", "Signal Level", config.NumBckgEnBins, off_en_edges)
    return;  

#-----------------------------------------------------------------------------------------------------------------------------------------

def minim(k):

    mini = sp.optimize.minimize(log_like, 10000, args=(k),method='Powell') 

    b_est = mini.x
    y_0 = log_like(b_est,k)
 
    print ('b_0: ',b_est,'y_0: ', y_0)

    sigma = 2

    root_fn = lambda b: log_like(b,k) - y_0 - 2*sigma

    b_low  = sp.optimize.brentq(root_fn, 1., b_est)
    b_high = sp.optimize.brentq(root_fn, b_est, 1e10)

    return b_est,b_low,b_high;   

#-----------------------------------------------------------------------------------------------------------------------------------------

def log_like(b_k,k):
#    summa = np.float128(0.)

    if (b_k < 0.):
        return 999999999.

    summa = np.sum(inter_func[k](log_ene_off[k]))
    y = -2.* (summa + numtotoff[k]*math.log(b_k) - b_k)
    print ('b_k: ',b_k, 'y: ',y)
    return y; 

#-----------------------------------------------------------------------------------------------------------------------------------------

#isto()
#interpolate()
#y = minim()
#print (y)
