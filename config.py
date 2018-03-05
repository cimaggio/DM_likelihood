#!/home/software/anaconda/bin/python

nsamples        = 2
SampleName = ["0306", "0307"]


# reconstructed energy limits
energy_cut      = 70
energy_upper_cut= 6500
energy_numbins  = 17    # from CreateAeff.C
energy_absmin   = 10  

# Background binning
NumBckgEnBins = 40
EnCutHighBG   = 1.15*energy_upper_cut  
EnCutLowBG    = 0.9 *energy_cut

# Initialization of array sizes 
EvtMaxOn      = 540000
EvtMaxOff     = 302000
BinsIntegrand = 5
BinsIntegrandNorm = 100
SigmaLimits   =  4

# Initialization of array types
ArrayDataType = 'double'

# Flux initialization
Phi    = 4500.    # should be in m^-2 s^-1 


