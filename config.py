#!/home/software/anaconda/bin/python

nsamples        = 2
SampleName = ["0306", "0307"]

GlobalDir = '/media/san2/astro/M15/'
ModelsDir = GlobalDir+'masses'

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
BinsIntegrand = 3
BinsIntegrandNorm = 10
SigmaLimits   =  4

# Initialization of array types
ArrayDataType = 'double'

# Flux initialization
Mode    = 'ann'       # possible modes: ann or dec
ModelInterpolation = 'linear'
Jfactor = 1.4e25      # J-factor from (H.E.S.S. for M15)
tauDMm1 = 4.1e17      # age of Universe in seconds
refsv   = 3.0e-26     # thermal relic cross-section for reference

