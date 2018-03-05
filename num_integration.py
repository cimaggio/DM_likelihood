#
# Library for numerical integration
#
import numpy as np
from math import fabs, pow

def DGauss(func, EtMin, EtMax, Eps=1.0e-5):
    """Integral implemented by M.Martinez in a Fortran code. 
       The inputs parameters are: 
       - the function that has to be integrated (func), 
       - the extremes of the integration range (EtMin, EtMax) 
       - and the precision required for the integral (Eps).

       Does not work with plain or (partly) linear functions!

    """

    numprecision = 15
    precref = pow(10.,-numprecision)
    if (Eps < precref):
        print ('Required precision of', Eps, ' cannot be accomplished, will switch to default: ',precref)
        Eps = precref

    
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

    dgauss = np.float128(0.)
    
    if EtMax == EtMin: 
        return dgauss 
    else:
        const = 0.005/(EtMax-EtMin)
        BB    = EtMin
        control = 0        
        
        while 1:
            if (control == 0):
                AA = BB
                BB = EtMax
                
            C1 = 0.5*(BB+AA)
            C2 = 0.5*(BB-AA)
            S8  = 0.0
            S16 = 0.0

            for i in range (0, 12):
           
                if (i<4):
                    U = C2*X[i]
                    S8 = S8+W[i]*(func(C1+U)+func(C1-U))
                    print ('S8: ',S8,'      C1+U=',C1+U,' C1-U=',C1-U, ' func1:',func(C1+U), 'func2:',func(C1-U))
                else:
                    U = C2*X[i]
                    S16 = S16+W[i]*(func(C1+U)+func(C1-U))
                    print ('S16: ',S16,'     C1+U=',C1+U,' C1-U=',C1-U, ' func:',func(C1+U), 'func2:',func(C1-U))
          
            S8 = C2*S8
            S16 = C2*S16

 #           print ('AA=',AA,' BB=',BB, 'S8=',S8, 'S16=',S16,' dgauss=',dgauss, ' S16-18=', S16-S8)
              
            if (fabs(S16-S8)<=Eps*(fabs(S16))):
                dgauss = dgauss+S16

#                print ('\nBB=',BB,'\n')
#                print ('CHECK3=',round((BB-EtMax)/EtMax,numprecision),'\n')
                if (round(fabs((BB-EtMax)/EtMax),numprecision) != 0.0):
                    control = 0
                    continue
                else:
                    print('\n')
                    print('The integral is: ', dgauss)
                    break
            else:
                BB = C1
#                print ('\nCHECK=',fabs(const*C2))
#                print ('CHECK2=',round(const*C2,numprecision),'\n')
                if (round(fabs(const*C2),numprecision) != 0.0):
                    control = 1
                    continue  #C1
                else:
                    print('\n')
                    print('The integral is: ', dgauss)
                    break

        return dgauss
