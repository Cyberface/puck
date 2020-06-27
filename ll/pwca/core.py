

#
import pickle
from positive import alert,magenta
from numpy import loadtxt,load
from os.path import exists

# Always load catalog list for calibration runs 
pwca_catalog_path = '/Users/book/KOALA/puck/ll/data/pwca_catalog.pickle'
pwca_catalog = pickle.load( open( pwca_catalog_path, "rb" ) )
alert('Catalog of calibration runs stored to %s'%magenta('"pwca.pwca_catalog"'),fname='pwca.core')

# Always load curated metadata for calibration runs 
metadata_dict_path = '/Users/book/KOALA/puck/ll/data/metadata_dict.pickle'
metadata_dict = load(metadata_dict_path)
alert('Metadata dictionary for calibration runs stored to %s'%magenta('"pwca.metadata_dict"'),fname='pwca.core')

#


# Function to determine version2 data fitting region
def determine_data_fitting_region( data, fmin=0.03, fmax=0.12 ):
    '''
    Given version2 data array, determine fitting region dynamically.
    This function assumes clean data within the default or input values of fmin and fmax, and then uses the phase derivate to determine new fmin and fmax values that are ultimately used to define a fitting region.
    '''
    
    # Import usefuls
    from numpy import argmin,log
    from positive import smooth,find
    
    # Extract data 
    f,amp_td,amp_fd,dphi_td,dphi_fd = data

    # Use default domain bounds to determine a mask
    mask = (f>=fmin) & (f<=fmax)
    
    # Determine the minimum dphi
    # Smooth dphi using postiive's savgol filter
    x = log(f[mask])
    y = smooth(dphi_td[mask]).answer
    knot = argmin(y)
    y_knot = y[knot]
    
    # Determine new fmin and max using heuristic 
    f_knot = f[mask][knot]
    new_fmin = f_knot * 0.4# 0.325
    new_fmax = f_knot + 0.020 # 0.025 
    
    #
    new_mask = (f>=new_fmin) & (f<=new_fmax)
    new_data = data.T[new_mask,:]
    new_data[:,-2:] -= y_knot
    
    #
    new_knot = find(f>=fmin)[0]+knot
    
    #
    return new_data,new_knot,new_fmin,new_fmax,f_knot
    

#
def select_scenty_metadata( sceo ):
    
    '''
    Given nrutils' scentry object, collect metedata useful for generating model waveforms 
    '''
    
    #
    from numpy.linalg import norm
    from numpy import arccos,dot,pi,array
    from positive.physics import calc_chi_eff,calc_chi_p
    
    #
    a = sceo
    
    #
    l = a.L/norm(a.L)
    if (abs(a.m1-a.m2)<1e-3) and (norm(a.X1)<norm(a.X2)):
        a.X1,a.X2 = [ array(k) for k in (a.X2,a.X1) ]
        a.m1,a.m2 = [ float(k) for k in (a.m2,a.m1) ]

    #
    m1,m2 = [ k/(a.m1+a.m2) for k in (a.m1,a.m2) ] 
    eta = m1*m2/(m1+m2)
    
    #
    X1,X2,L,S = a.X1,a.X2,a.L,a.S
    
    #
    l = L/norm(L)
    s = S/norm(S)
    
    #
    theta = arccos( dot( l, s ) ) * 180/pi
    
    #
    chi1 = dot(X1,l)
    chi2 = dot(X2,l)
    
    #
    chi_p   = calc_chi_p(   m1,X1, m2,X2, L )
    chi_eff = calc_chi_eff( m1,X1, m2,X2, L )
    
    #
    delta = (m1-m2)/(m1+m2)
    
    #
    return theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2 