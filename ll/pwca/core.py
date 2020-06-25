

#
import pickle
from positive import alert,magenta
from numpy import load
from os.path import exists

# Always load catalog list for calibration runs 
pwca_catalog_path = '/Users/book/KOALA/puck/ll/data/pwca_catalog.pickle'
pwca_catalog = pickle.load( open( pwca_catalog_path, "rb" ) )
alert('Catalog of calibration runs stored to %s'%magenta('"pwca.pwca_catalog"'),fname='pwca.core')

# Always load curated metadata for calibration runs 
metadata_dict_path = '/Users/book/KOALA/puck/ll/data/metadata_dict.pickle'
metadata_dict = load(metadata_dict_path)
alert('Metadata dictionary for calibration runs stored to %s'%magenta('"pwca.metadata_dict"'),fname='pwca.core')

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
    
    # Determine new fmin and max using heuristic 
    f_knot = f[mask][knot]
    new_fmin = f_knot * 0.325
    new_fmax = f_knot + 0.025 
    
    #
    new_mask = (f>=new_fmin) & (f<=new_fmax)
    new_data = data.T[new_mask,:]
    
    #
    new_knot = find(f>=fmin)[0]+knot
    
    #
    return new_data,new_knot,new_fmin,new_fmax,f_knot