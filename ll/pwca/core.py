

#
import pickle
from positive import alert,magenta

# Always load metadata for calibration runs 
metadata_path = '/Users/book/KOALA/puck/ll/data/pwca_catalog.sce'
pwca_catalog = pickle.load( open( metadata_path, "rb" ) )
alert('Metadata for calibration runs stored to %s'%magenta('"pwca.pwca_catalog"'),fname='pwca.core')

# Function to collect metadata for run 
def collect_file_metadata(f,catalog=None):
    
    #
    import pickle
    
    #
    file_name = f.split('/')[-1].split('.')[0]
    
    #
    if catalog is None:
        catalog = pwca_catalog
    
    #
    A = scsearch( keyword=file_name, verbose=True, catalog=catalog )
    
    #
    a = A[0]
    
    #
    m1,m2 = a.m1,a.m2 
    eta = m1*m2/(m1+m2)
    theta = None
    
    #
    return A[0]

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
    new_fmin = f[mask][knot] * 0.325
    new_fmax = f[mask][knot] + 0.025 
    
    #
    new_mask = (f>=new_fmin) & (f<=new_fmax)
    new_data = data.T[new_mask,:]
    
    #
    new_knot = find(f>=fmin)[0]+knot
    
    #
    return new_data,new_knot,new_fmin,new_fmax