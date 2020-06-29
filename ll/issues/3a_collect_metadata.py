#!/usr/bin/env python2.7

# Setup python environment
from numpy import *
from positive import *
from nrutils import scsearch
from glob import glob
from numpy.linalg import norm
from pwca import determine_data_fitting_region,pwca_catalog
import pickle
import pwca

#
package_dir = parent( pwca.__path__[0] )
data_dir = package_dir + 'data/version2/'
files = glob( datadir+'q*.txt' )

#
alert('We have %s files to consider.'%red(str(len(files))))

#
metadata = []
simnames = []

#
alert('Using %s and txt files at %s to gather metadata...'%(magenta('pwca_catalog'),magenta(datadir)))
for f in files:

    #
    file_name = f.split('/')[-1].split('.')[0]

    #
    A = scsearch( keyword=file_name, verbose=not True, catalog=pwca_catalog )

    #
    a = A[0]
    
    #
    l = a.L/norm(a.L)
    if ('q1a' in file_name) and (norm(a.X1)<norm(a.X2)):
        warning(magenta('Flipping 1-2 labels in euqal mass case'),fname=file_name)
        a.X1,a.X2 = [ array(k) for k in (a.X2,a.X1) ]
        a.m1,a.m2 = [ float(k) for k in (a.m2,a.m1) ]
        print '\t * chi1_l = ',dot(a.X1,l)
        print '\t * chi2_l = ',dot(a.X2,l)

    #
    m1,m2 = [ k/(a.m1+a.m2) for k in (a.m1,a.m2) ] 
    eta = m1*m2/(m1+m2)
    
    #
    X1,X2,L,S = a.X1,a.X2,a.L,a.S
    
    #
    a1,a2 = norm(a.X1),norm(a.X2)
    
    #
    l = L/norm(L)
    s = S/norm(S)
    
    # NOTE that theta is in radians
    theta = arccos( dot( l, s ) ) 
    
    #
    chi1 = dot(X1,l)
    chi2 = dot(X2,l)
    
    #
    chi_p   = calc_chi_p(   m1,X1, m2,X2, L )
    chi_eff = calc_chi_eff( m1,X1, m2,X2, L )
    
    #
    delta = (m1-m2)/(m1+m2)
    
    #
    simnames.append(file_name)
    metadata.append( [ theta,
                       m1,
                       m2,
                       eta,
                       delta,
                       chi_eff,
                       chi_p,
                       chi1,
                       chi2,
                       a1,
                       a2 ] )

#
print 'Done.'

#
metadata_array = array(metadata)

#
keys = [ 'theta', 'm1', 'm2', 'eta', 'delta', 'chi_eff', 'chi_p', 'chi1', 'chi2' ]

#
metadata_dict = {}
metadata_dict = { keys[k]: metadata_array[:,k] for k in range(len(keys)) }
metadata_dict['simname'] = simnames

#
metadata_dict['array_data_columns'] = keys
metadata_dict['array_data'] = metadata_array

#
metadata_path = '/Users/book/KOALA/puck/ll/data/metadata_dict.pickle'
alert('Saving metadata dictionary to %s'%magenta(metadata_path))
pickle.dump( metadata_dict, open( metadata_path, "wb" ) )

#
alert('Saving complete.')
    
        

