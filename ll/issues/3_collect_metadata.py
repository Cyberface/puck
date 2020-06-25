#!/usr/bin/env python2.7

# Setup python environment
from numpy import *
from positive import *
from nrutils import scsearch
from glob import glob
from numpy.linalg import norm
from pwca import determine_data_fitting_region,pwca_catalog

#
datadir = '/Users/book/KOALA/puck/ll/data/version2/'
files = glob( datadir+'*.txt' )

#
metadata = []

#
alert('Using %s and txt files at %s to gather metadata...'%(magenta('pwca_catalog'),magenta(datadir)))
for f in files:
    
    #
    print '.',

    #
    file_name = f.split('/')[-1].split('.')[0]

    #
    catalog = pwca_catalog

    #
    A = scsearch( keyword=file_name, verbose=not True, catalog=catalog )

    #
    a = A[0]

    #
    m1,m2 = a.m1,a.m2 
    eta = m1*m2/(m1+m2)

    #
    if 'th' in file_name:
        foo = file_name.split('th')[-1]
    elif 't' in file_name:
        foo = file_name.split('Dit')[0]
        foo = foo.split('t')[-1]
        foo = foo.lower().split('d')[0]
        foo = foo.split('_')[0]
    #
    theta = int(foo)
    
    #
    X1,X2,L = a.X1,a.X2,a.L
    
    #
    chi1 = norm(X1)
    chi2 = norm(X2)
    if abs(m1-m2)<1e-4:
        if (chi1<1e-4) and (chi2>1e-4):
            chi2 = 0
            chi1 = chi2
            X1,X2 = [ array(k) for k in (X2,X1) ]
            m1=m2=0.5
    
    #
    chi_p   = calc_chi_p(   m1,X1, m2,X2, L )
    chi_eff = calc_chi_eff( m1,X1, m2,X2, L )
    
    #
    delta = (m1-m2)/(m1+m2)
    
    #
    metadata.append( [ file_name,
                      theta,
                      eta,
                      delta,
                      chi_eff,
                      chi_p,
                      chi1,
                      chi2 ] )

#
print 'Done.'

#
metadata_array = array(metadata)
metadata_dict = {}

#
keys = [ 'simname', 'theta', 'eta', 'delta', 'chi_eff', 'chi_p', 'chi1', 'chi2' ]
metadata_dict = { keys[k]:metadata_array[k] for k in range(len(keys)) }

#
metadata_dict['array_data'] = metadata_array

#
metadata_path = '/Users/book/KOALA/puck/ll/data/metadata_dict'
alert('Saving metadata dictionary to %s'%magenta(metadata_path+'.npz'))
savez(metadata_path,metadata_dict)

#
alert('Saving complete.')
    
        

