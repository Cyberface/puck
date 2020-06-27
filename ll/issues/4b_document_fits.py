#!/usr/bin/env python2.7

# Setup python environment
from matplotlib.pyplot import *
from numpy import *
from positive import *
from nrutils import scsearch, gwylm
from glob import glob
import pwca
from pwca import determine_data_fitting_region,pwca_catalog,metadata_dict

# --------------------------------------- #
# Preliminaries
# --------------------------------------- #

#Load parameter space fit data
alert('Loading parameter space fit data.')
datadir = '/Users/book/KOALA/puck/ll/data/version2/'
foo_path = datadir+'parameter_space_fits.pickle'
foo = pickle.load( open( foo_path, "rb" ) )

# --------------------------------------- #
# Plot ans save fits 
# --------------------------------------- #

alert('Generate and save diagnostic plots ...')
for k in foo:
    if ('mu' in k)or('nu' in k):
        # Generate plot
        labels={'python':[k,('u', 'eta', 'chi_eff', 'chi_p'),''],'latex':[k,(r'\cos(\theta)', r'\eta', r'\chi_s', r'\chi_p'),'']}
        fig = foo[k].plot(labels=labels,size_scale=1.2)
        # Save figure 
        figure_path = datadir + k+'_fit_diagnostic.pdf'
        alert('Saving diagnostic plot to %s'%magenta(figure_path))
        savefig( figure_path, pad_inches=0 )
    
   
# --------------------------------------- #
# Generate fit python code 
# --------------------------------------- #

#
code_string = ['\n\n#\ndef generate_model_params(theta,eta,chi_eff,chi_p):\n\n',
               '\t\'\'\'\n\tHola, soy un codigo escribido por codigo. \n\t~londonl@mit.edu/pilondon2@gmail.com 2020\n\t\'\'\'  \n\n',
               '\t# Import usefuls\n',
               '\tfrom numpy import cos\n\n',
               '\t# Preliminaries\n',
               '\tu = cos(theta)\n','\tu2 = u*u\n', '\tu3 = u2*u\n', 
               '\tu3 = u2*u\n', '\teta2 = eta*eta\n', '\teta3 = eta2*eta\n', 
               '\tchi_eff2 = chi_eff*chi_eff\n', '\tchi_eff3 = chi_eff2*chi_eff\n', 
               '\tchi_p2 = chi_p*chi_p\n', 
               '\tchi_p3 = chi_p2*chi_p\n\n' 
              ]

# determine list of fitted variables and sort
fit_var = sort( [ k for k in foo.keys() if ('mu' in k)or('nu' in k) ] )

#
for k in fit_var:
    
    # Store python code for fit
    code_string.append( '\t# %s\n'%k )

    #
    this_code_string = foo[k].__str_python__()
    this_code_string = this_code_string.replace('lambda u,eta,chi_eff,chi_p: ','')
    this_code_string = this_code_string.replace('chi_p*chi_p','chi_p2')
    this_code_string = this_code_string.replace('chi_p2*chi_p','chi_p3')
    this_code_string = this_code_string.replace('u*u','u2')
    this_code_string = this_code_string.replace('u2*u','u3')
    this_code_string = this_code_string.replace('chi_eff*chi_eff','chi_eff2')
    this_code_string = this_code_string.replace('chi_eff2*chi_eff','chi_eff3')
    this_code_string = this_code_string.replace('eta*eta','eta2')
    this_code_string = this_code_string.replace('eta2*eta','eta3')

    #
    code_string.append( '\t'+this_code_string+'\n\n' )
        
#
code_string.append( '\t#\n' )
code_string.append( '\treturn %s\n'%(','.join(fit_var)) )

# Write fit equations to file 
codedir = '/Users/book/KOALA/puck/ll/pwca/'
code_path = codedir+'parameter_space_fits.py'
alert('Write fit equations to file at %s'%magenta(code_path))
f = open(code_path,'w+')
f.writelines(code_string)
f.close()

#
alert('All done.')
