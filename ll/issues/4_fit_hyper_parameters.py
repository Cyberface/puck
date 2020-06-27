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

#
alert('Loading parameter space data.')
datadir = '/Users/book/KOALA/puck/ll/data/version2/'

# Load and unpuack physical parameter space
raw_domain = loadtxt(datadir+'fit_intial_binary_parameters.txt')
theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2 = raw_domain.T


# Define desired model domain variables and array 
u = cos(theta)
model_domain = array( [ u, eta, chi_eff, chi_p ] ).T

# Load and unpuack physical parameter space -- dphi
dphi_range = loadtxt(datadir+'fit_opt_dphase_parameters.txt')
nu4,nu5,nu6 = dphi_range.T

# Load and unpuack physical parameter space -- amp
amp_range = loadtxt(datadir+'fit_opt_amplitude_parameters.txt')
mu1, mu2, mu3, mu4 = amp_range.T

#
foo = {}

# --------------------------------------- #
# Fit dphase parameters 
# --------------------------------------- #
'''NOTE that function setting have been determined by manual play'''
alert('Fitting dphase parameters ...')

# nu4
foo['nu4'] = gmvpfit( model_domain, nu4,fitatol=0.002,verbose=True,maxdeg_list=[3,2,1,1],center=True)
# nu5
foo['nu5'] = gmvpfit( model_domain, nu5,fitatol=0.002,verbose=True,maxdeg_list=[3,3,1,1],center=True)
# nu6
foo['nu6'] = gmvpfit( model_domain, nu6,fitatol=0.002,verbose=True,maxdeg_list=[3,3,1,1],center=False)

# --------------------------------------- #
# Fit amplitude parameters 
# --------------------------------------- #
'''NOTE that function setting have been determined by manual play'''
alert('Fitting amplitude parameters ...')

# mu1 
foo['mu1'] = gmvpfit( model_domain, mu1,fitatol=0.004,verbose=True,maxdeg_list=[2,3,1,2],center=False,temper=False)
# mu2
foo['mu2'] = gmvpfit( model_domain, mu2,fitatol=0.004,verbose=True,maxdeg_list=[2,2,2,2],center=True)
# mu3
foo['mu3'] = gmvpfit( model_domain, mu3,fitatol=0.004,verbose=True,maxdeg_list=[2,2,2,2],center=True)
# mu4
foo['mu4'] = gmvpfit( model_domain, mu4,fitatol=0.004,verbose=True,maxdeg_list=[2,2,2,2],center=False)

# --------------------------------------- #
# Saving fit data
# --------------------------------------- #
data_path = datadir+'parameter_space_fits.pickle'
alert('Saving '+yellow('parameter_space_fits')+' to %s'%magenta(data_path))
pickle.dump( foo, open( data_path, "wb" ) )

#
alert('All done!')