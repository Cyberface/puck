#!/usr/bin/env python2.7

# Setup ipython environment
from matplotlib.pyplot import *
from numpy import *
from positive import *
from nrutils import scsearch, gwylm
from glob import glob
import pwca
from pwca import determine_data_fitting_region,pwca_catalog,metadata_dict


#
from numpy.linalg import norm
from scipy.optimize import curve_fit

#
datadir = '/Users/book/KOALA/puck/ll/data/version2/'
files = glob( datadir+'q*.txt' )

#
fig,ax = subplots( len(files), 2, figsize=2.5*array([ 2.5*2/(0.618), 1.5*len(files) ]) )
ax = ax.flatten()

#
tight_layout(1,2,4)

#
foo = {}

#
p = 0
amp_popt_array = zeros( (len(files), 4) )
dphi_popt_array = zeros( (len(files),3) )
amp_pcov_list = []
dphi_pcov_list = []
physical_param_array = zeros( (len(files), 11) )
alert('Plotting ...')
for j,f_ in enumerate(files):

    #
    simname = f_.split('/')[-1].split('.')[0]
    
    # Find index location of metadata for simname 
    k = list( metadata_dict['simname'] ).index(simname)
    
    # Load data for this case
    raw_data = loadtxt(f_).T
    data,_,fmin,fmax,fknot = determine_data_fitting_region(raw_data)
    
    # Collect params for this case 
    metadata = metadata_dict['array_data'][k,:]
    
    #
    f,amp_td,amp_fd,dphi_td,dphi_fd = data.T
    theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2 = metadata_dict['array_data'][k]
    
    #
    physical_param_array[j,:] = metadata_dict['array_data'][k]
    
    # PHASE DERIVATIVE
    # ---
    
    #
    template_phi = pwca.ansatz.template_dphi_mrd(m1, m2, chi1, chi2, chi_p)
    phenomd_dphi = template_phi(f)
    
    #
    dphi_popt, dphi_pcov = curve_fit(template_phi, f, dphi_td,p0=[0,0,0])
    best_fit_dphi = template_phi(f,*dphi_popt)
    
    #
    dphi_popt_array[j,:] = dphi_popt
    dphi_pcov_list.append( dphi_pcov )
    
    # AMPLITUDE
    # ---
    
    #
    amp_scale = f ** (7.0/6)
    inv_amp_scale = f ** (-7.0/6)
    
    #
    template_amp = pwca.ansatz.template_amp_mrd(m1, m2, chi1, chi2, chi_p)
    phenomd_amp = template_amp(f) * inv_amp_scale
    
    #
    log_scaled_amp_fd = log( amp_fd * amp_scale )
    log_template_amp = lambda F,*args: log( template_amp(F,*args) )
    amp_popt, amp_pcov = curve_fit(log_template_amp, f, log_scaled_amp_fd,p0=[0,0,0,0])
    best_fit_amp = template_amp(f,*amp_popt) * inv_amp_scale
    
    #
    amp_popt_array[j,:] = amp_popt
    amp_pcov_list.append( amp_pcov )
    
    # PLOTTING
    # ---
    
    #subplot(1,2,1)
    sca(ax[p]); p+=1
    plot( f, phenomd_dphi, label='PhenomD', ls='--',alpha=0.9,color='k',lw=2 )
    plot( f, dphi_td, label='NR:Precessing', color='k', alpha=0.15, lw=6 )
    plot( f, best_fit_dphi, label='Best Fit', color='r', ls='-',lw=2 )
    title(simname,size=12,loc='left')
    legend(ncol=3,loc=1)
    ylabel(r'$\frac{d}{df}\arg(\tilde{h}_{22})$')
    #
    title(f_.split('/')[-1].split('.')[0],loc='left',size=12)
    if (p==(2*len(files)-1)) or (p==1):  xlabel('$fM$')
    
    #subplot(1,2,2)
    sca(ax[p]); p+=1
    plot( f, phenomd_amp, label='PhenomD', ls='--',alpha=0.9,color='k',lw=2 )
    plot( f, amp_fd, label='NR:Precessing', color='k', alpha=0.15, lw=6 )
    plot( f, best_fit_amp, label='Best Fit', color='r', ls='-',lw=2 )
    title(simname,size=12,loc='left')
    yscale('log')
    legend(ncol=3,loc=3)
    ylabel(r'$|\tilde{h}_{22}(f)|$')
    #
    title(f_.split('/')[-1].split('.')[0],loc='left',size=12)
    if (p==(2*len(files))) or (p==2):  xlabel('$fM$')
    
    #
    print '.',
        
#
print ''
alert('Done.')
file_path = datadir+'waveform_fit_diagnostic.pdf'
alert('Saving batch plot to %s'%magenta(file_path))
savefig(file_path,pad_inches=2, bbox_inches = "tight")

# SAVE FIT DATA
# --

# Initial binary parameters
data_path = datadir+'fit_intial_binary_parameters.txt'
alert('Saving %s to %s'%( magenta('physical_param_array'), magenta(data_path)) )
savetxt( data_path, physical_param_array, header='see "issues/3_collect_metadata.py"; columns are theta, m1, m2, eta, delta, chi_eff, chi_p, chi1, chi2, a1, a2' )

# Amplitude fit parameters
data_path = datadir+'fit_opt_amplitude_parameters.txt'
alert('Saving %s to %s'%( magenta('amp_popt_array'), magenta(data_path)) )
savetxt( data_path, amp_popt_array, header='see "template_amp_mrd()" in ansatz.py; columns are mu1, mu2, mu4' )

# Phase derivative fit parameters
data_path = datadir+'fit_opt_dphase_parameters.txt'
alert('Saving %s to %s'%( magenta('dphi_popt_array'), magenta(data_path)) )
savetxt( data_path, dphi_popt_array, header='see "template_dphi_mrd()" in ansatz.py; columns are nu4, nu5, nu6' )

#
data_path = datadir+'fit_pcov_dphase.pickle'
alert('Saving dphi_pcov_list to %s'%magenta(data_path))
pickle.dump( dphi_pcov_list, open( data_path, "wb" ) )

#
data_path = datadir+'fit_pcov_amp.pickle'
alert('Saving amp_pcov_list to %s'%magenta(data_path))
pickle.dump( amp_pcov_list, open( data_path, "wb" ) )

#
alert('Fitting complete.')
alert('Plotting complete.')
alert('Saving complete.')