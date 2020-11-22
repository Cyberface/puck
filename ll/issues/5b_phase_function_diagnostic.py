#!/usr/bin/env python2.7

# Setup python environment
from matplotlib.pyplot import *
from numpy import *
from positive import *
from nrutils import scsearch, gwylm
from pwca import *
from glob import glob
import pwca

# Display basic info about calibration runs
scsearch( catalog=pwca_catalog, verbose=True )


# Preliminaries 
# ---

# Define data location
package_dir = parent( pwca.__path__[0] )
datadir = package_dir + 'data/version2/'

# Load select waveforms
files = glob( datadir+'q*.txt' )
simnames = [ f_.split('/')[-1].split('.t')[0] for f_ in files ]

# Setup plot
# ---
fig,ax = subplots( len(simnames), 1, figsize = 16*array([ 1, len(simnames)*0.618/2 ]) )
ax = ax.flatten()

#
tight_layout(1,2,1)

# Plot diagnostics for all calibration cases
ax_id = 0
for a in pwca_catalog:

    # Find index corresponding to chosen case
    k = simnames.index( a.simname )
    # Load data for this case
    raw_data = loadtxt(files[k]).T
    data,_,fmin,fmax,fknot = determine_data_fitting_region(raw_data)

    # Load and unpuack OPTIMAL physical parameter space -- dphi
    dphi_range = loadtxt(datadir+'fit_opt_dphase_parameters.txt')
    opt_nu2,opt_nu4,opt_nu5,opt_nu6 = dphi_range[k,:]
    # Load and unpuack OPTIMAL physical parameter space -- amp
    amp_range = loadtxt(datadir+'fit_opt_amplitude_parameters.txt')
    # opt_mu0, opt_mu1, opt_mu2, opt_mu3, opt_mu4 = amp_range[k,:]

    # extract useful params from scentry object
    theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2 = select_scenty_metadata(a)
    # generate model parameters 
    mu0,mu1,mu2,mu3,mu4,nu4,nu5,nu6,zeta2 = generate_model_params(theta,eta,a1)

    # generate template functions
    _, template_dphi, _ = pwca.template_amp_phase(m1, m2, chi1, chi2, chi_p)
    
    # Load and unpuack physical parameter space
    raw_domain = loadtxt(datadir+'fit_intial_binary_parameters.txt')
    param_test_quantity = sum(raw_domain[k,:]-array([theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2]))==0
    if param_test_quantity:
        alert(bold(green('CHECK PASSED: '))+'Generated physical parameters are identical to calibration ones.')
    else:
        error(bold(red('CHECK PASSED: '))+'Generated physical parameters are not identical to calibration ones.')

    #
    raw_f,raw_amp_td,raw_amp_fd,raw_dphi_td,raw_dphi_fd,raw_phi_td,raw_phi_fd = raw_data
    adjusted_raw_dphi_td = raw_dphi_td-min( raw_dphi_td[ (raw_f>=fmin) & (raw_f<=fmax)  ])
    f,amp_td,amp_fd,dphi_td,dphi_fd,phi_td,phi_fd = data.T
    #
    raw_positive_mask = raw_f>0
    raw_fp = raw_f[raw_positive_mask]
    
    # ------------------------------------ # 
    # GENERATE END-MODEL PHASE
    model_phi = pwca_phi( raw_fp, m1, m2, chi1, chi2, chi_p, nu4, nu5, nu6, zeta2 )
    # ------------------------------------ # 
    
    # 
    d_model_phi = -spline_diff( raw_fp, model_phi )
    adjusted_d_model_phi = d_model_phi - min(d_model_phi[(raw_fp>=fmin) & (raw_fp<=fmax)])

    #
    phenomd_dphi   = template_dphi( raw_fp )
    model_dphi     = template_dphi( raw_fp, nu4, nu5, nu6, zeta2 )
    
    #
    sca(ax[ax_id]); ax_id += 1
    
    #
    plot( raw_f, adjusted_raw_dphi_td, label='Symmetrised NR', color='k', alpha=0.5, lw=2 )
    #
    plot( raw_fp, model_dphi, ls='-', color='c', lw=2, label='End Model' )
    plot( raw_fp, phenomd_dphi, color='k', ls=':', alpha=1, label='PhenomD', lw=2 )
    
    #
    plot( raw_fp, adjusted_d_model_phi, color='r', ls='--', label='*Derivative of End Model Phase', lw=2 )
    
    #
    ylim(lim(dphi_td,dilate=0.1)+array([0,300]))
    xlim(lim(f,dilate=0.5))

    #
    axvspan( min(xlim()), fmin, alpha=0.05, color='k')
    axvspan( fmax, max(xlim()), alpha=0.05, color='k')
    
    #
    legend(ncol=2,loc=1)
    ylabel(r'$\frac{d}{df}\arg(\tilde{h}_{22})$')
    if (ax_id==len(pwca_catalog)) or (ax_id==1):
        xlabel('$fM$')
    title(a.simname,loc='left',size=12)

#
file_path = datadir+'phase_function_diagnostics.pdf'
alert('Saving batch plot to %s'%magenta(file_path))
savefig(file_path,pad_inches=2, bbox_inches = "tight")
alert('Done.')