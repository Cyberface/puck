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
package_dir = parent( pwca.__path__[0] )
datadir = package_dir + 'data/version2/'
files = glob( datadir+'q*.txt' )

#
fig,ax = subplots( len(files), 2, figsize=2.5*array([ 2.5*2/(0.618), 1.5*len(files) ]) )
ax = ax.flatten()

#
tight_layout(1,2,4)

#
alert(yellow('NOTE that we fit both phase and phase derivate to provide a validation option.'),header=True)

#
foo = {}

#
p = 0
amp_popt_array  = zeros( (len(files),2) )
phi_popt_array  = zeros( (len(files),4) )
dphi_popt_array = zeros( (len(files),4) )
amp_pcov_list, dphi_pcov_list, phi_pcov_list = [],[],[]
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
    f,amp_td,amp_fd,dphi_td,dphi_fd,phi_td,phi_fd = data.T
    theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2 = metadata_dict['array_data'][k]
    
    #
    physical_param_array[j,:] = metadata_dict['array_data'][k]
    
    # GENERATE TEMPLATE FUNCTIONS
    # ---
    template_amp, template_dphi, template_phi = pwca.template_amp_phase(m1, m2, chi1, chi2, chi_p)
    
    # DEFINE DATA TO USE FOR CALIBRATION
    # --- 
    # CALIBRATION_PHI  = phi_fd
    # CALIBRATION_DPHI = dphi_fd
    # CALIBRATION_AMP  = amp_fd
    CALIBRATION_PHI  = phi_td
    CALIBRATION_DPHI = dphi_td
    CALIBRATION_AMP  = amp_td
    

    # PHASE 
    # ---

    # NOTE that the td phase is used to exact consistency with the models of coprecessing angles
    phenomd_phi = template_phi(f)
    phi_popt, phi_pcov = curve_fit(template_phi, f, CALIBRATION_PHI,p0=[0,0,0,0])
    
    #
    phi_popt_array[j,:] = phi_popt
    phi_pcov_list.append( phi_pcov )
    
    # PHASE DERIVATIVE
    # ---

    # NOTE that the td phase derivative  is used to exact consistency with the models of coprecessing angles
    phenomd_dphi = template_dphi(f)
    dphi_popt, dphi_pcov = curve_fit(template_dphi, f, CALIBRATION_DPHI,p0=[0,0,0,0])
    best_fit_dphi = template_dphi(f,*dphi_popt)
    # Get the phase fit version of dphi
    best_fit__phi  = template_dphi(f,*phi_popt)
    
    #
    dphi_popt_array[j,:] = dphi_popt
    dphi_pcov_list.append( dphi_pcov )
    
    # AMPLITUDE
    # ---
    
    #
    amp_scale = f ** (7.0/6)
    inv_amp_scale = f ** (-7.0/6)
    
    #
    log_scaled_template_amp = lambda X,MU2,MU4: log(  template_amp(X,MU2,MU4)*amp_scale  )
    phenomd_amp = template_amp(f)
    
    # NOTE that the td amplitude is used to exact consistency with the models of coprecessing angles
    scaled_amp = CALIBRATION_AMP * amp_scale
    log_scaled_amp = log(scaled_amp)
    log_scaled_amp_popt, log_scaled_amp_pcov = curve_fit(log_scaled_template_amp, f, log_scaled_amp,p0=[0,0])
    best_fit_amp = exp(log_scaled_template_amp(f,*log_scaled_amp_popt)) * inv_amp_scale
    
    #
    amp_popt_array[j,:] = log_scaled_amp_popt
    amp_pcov_list.append( log_scaled_amp_pcov )
    
    # PLOTTING
    # ---
    
    #subplot(1,2,1)
    sca(ax[p]); p+=1
    plot( f, phenomd_dphi, label='PhenomD', ls='--',alpha=0.9,color='k',lw=2 )
    plot( f, dphi_td, label='NR:Precessing (TD)', color='k', alpha=0.15, lw=6 )
    plot( f, dphi_fd, label='NR:Precessing (FD)', color='k', alpha=0.30, lw=6, ls=':' )
    plot( f, best_fit_dphi, label='Best Fit (dphi fit)', color='r', ls='-',lw=2 )
    # plot( f, best_fit__phi, label='Best Fit (phi fit)', color='b', ls='-',lw=2 )
    title(simname,size=12,loc='left')
    legend(ncol=2,loc=1)
    ylabel(r'$\frac{d}{df}\arg(\tilde{h}_{22})$')
    xscale('log')
    #
    title(f_.split('/')[-1].split('.')[0],loc='left',size=12)
    if (p==(2*len(files)-1)) or (p==1):  xlabel('$fM$')
    
    #subplot(1,2,2)
    sca(ax[p]); p+=1
    plot( f, phenomd_amp, label='PhenomD', ls='--',alpha=0.9,color='k',lw=2 )
    plot( f, amp_td, label='NR:Precessing (TD)', color='k', alpha=0.15, lw=6 )
    plot( f, amp_fd, label='NR:Precessing (FD)', color='k', alpha=0.30, lw=6, ls=':' )
    plot( f, best_fit_amp, label='Best Fit (TD)', color='r', ls='-',lw=2 )
    title(simname,size=12,loc='left')
    xscale('log')
    yscale('log')
    legend(ncol=2,loc=3)
    #ylim(1e-1,2e1)
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
savetxt( data_path, amp_popt_array, header='see "template_amp_phase()" in ansatz.py; columns are mu1, mu2, mu4' )

# Phase derivative fit parameters
data_path = datadir+'fit_opt_dphase_parameters.txt'
alert('Saving %s to %s'%( magenta('dphi_popt_array'), magenta(data_path)) )
savetxt( data_path, dphi_popt_array, header='see "template_amp_phase())" in ansatz.py; columns are nu4, nu5, nu6' )

# Phase fit parameters
data_path = datadir+'fit_opt_phase_parameters.txt'
alert('Saving %s to %s'%( magenta('dphi_popt_array'), magenta(data_path)) )
savetxt( data_path, phi_popt_array, header='see "template_amp_phase()" in ansatz.py; columns are nu4, nu5, nu6' )

#
data_path = datadir+'fit_pcov_dphase.pickle'
alert('Saving dphi_pcov_list to %s'%magenta(data_path))
pickle.dump( dphi_pcov_list, open( data_path, "wb" ) )

#
data_path = datadir+'fit_pcov_phase.pickle'
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