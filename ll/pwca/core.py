

#
import pickle
from positive import alert,magenta,parent
from numpy import loadtxt,load
from os.path import exists
import pwca
import json

#
package_dir = parent( pwca.__path__[0] )
data_dir = package_dir + 'data/'

# Always load catalog list for calibration runs 
pwca_catalog_path = data_dir+'pwca_catalog.pickle'
pwca_catalog = pickle.load( open( pwca_catalog_path, "rb" ) )
alert('Catalog of calibration runs stored to %s'%magenta('"pwca.pwca_catalog"'),fname='pwca.core')

# Always load curated metadata for calibration runs 
metadata_dict_path = data_dir+'metadata_dict.pickle'
metadata_dict = load(metadata_dict_path,allow_pickle=True)
alert('Metadata dictionary for calibration runs stored to %s'%magenta('"pwca.metadata_dict"'),fname='pwca.core')

#
catalog_paper_md_path = data_dir+'catalog_paper_metadata.json'
with open(catalog_paper_md_path, 'r') as f:
    catalog_paper_metadata = json.load(f)
alert('Metadata dictionary for Ed\'s catalog paper stored to %s'%magenta('"pwca.catalog_paper_metadata"'),fname='pwca.core')

# Function to determine version2 data fitting region
def determine_data_fitting_region( data, fmin=0.03, fmax=0.12 ):
    '''
    Given version2 data array, determine fitting region dynamically.
    This function assumes clean data within the default or input values of fmin and fmax, and then uses the phase derivate to determine new fmin and fmax values that are ultimately used to define a fitting region.
    '''
    
    # Import usefuls
    from numpy import argmin,log
    from positive import smooth,find,lim
    from matplotlib.pyplot import figure,plot,show,axhline,xlim,ylim
    
    # Extract data 
    f,amp_td,amp_fd,dphi_td,dphi_fd,phi_td,phi_fd = data

    # Use default domain bounds to determine a mask
    mask = (f>=fmin) & (f<=fmax)
    
    # Determine the minimum dphi
    # Smooth dphi using postiive's savgol filter
    x = log(f[mask])
    y = smooth(dphi_td[mask]).answer
    knot = argmin(y)
    y_knot = y[knot]
    data[3] = dphi_td - smooth(dphi_td[mask]).answer[knot] + y_knot
    data[4] = dphi_fd - smooth(dphi_fd[mask]).answer[knot] + y_knot
    
    # Determine new fmin and max using heuristic 
    f_knot = f[mask][knot]
    new_fmin = max(f_knot * 0.22,0.018) # 0.5 # 0.325
    new_fmax = f_knot + 0.020 # 0.025 
    
    #
    new_mask = (f>=new_fmin) & (f<=new_fmax)
    new_data = data.T[new_mask,:]
    
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
    return theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2 
    
# Given underlying physical parameters, calculate ones useful form modeling
def parama_party( eta,theta,a1 ):
    '''
    PARAMA-PARTY:
    If L || z and m1>m2 and q=m1/m2, then 

    S2 = 0
    S1 = m1**2 a1 * exp( 1j * theta ) = Sz + 1j*Sperp
    X1 = X = S1/m1**2

    chi_eff = m1*a1*cos(theta)/(m1+m2) = a1*cos(theta)*/(1+1.0/q)

    A1 = 2 + (3*m2)/(2*m1)
    A2 = 2 + (3*m1)/(2*m2)
    B1 = A1 * a1*sin(theta)
    B2 = 0
    chi_p = max( B1,B2 ) / ( A1 * m1*m1 )
    L = L

    '''
    
    #
    from positive import eta2m1m2
    from numpy import cos,sin,maximum
    
    #
    m1,m2 = eta2m1m2(eta)
    
    #
    q = m1/m2
    chi_eff = m1*a1*cos(theta)/(m1+m2)
    
    #
    A1 = 2 + (3.0*m2)/(2.0*m1)
    Norm_S1_perp = abs(a1*sin(theta)*m1*m1)
    B1 = A1 * Norm_S1_perp 
    chi_p = maximum( B1,0 ) / ( A1 * m1*m1 )
    
    '''
    NOTE that the above is basically a1*sin(theta), but we retain the additional steps to illustrate consistency with the more general formula
    '''
    
    #
    return chi_eff, chi_p
     
# Advanced gloss atop mvpolyfit.plot and mvrfit.plot
def advanced_gmvx_plot( fit_object ):
    
    '''
    Advanced gloss atop mvpolyfit.plot and mvrfit.plot
    '''
    
    from matplotlib.pyplot import subplots, plot, xlabel, ylabel, title, sca, gca, figaspect, tight_layout
    from numpy import cos,sin,array,around,ones_like,sort,pi,linspace
    from positive import eta2q,q2eta,eta2delta
    from glob import glob
    from pwca import determine_data_fitting_region,pwca_catalog,metadata_dict
    
    # Load and unpuack physical parameter space
    raw_domain = loadtxt(data_dir+'version2/fit_intial_binary_parameters.txt')
    theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2 = raw_domain.T


    # Define desired model domain variables and array 
    u = cos(theta)
    v = sin(theta)
    q = 1.0/eta2q(eta)

    # --------------------------------------- #
    # Plot ans save fits 
    # --------------------------------------- #

    # Collect set of unique a1 values
    a1_point = around(a1,2)
    a1_set = array(sort(list( set(a1_point) )))

    # Collect set of unique angle values
    degree_point = (theta*180/pi).astype(int)
    theta_point = degree_point*pi/180
    theta_set = array(sort(list( set(theta_point) )))

    # Collect set of unique mass-ratio values
    q_point = around(array([eta2q(n) for n in eta]),2)
    q_set = array(sort(list( set(q_point) )))

    # Collect set of unique eta values
    eta_point = q2eta( q_point )
    eta_set = q2eta(q_set)

    # Summary figure for internal diagnostics 
    summary_fig = fit_object.plot(size_scale=1.5)
    ax = summary_fig.axes
    
    #
    sca(ax[0])
    title(fit_object.labels['python'][0])
    
    #
    num_figs = len(a1_set)*len(theta_set)
    eta_set_figs,set_fig_ax = subplots( len(a1_set), len(theta_set), figsize=5*array([ len(theta_set),len(a1_set) ]) )
    set_fig_ax = set_fig_ax.flatten();
    tight_layout(4,4)
    ax_counter = 0

    #
    for _a1 in a1_set:
        for _theta in theta_set:

            #
            theta_mask = (_theta==theta_point)
            a1_mask = (_a1==a1_point)
            mask = a1_mask & theta_mask

            #
            _eta = eta_point[mask]
            _u = cos(_theta) 

            #
            case_eta   = linspace( min(_eta),max(_eta),1000 ) 
            case_delta = eta2delta( case_eta )
            case_q     = 1.0/eta2q(case_eta)  
            case_theta = _theta * ones_like(case_eta)
            case_u     = cos(case_theta)
            case_a1    = _a1    * ones_like(case_eta)

            #
            case_domain = array([case_u,case_eta,case_delta,case_a1]).T
            case_range = fit_object.eval(case_domain)
            opt_range  = fit_object.eval(fit_object.domain[mask,:])

            #
            sca(ax[0])
            ax[0].plot3D( case_u, case_eta, case_range, lw=1, alpha=1, color = 'tab:blue' if _a1==a1_set[0] else 'red' )
            
            #
            sca( set_fig_ax[ax_counter] ); ax_counter += 1
            plot( eta[mask], fit_object.range[mask] if hasattr(fit_object,'range') else fit_object.scalar_range[mask], marker='o',ls='none'  )
            plot( eta[mask], opt_range, marker='o',ms=10,mfc='none', color='r',ls='none'  )
            plot( case_eta, case_range, ls='-', color='r' )
            title( r'$a_1=%1.2f$, $\theta=%1.2f$'%(_a1,round(_theta*180.0/pi,0)) )
            xlabel(r'$\eta$')
            ylabel(r'$\%s$'%fit_object.labels['python'][0])
            

    
    #
    num_figs = len(a1_set)*len(eta_set)
    theta_set_figs,set_fig_ax = subplots( len(a1_set), len(eta_set), figsize=5*array([ len(eta_set),len(a1_set) ]) )
    set_fig_ax = set_fig_ax.flatten();
    tight_layout(4,4)
    ax_counter = 0

    #
    for _a1 in a1_set:
        for _eta in eta_set:

            #
            eta_mask = (_eta==eta_point)
            a1_mask = (_a1==a1_point)
            mask = a1_mask & eta_mask

            #
            _theta = theta_point[mask]
            _u = cos(_theta) 

            #
            case_theta   = linspace( min(_theta),max(_theta),1000 ) # 
            case_u     = cos(case_theta)
            case_eta   = _eta * ones_like(case_theta)
            case_delta = eta2delta( case_eta )
            case_a1    = _a1  * ones_like(case_theta)

            #
            case_domain = array([case_u,case_eta,case_delta,case_a1]).T
            case_range = fit_object.eval(case_domain)
            opt_range  = fit_object.eval(fit_object.domain[mask,:])

            #
            sca(ax[0])
            ax[0].plot3D( case_u, case_eta, case_range, lw=1, alpha=1, color = 'tab:blue' if _a1==a1_set[0] else 'red' )
            
            #
            sca( set_fig_ax[ax_counter] ); ax_counter += 1
            plot( cos(theta[mask]), fit_object.range[mask] if hasattr(fit_object,'range') else fit_object.scalar_range[mask], marker='o',ls='none',color='r'  )
            plot( cos(theta[mask]), opt_range, marker='o',ms=10,mfc='none', color='b',ls='none'  )
            plot( cos(case_theta), case_range, ls='-', color='b' )
            title( r'$a_1=%1.2f$, $q=%1.2f$'%(_a1,eta2q(_eta)) )
            xlabel(r'$\cos(\theta)$')
            ylabel(r'$\%s$'%fit_object.labels['python'][0])
            
    #
    return summary_fig, eta_set_figs, theta_set_figs
 
#
def calc_effective_a1_theta_helper( eta,chi_eff,chi_p ):
    
    # Import usefuls
    from numpy import dot, arctan2, sqrt
    from positive import eta2m1m2
    
    # Convert to effective single-spin system
    # NOTE that the formulas below result from taking single-spin equations for chi_p(theta,a1) and chi_s(theta,a1), and inverting them to get theta and a1
    # ---
    # >>> Preliminaries
    m1,m2 = eta2m1m2(eta) # NOTE m1>m2 convention
    if m1<m2: error('m1<m2, check positive\'s eta2m1m2 as it should output m1>=m2')
    M = m1+m2
    M2 = M*M
    chi_eff2 = chi_eff*chi_eff
    chi_p2 = chi_p*chi_p
    m1_2 = m1*m1
    effective_a1_times_m1 = sqrt( chi_p2*m1_2 + chi_eff2*M2 )
    # >>> Compute effective spin1 (dimensionless) magnitude
    a1 = effective_a1_times_m1 / m1
    # >>> and lastly compute effective angle between L and S
    # NOTE (tangentially) that mathematica uses a different order for its two-input ArcTan function
    theta = arctan2(  chi_p*m1/effective_a1_times_m1, chi_eff*M/effective_a1_times_m1  )
    
    #
    return a1,theta

#
def calc_effective_a1_theta( m1, m2, X1, X2, L ):
    '''
    Given general binary parameters (quasi-cicular etc), compute effective single-spin parameters a1=|X1_single_spin| and theta="angle between L and X1+X2"
    '''
    
    # Import usefuls
    from positive import calc_chi_p, calc_chi_eff
    from numpy import array, dot, cos, arccos, arctan2
    from numpy.linalg import norm
    
    # Preliminaries
    M = m1+m2 
    m1,m2 = [ k/M for k in (m1,m2) ]
    
    # Compute relevant physical parameters 
    l = L/norm(L)
    chi1,chi2 = [ dot(k,l) for k in (X1,X2) ]
    chi_p     = calc_chi_p(m1,X1,m2,X2,L)
    chi_eff   = calc_chi_eff(m1,X1,m2,X2,L)
    eta = m1*m2/(m1+m2)**2
    
    # Convert to effective single-spin system
    # NOTE that the formulas below result from taking single-spin equations for chi_p(theta,a1) and chi_s(theta,a1), and inverting them to get theta and a1
    # ---
    a1,theta = calc_effective_a1_theta_helper( eta,chi_eff,chi_p )
    
    #
    return a1,theta 

#
def pwca_phi_helper(f, theta, eta, a1, chi1, chi2, chip, nu4, nu5, nu6, zeta2, phi0=0, Mtotal=None, fref=None, phiRef=None):
    
    """
    Same as pwca_phi but with input being actual parameters used for model.
    Generate phase of merger-ringown FD waveform.
    Adapted from equation 14 arXiv:1508.07253, and made to be consistent with template_dphi_mrd() in the pwca package.
    londonl@mit.edu/pilondon2@gmail.com 2020
    """
    
    # Import usefuls
    from numpy import exp,sqrt,pi,ndarray,arctan,cos
    import pwca, phenom
    from positive import eta2m1m2
    
    # Generate model parameters
    if (nu4 is None) or (nu5 is None) or (nu6 is None) or (zeta2 is None):
        _,_,_,_,_,nu4,nu5,nu6,zeta2 = pwca.generate_model_params(theta,eta,a1)
    
    # NOTE that the minus sign signals the phase convention used internally
    m1,m2 = eta2m1m2(eta)
    _,_,template_phi = pwca.template_amp_phase(m1, m2, chi1, chi2, chip, Mtotal=Mtotal, fref=fref, phiRef=phiRef)
#     phi = phi0  -  template_phi( f, nu4, nu5, nu6, zeta2 )
    phi = template_phi( f, nu4, nu5, nu6, zeta2)
    
    return phi


#
def pwca_phi( f, m1, m2, chi1, chi2, chip, nu4, nu5, nu6, zeta2, phi0=0,Mtotal=None, fref=None, phiRef=None ):
    """
    Generate phase of merger-ringown FD waveform.
    Adapted from equation 14 arXiv:1508.07253, and made to be consistent with template_dphi_mrd() in the pwca package.
    londonl@mit.edu/pilondon2@gmail.com 2020
    """
    
    # Import usefuls
    from numpy import exp,sqrt,pi,ndarray,arctan
    import pwca, phenom
    
    #
    chi = (m1*chi1 + m2*chi2)/(m1+m2)
    eta = m1*m2 / ( m1+m2 )**2
    
    # Convert to effective single-spin system
    # NOTE that the function below result from taking single-spin equations for chi_p(theta,a1) and chi_s(theta,a1), and inverting them to get theta and a1
    # ---
    a1,theta = calc_effective_a1_theta_helper( eta,chi,chip )
    
    # Generate model parameters
    if (nu4 is None) or (nu5 is None) or (nu6 is None) or (zeta2 is None):
        _,_,_,_,_,nu4,nu5,nu6,zeta2 = pwca.generate_model_params(theta,eta,a1)
    
    # NOTE that the minus sign signals the phase convention used internally
    _,_,template_phi = pwca.template_amp_phase(m1, m2, chi1, chi2, chip, Mtotal=Mtotal, fref=fref, phiRef=phiRef)
#     phi = phi0  -  template_phi( f, nu4, nu5, nu6, zeta2 )
    phi = template_phi( f, nu4, nu5, nu6, zeta2)
    
    #
    return phi
       
#
def generate_pwca_waveform_helper( f, theta, eta, a1, chi1, chi2, chi_p, Mtotal=None, fref=None, phiRef=None ):
    '''
    DESCRIPTION
    ---
    Same as generate_pwca_waveform but with input being actual parameters used for model.
    
    USAGE
    ---
    fd_waveform_array = pwca_waveform( f, m1, m2, X1, X2, L )
    
    f,       frequency array of desired values; frequencies are geometric: M*f_hz with M=1
    theta,
    eta,
    a1,
    chip,       the total in-plane spin or chi_p (either may be used with care in the correct context)
    '''
    
    # Import usefuls
    from numpy import dot,exp,cos
    from positive import eta2m1m2
    from numpy.linalg import norm
    from positive import calc_chi_p,calc_chi_eff
    from pwca import template_amp_phase, pwca_phi, calc_effective_a1_theta, generate_model_params
    
    # Generate template amplitude function -- input chi_p along with PhenomD parameters
    m1,m2 = eta2m1m2(eta)
    template_amp,_,_  = template_amp_phase(  m1, m2, chi1, chi2, chi_p, Mtotal=Mtotal, fref=fref, phiRef=phiRef )
    
    # Generate model parameters
    mu2,mu4,nu4,nu5,nu6,zeta2 = generate_model_params(theta,eta,a1)
    
    # Evaluate phase model 
    model_phi     = pwca_phi_helper( f, theta, eta, a1, chi1, chi2, chi_p, nu4, nu5, nu6, zeta2, Mtotal=Mtotal, fref=fref, phiRef=phiRef )
    # Evaluate amplitude model
    scale_factor = 1
    model_amp     = template_amp( f, mu2, mu4 ) * scale_factor
    
    # Compute complex waveform
    # NOTE minus sign added to be consistent with external phase conventions
    y = model_amp * exp( -1j * model_phi )
    
    #
    return y

#
def generate_pwca_waveform( f, m1, m2, X1, X2, L, Mtotal=None, fref=None, phiRef=None ):
    '''
    DESCRIPTION
    ---
    Calculate waveform for BBH coprecessing merger-ringdown given inputs described below
    
    USAGE
    ---
    fd_waveform_array = pwca_waveform( f, m1, m2, X1, X2, L )
    
    f,       frequency array of desired values; frequencies are geometric: M*f_hz with M=1
    m1,      mass of larger component, m1+m2=1
    m2,      mass of smaller component, m1+m2=1
    X1,      initial spin vector for larger component
    X2,      initial spin vector for smaller component
    L,       initial orbital angular momentum direction
    '''
    
    # Import usefuls
    from numpy import dot,exp
    from numpy.linalg import norm
    from positive import calc_chi_p,calc_chi_eff
    from pwca import template_amp_phase, pwca_phi, calc_effective_a1_theta, generate_model_params
    
    # Compute physical parameters for model 
    l = L/norm(L)
    chi1,chi2 = [ dot(k,l) for k in (X1,X2) ]
    chi_p     = calc_chi_p(m1,X1,m2,X2,L)
    chi_eff   = calc_chi_eff(m1,X1,m2,X2,L)
    # Convert to effective single-spin system
    # NOTE that the function below result from taking single-spin equations for chi_p(theta,a1) and chi_s(theta,a1), and inverting them to get theta and a1
    # ---
    a1,theta = calc_effective_a1_theta( m1, m2, X1, X2, L )
    
    # Generate template amplitude function -- input chi_p along with PhenomD parameters
    template_amp,_,_  = template_amp_phase(  m1, m2, chi1, chi2, chi_p, Mtotal=Mtotal, fref=fref, phiRef=phiRef )
    
    # Generate model parameters
    eta = m1*m2/((m1+m2)**2)
    mu2,mu4,nu4,nu5,nu6,zeta2 = generate_model_params(theta,eta,a1)
    
    # Evaluate phase model 
    model_phi     = pwca_phi( f, m1, m2, chi1, chi2, chi_p, nu4, nu5, nu6, zeta2, Mtotal=Mtotal, fref=fref, phiRef=phiRef )
    
    # Evaluate amplitude model
    scale_factor = 1
    model_amp     = template_amp( f, mu2, mu4 ) * scale_factor
    
    # Compute complex waveform
    # NOTE minus sign added to be consistent with external phase conventions
    y = model_amp * exp( -1j * model_phi )
    
    #
    return y
    
#
def __generate_modified_phenomd_helper__(f, m1, m2, chi1,chi2,chi_p, fref=None, phiRef=None):
    
    # Import usefuls
    import pwca, phenom
    from numpy import dot,exp,array
    from numpy.linalg import norm
    from phenom.utils.utils import UsefulPowers
    from positive import calc_chi_p,calc_chi_eff
    
    # # Compute physical parameters for model 
    # l = L/norm(L)
    # chi1,chi2 = [ dot(k,l) for k in (X1,X2) ]
    # chi_p     = calc_chi_p(m1,X1,m2,X2,L)
    
    # Create PhenomD instance (PhenomD Object -- PDO) with appropriate final spin function input 
    pdo = phenom.PhenomD(m1=m1, m2=m2, chi1z=chi1, chi2z=chi2,chip=chi_p,finspin_func = "FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH",fRef=fref,phiRef=phiRef)
    
    #
    pdo_amp  = array( [pdo.IMRPhenomDAmplitude(k, pdo.model_pars, UsefulPowers(k)) for k in f] )
    
    #
    pdo_phi  = array( [pdo.IMRPhenomDPhase(k,pdo.p['eta'], pdo.model_pars, UsefulPowers(k)) - ((pdo.model_pars['t0'])*(k-pdo.model_pars['MfRef']) + pdo.model_pars['phi_precalc']) for k in f] )
    
    #
    y = pdo_amp * exp( -1j * pdo_phi)
    
    #
    return y
    
#
def generate_modified_phenomd(f, m1, m2, X1, X2, L, fref=None, phiRef=None):
    
    # Import usefuls
    import pwca, phenom
    from numpy import dot,exp,array
    from numpy.linalg import norm
    from phenom.utils.utils import UsefulPowers
    from positive import calc_chi_p,calc_chi_eff
    
    # Compute physical parameters for model 
    l = L/norm(L)
    chi1,chi2 = [ dot(k,l) for k in (X1,X2) ]
    chi_p     = calc_chi_p(m1,X1,m2,X2,L)
    
    #
    y = __generate_modified_phenomd_helper__( f, m1, m2, chi1, chi2, chi_p, fref=fref, phiRef=phiRef )
    
    #
    return y