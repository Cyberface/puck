
#
def template_amp_mrd( m1, m2, chi1, chi2, chip ):
    """
    amplitude merger-ringdown ansatz
    
    equation 19: (arxiv:1508.07253)
    
    londonl@mit.edu, pilondon2@gmail.com, 2020
    """
    
    #
    from numpy import exp,sqrt,pi,log
    import pwca, phenom
    
    #
    chi = (m1*chi1 + m2*chi2)/(m1+m2)
    eta = m1*m2 / ( m1+m2 )**2
    
    # Given eta and chi, calculate PhenomD fit params 

    # Amplitude parameters
    gamma1 = pwca.d.gamma1(eta,chi)
    gamma2 = pwca.d.gamma2(eta,chi)
    gamma3 = pwca.d.gamma3(eta,chi)

    # PhenomD remnant parameters
    final_spin = phenom.remnant.FinalSpin0815_s(eta,chi)
    fring = phenom.remnant.fring(eta,chi1,chi2,final_spin)
    fdamp = phenom.remnant.fdamp(eta,chi1,chi2,final_spin)
    
    # Leading order PN amplitude 
    amp0 = (sqrt(2.0/3.0)*sqrt(eta)) / pi**(1.0/6.0)
    
    #
    def template(f, mu1=0, mu2=0, mu3=0, mu4=0):
        '''
        f:      frequency (code units)
        mu1:    parameter modifying PhenomD gamma1
        mu2:    parameter modifying PhenomD gamma2
        mu3:    parameter modifying PhenomD gamma3
        mu4:    parameter modifying PhenomD ringdown central frequency
        '''

        # Define new paremeters
        new_gamma1 = gamma1 + ( chip * mu1 )
        new_gamma2 = gamma2 + ( chip * mu2 )
        new_gamma3 = gamma3
        new_fdamp = fdamp + chip*mu3
        new_fring = fring + chip*mu4
        new_dfring = f - new_fring

        #
        new_gamma3_fdamp = new_gamma3 * new_fdamp

        #
        part1 = new_gamma1 * new_gamma3_fdamp / ( new_dfring**2 + new_gamma3_fdamp**2 )
        part2 = exp(  - new_gamma2 * new_dfring / (new_gamma3_fdamp)  )

        #
        template_amplitude = part1 * part2
        
        #
        template_amplitude = template_amplitude*(sqrt(2.0/3.0)*sqrt(eta)) / pi**(1.0/6.0)
        # template_amplitude = template_amplitude*(sqrt(2.0/3.0)*sqrt(eta)) / pi**(1.0/6.0)*(f**(-7.0/6))

        #
        return template_amplitude

    #
    return template


#
def template_dphi_mrd( m1, m2, chi1, chi2, chip ):
    """
    phase merger-ringdown ansatz
    derivative of Ansatz for the merger-ringdown phase Equation 14 arXiv:1508.07253
    londonl@mit.edu, pilondon2@gmail.com, 2020
    """
    
    #
    from numpy import exp,sqrt,pi,ndarray
    import pwca, phenom
    
    #
    chi = (m1*chi1 + m2*chi2)/(m1+m2)
    eta = m1*m2 / ( m1+m2 )**2
    
    # Given eta and chi, calculate PhenomD fit params 

    # Phase parameters
    alpha1 = pwca.d.alpha1(eta,chi)
    alpha2 = pwca.d.alpha2(eta,chi)
    alpha3 = pwca.d.alpha3(eta,chi)
    alpha4 = pwca.d.alpha4(eta,chi)
    alpha5 = pwca.d.alpha5(eta,chi)

    # PhenomD remnant parameters
    final_spin = pwca.d.FinalSpin0815_s(eta,chi)
    fring = pwca.d.fring(eta,chi1,chi2,final_spin)
    fdamp = pwca.d.fdamp(eta,chi1,chi2,final_spin)
    
    #
    #def template(f, nu1=0, nu2=0, nu3=0, nu4=0, nu5=0, nu6=0):
    def template(f, nu4=0, nu5=0, nu6=0):
        '''
        f:      frequency (code units)
        nu1:    parameter modifying PhenomD alpha1
        nu2:    parameter modifying PhenomD alpha2
        nu3:    parameter modifying PhenomD alpha3
        nu4:    parameter modifying PhenomD alpha4
        nu5:    parameter modifying PhenomD alpha5
        nu6:    parameter modifying PhenomD ringdown fdamp
        '''
        
        #
        nu1 = 0
        nu2 = 0
        nu3 = 0
        #nu4 = 0
        #nu5 = 0
        #nu6 = 0

        # Define new paremeters -- Late inspiral, Plunge
        new_alpha1 = alpha1 + ( chip * nu1 )
        new_alpha2 = alpha2 + ( chip * nu2 )
        new_alpha3 = alpha3 + ( chip * nu3 )
        
        #Define new paremeters --  Merger
        new_alpha4 = alpha4 + ( chip * nu4 )
        new_fring = fring + chip*nu5
        new_fdamp = fdamp + chip*nu6
        new_dfring = f - alpha5*new_fring

        #
        part1 = new_alpha1 + new_alpha2*f**(-2.0) + new_alpha3*f**(-0.25) 
        part2 = new_alpha4 / (  new_fdamp*(1 + (new_dfring**2)/(new_fdamp**2))  )

        #
        template_dphi = part1 + part2
        template_dphi *= -1.0/eta
        
        #
        return template_dphi-min(template_dphi[(f>0.03)*f<(0.12)]) if isinstance(f,ndarray) else template_dphi

    #
    return template
 
 
#
def calc_effective_a1_theta_helper( eta,chi_eff,chi_p ):
    
    # Import usefuls
    from numpy import dot, arctan2, sqrt
    from positive import eta2m1m2
    
    # Convert to effective single-spin system
    # NOTE that the formulas below result from taking single-spin equations for chi_p(theta,a1) and chi_s(theta,a1), and inverting them to get theta and a1
    # ---
    # >>> Premilinaries
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
    l = L/norm(l)
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
def pwca_phi_mrd( f, m1, m2, chi1, chi2, chip, phi0=0,nu4=None,nu5=None,nu6=None ):
    """
    Generate phase of merger-ringown FD waveform.
    Adapted from equation 14 arXiv:1508.07253, and made to be consistent with template_dphi_mrd() in the pwca package.
    londonl@mit.edu/pilondon2@gmail.com 2020
    """
    
    # sqrootf = sqrt(Mf)
    # fpow1_5 = Mf * sqrootf
    # # // check if this is any faster: 2 sqrts instead of one pow(x,0.75)
    # fpow0_75 = sqrt(fpow1_5); # pow(f,0.75)

    # ans = -(model_pars['alpha2']/Mf) \
    #         + (4.0/3.0) * (model_pars['alpha3'] * fpow0_75) \
    #         + model_pars['alpha1'] * Mf \
    #         + model_pars['alpha4'] * rholm * arctan( (Mf - model_pars['alpha5'] * model_pars['fRD']) / (rholm * model_pars['fDM'] * taulm) )
    
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
    _,_,_,_,nu4,nu5,nu6 = pwca.generate_model_params(theta,eta,a1)

    # Phase parameters
    alpha1 = pwca.d.alpha1(eta,chi)
    alpha2 = pwca.d.alpha2(eta,chi)
    alpha3 = pwca.d.alpha3(eta,chi)
    alpha4 = pwca.d.alpha4(eta,chi)
    alpha5 = pwca.d.alpha5(eta,chi)

    # PhenomD remnant parameters
    final_spin = pwca.d.FinalSpin0815_s(eta,chi)
    fring = pwca.d.fring(eta,chi1,chi2,final_spin)
    fdamp = pwca.d.fdamp(eta,chi1,chi2,final_spin)
        
    # Start the precession parama party
    nu1 = nu2 = nu3 = 0
    if (nu4 is None) or (nu5 is None) or (nu6 is None):
        # Generate model parameters
        a1,theta = calc_effective_a1_theta_helper( eta,chi,chip )
        _,_,_,_,nu4,nu5,nu6 = generate_model_params(theta,eta,a1)
    
    # Define new paremeters -- Late inspiral, Plunge
    new_alpha1 = alpha1 + ( chip * nu1 )
    new_alpha2 = alpha2 + ( chip * nu2 )
    new_alpha3 = alpha3 + ( chip * nu3 )
    
    #Define new paremeters --  Merger
    new_alpha4 = alpha4 + ( chip * nu4 )
    new_fring = fring + chip*nu5
    new_fdamp = fdamp + chip*nu6
    new_dfring = f - alpha5*new_fring
    
    # NOTE that part1 is identical to the PhenomD equivalent
    f4by3 = 1.3333333333333333
    part1 = new_alpha1*f  -  new_alpha2/f  +  f4by3*new_alpha3*(f**0.75)
    # NOTE that we use arctan here not arctan2 because the denominator is always positive
    part2 = new_alpha4 * arctan( new_dfring / new_fdamp )
    
    # NOTE that the minus sign signals the phase convention used internally
    phi = phi0  -  (part1 + part2)/eta
    
    #
    return phi
    
#
def calc_fd_copreessing_merger_waveform( f, m1, m2, X1, X2, L ):
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
    from positive import calc_chi_p
    from pwca import template_amp_mrd,template_dphi_mrd
    
    # Compute physical parameters for model 
    l = L/norm(l)
    chi1,chi2 = [ dot(k,l) for k in (X1,X2) ]
    chi_p     = calc_chi_p(m1,X1,m2,X2,L)
    chi_eff   = calc_chi_eff(m1,X1,m2,X2,L)
    # Convert to effective single-spin system
    # NOTE that the function below result from taking single-spin equations for chi_p(theta,a1) and chi_s(theta,a1), and inverting them to get theta and a1
    # ---
    a1,theta = calc_effective_a1_theta( m1, m2, X1, X2, L )
    
    # Generate template amplitude function -- input chi_p along with PhenomD parameters
    template_amp  = template_amp_mrd(  m1, m2, chi1, chi2, chi_p )
    
    # Generate model parameters
    mu1,mu2,mu3,mu4,nu4,nu5,nu6 = generate_model_params(theta,eta,a1)
    
    # Evaluate phase model 
    model_phi     = pwca_phi_mrd( raw_fp, m1, m2, chi1, chi2, chip, a1, phi0=0,nu4=None,nu5=None,nu6=None )
    # Evaluate amplitude model
    scale_factor = raw_fp ** (-7.0/6.0)
    model_amp     = template_amp( raw_fp, mu1, mu2, mu3, mu4 ) * scale_factor
    
    #
    y = model_amp * exp( 1j * model_phi )
    
    #
    return y