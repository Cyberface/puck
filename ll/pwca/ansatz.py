
#
def template_amp_phase(m1, m2, chi1, chi2, chip):
    '''
    Given system masses and L-parallel dimensionless spins, output amplitude and phase methods. 
    
    USAGE
    ---
    template_amp, template_dphi, template_phi = template_amp_phase(m1, m2, chi1, chi2, chip)
    
    '''
    
    #
    from phenom.utils.utils import UsefulPowers
    from numpy import array
    import pwca, phenom
    
    # Create PrototypePhenomDCoprecess instance (PhenomD Object -- PDO)
    pdo = phenom.PrototypePhenomDCoprecess(m1=m1, m2=m2, chi1z=chi1, chi2z=chi2)
    
    #
    def template_amp( f, mu0=0, mu1=0, mu2=0, mu3=0, mu4=0 ):
        
        # Recompute internal parameters -- including connection coefficients for continuous inspiral to merger-ringdown
        pdo.model_pars = pdo.compute_model_parameters( pdo.p, pdo.__finspin_func_string__, mu0=mu0, mu1=mu1, mu2=mu2, mu3=mu3, mu4=mu4, chip=chip )
        
        #
        pdo_amp  = array( [pdo.IMRPhenomDAmplitude(k, pdo.model_pars, UsefulPowers(k)) for k in f] )
        
        #
        return pdo_amp
    
    #
    # def template_dphi( f, nu4=0, nu5=0, nu6=0 ):
    def template_dphi( f, nu4=0, nu5=0, nu6=0, zeta2=0 ):
        
        # Recompute internal parameters -- including connection coefficients for continuous inspiral to merger-ringdown
        pdo.model_pars = pdo.compute_model_parameters( pdo.p, pdo.__finspin_func_string__, nu4=nu4, nu5=nu5, nu6=nu6, zeta2=zeta2, chip=chip )
        
        #
        pdo_dphi  = array( [pdo.IMRPhenomDPhaseDerivFrequencyPoint(k,pdo.p['eta'], pdo.model_pars) for k in f] )
        
        #
        shift = max( pdo_dphi[ (f>=0.03) & (f<=0.12) ] )
        pdo_dphi -= shift
        pdo_dphi *= -1
        
        #
        return pdo_dphi
    
    #
    def template_phi( f, nu4=0, nu5=0, nu6=0, zeta2=0 ):
        
        #
        from numpy import argmax
        from positive import spline_diff
        
        # Recompute internal parameters -- including connection coefficients for continuous inspiral to merger-ringdown
        pdo.model_pars = pdo.compute_model_parameters( pdo.p, pdo.__finspin_func_string__, nu4=nu4, nu5=nu5, nu6=nu6, zeta2=zeta2, chip=chip )
        
        #
        pdo_phi  = array( [pdo.IMRPhenomDPhase(k,pdo.p['eta'], pdo.model_pars, UsefulPowers(k)) for k in f] )
        
        #
        dphi = spline_diff(f,pdo_phi)
        shift_id = sum(f<0.03) + argmax( dphi[ (f>=0.03) & (f<=0.12) ] )
        shift = dphi[shift_id]
        pdo_phi -= f*shift
        pdo_phi -= pdo_phi[shift_id]
        pdo_phi *= -1
        
        #
        return pdo_phi
        
    #
    return template_amp, template_dphi, template_phi

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
 
 