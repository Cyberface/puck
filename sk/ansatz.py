import numpy as np

def ampmrd_func(f, gamma1, gamma2, gamma3, gamma4, fdamp, frd):
    """
    amplitude merger-ringdown ansatz
    
    equation 19: (arxiv:1508.07253)
    
    added gamma4 - to modify frd
    """
    
    g3fd = gamma3*fdamp
    dfrd = f-(frd*(1+gamma4))
#     dfrd = f-(frd + gamma4)
    
    part1 = gamma1 * g3fd / ( dfrd**2 + g3fd**2 )   
    part2 = np.exp(-gamma2 * dfrd / (gamma3*fdamp))

    return part1 * part2


def dphimrd_func(f, a1, a2, a3, a4, a5, a6, fdamp, frd):
    """
    phase derivative merger-ringdown ansatz
        
    added a6 to modify fdamp
    """
    part1 = a1 + a2*f**(-2) + a3*f**(-1/4.)
    part2 = a4 * (a6*fdamp) / ((a6*fdamp)**2 + (f-a5*frd)**2)
    return part1 + part2