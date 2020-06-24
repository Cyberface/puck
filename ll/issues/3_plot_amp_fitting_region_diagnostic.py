#!/usr/bin/env python2.7

# Setup python environment
from matplotlib.pyplot import *
from numpy import *
from positive import *
from glob import glob
from pwca import determine_data_fitting_region


# Preliminaries
# --

#
datadir = '/Users/book/KOALA/puck/ll/data/version2/'
files = glob( datadir+'*.txt' )

#
data = []
for f in files:
    alert('Loading %s'%red(f))
    data.append(loadtxt(f).T)

#
fig,ax = subplots( len(data), 1, figsize=6*figaspect(len(data)*0.618/4) )
ax = ax.flatten()

#
tight_layout(1,2,20)

#
alert('Plotting ...')
for k in range(len(data)):
    
    #
    print '.',

    #
    sca( ax[k] )

    #
    f,amp_td,amp_fd,dphi_td,dphi_fd = data[k]
    
    #
    new_data,new_knot,new_fmin,new_fmax = determine_data_fitting_region(data[k])
    plot( f[new_knot], amp_fd[new_knot], color='k', mfc='none', marker='o', ms=20, mew=4, alpha = 0.15  )
    axvline( new_fmin, color='k', ls='-',lw=8,alpha = 0.15 )
    axvline( new_fmax, color='k', ls='-',lw=8,alpha = 0.15 )

    #
    fmin,fmax = 0.03, 0.12
    mask = (f>=fmin) & (f<=fmax)
    
    #
    x = log(f[mask])
    y = smooth(dphi_td[mask]).answer
    knot = argmin(y)
    plot( f[mask][knot], amp_fd[mask][knot], color='b', mfc='none', marker='o', ms=10, mew=2  )

    #
    plot( f, amp_td, alpha=0.5, color='orange', label=r'cp-$\psi_4$-td',lw=2 )
    plot( f, amp_fd, color='k', ls='-', label=r'cp-$\psi_4$-fd',lw=2 )

    #
    xlim(0.002,0.2)
    ylim(1e-4,1e2)

    #
    axvline( fmin, color='k', ls=':',lw=1 )
    axvline( fmax, color='k', ls=':',lw=1 )
    
    fmin_new = exp(x[knot]) * 0.325
    fmax_new = exp(x[knot]) + 0.025 # * 1.315
    axvline( fmin_new, color='b', ls='--',lw=2 )
    axvline( fmax_new, color='b', ls='--',lw=2 )

    #
    xscale('log')
    yscale('log')
    
    #
    legend(ncol=2,loc=3)
    ylabel(r'$|\tilde{h}_{22}(f)|$')
    if k+1!=len(data):
        xticks([])
    else:
        xlabel('$fM$')
    title(files[k].split('/')[-1].split('.')[0],loc='left',size=12)

#
alert('Done.')
file_path = datadir+'amp_fitting_region_diagnostic.pdf'
alert('Saving batch plot to %s'%magenta(file_path))
savefig(file_path,pad_inches=2, bbox_inches = "tight")