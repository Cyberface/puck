#!/usr/bin/env python2.7

'''
#Load process and store NR data for coprecessing frame modeling.
londonl@mit.edu 2020

## Outline

1. Load simulation 
2. Compute TD Psi4 Co-Precessing Frame 
3. Symmetrize the co-precessing frame data 
4. Output the symmetrized FD amplitude and phase along with diagnostic plots
'''

# Setup python environment
from matplotlib.pyplot import *
from numpy import *
from positive import *
from nrutils import scsearch, gwylm
from glob import glob
import h5py
from os import path
import pickle
import pwca

# Preliminaries 
# --

# Define path for file IO
package_dir = parent( pwca.__path__[0] )
data_dir = package_dir + 'data/version2/'
mkdir(data_dir,verbose=True)

# Find and load data
# --

# Define simulations to load
simulation_keywords = ('q1a04t30_dPm2_T_96_552', 'q1a04t60_dPm1_T_96_552', 'q1a04t90_dP0_T_96_552', 'q1a04t120_dP0_T_96_552', 'q1a04t150_dP0_T_96_552',  'q1a08t30dPm25_T_96_408', 'q1a08t60dPm1.5_T_96_408', 'q1a08t90dPm1_T_96_408', 'q1a08t120dP0_T_96_408', 'q1a08t150dP0_T_96_408', 'q2a04t30dPm2_T_96_408', 'q2a04t60dPm1_T_96_408', 'q2a04t90dPm1_T_96_408', 'q2a04t120_T_96_408', 'q2a04t150_T_96_408', 'q2_a10_a28_ph0_th30', 'q2_a10_a28_ph0_th60', 'q2_a10_a28_ph0_th90', 'q2_a10_a28_ph0_th120', 'q2_a10_a28_ph0_th150', 'q4a04t30_T_96_360', 'q4a04t60dPm1.5D_T_96_360', 'q4a04t90_T_96_360', 'q4a04t120dP0D_T_96_360', 'q4a04t150_T_96_360', 'q4a08t30dPm5p5dRm47_T_96_360', 'q4a08t60dPm3dRm250_T_96_384', 'q4a08t90dPm1D_T_96_384', 'q4a08t120dP1_T_96_360', 'q4a08t150_T_96_360',  'q8a04t30dPm3_T_96_360', 'q8a04t60D_dPm1', 'q8a04t90dP0_T_96_360', 'q8a04t120dPp1_T_96_360', 'q8a04t150dP9_T_96_360', 'q8a08t30dPm9.35_r0.5_T_96_360', 'q8a08t60Ditm45dr075_96_360', 'q8a08t90dP0_T_96_384', 'q8a08t120dP2_r03_T_96_360', 'q8a08t150dP2_T_120_480')

# Find select simulations using scsearch
A = scsearch( keyword=simulation_keywords, notkeyword=('80_Points','ASJmodified','.0.1.0','q8precessing/q8a04t60D_dPm1/'), verbose= True, unique=True )

#
catalog_path = '/Users/book/KOALA/puck/ll/data/pwca_catalog.pickle'
alert('Saving scentry catalog list to %s'%magenta(catalog_path))
pickle.dump( A, open( catalog_path, "wb" ) )

# Let the people know.
alert('We have found %i simulations.'%len(A))

# Define loading parameters 
lmax = 2
pad = 1000
clean = True
dt = 0.5

# Load and prcess simulations 
# --

# For all sims 
for a in A:
    
    #
    txt_file_path = data_dir+'%s.txt'%a.simname
    if path.exists(txt_file_path):
        warning('It seems that %s already exists, so we\'re moving on ...'%magenta(txt_file_path),header=True)
        continue
    
    #
    alert('Processing: %s'%magenta(a.simname),header=True)
    
    # Load
    y_raw = gwylm(a,lmax=lmax,dt=dt,pad=pad,clean=clean,verbose=True)
    
    # Manage frames using dict defined below
    frame = {}
    frame['raw'] = y_raw

    # Put in initial J frame
    frame['init-j'] = y_raw.__calc_initial_j_frame__()

    # Compute TD adn FD coprecessing psi4 frames
    frame['cp-y-fd'] = frame['init-j'].__calc_coprecessing_frame__( transform_domain='fd', kind='psi4' )
    frame['cp-y-td'] = frame['init-j'].__calc_coprecessing_frame__( transform_domain='td', kind='psi4' )
    
    # Symmetrise data
    foo = dict(frame)
    for k in frame:
        if 'cp' in k:
            foo['sym-'+k] = frame[k].__symmetrize__()
    #
    frame = foo
    del foo
    
    # Produce diagnostic plots 
    def plot_amp_dphi(frame):

        #
        fig = figure( figsize=4*figaspect(0.8) )

        #
        kind = 'strain'

        #
        f = frame['cp-y-td'].f
        mask = abs(f)<0.1

        #
        ax1 = subplot(2,1,1)
        
        k = 'td'
        plot( frame['cp-y-'+k].f, frame['cp-y-'+k][2,2][kind].fd_amp, label=k+'-non-symmetrized', color='k', lw=2, alpha=0.15 )
        plot( frame['sym-cp-y-'+k].f, frame['sym-cp-y-'+k][2,2][kind].fd_amp, label=k+'-symmetrized',color='k', alpha=1, lw=1 )
        
        k = 'fd'
        plot( frame['cp-y-'+k].f, frame['cp-y-'+k][2,2][kind].fd_amp, label=k+'-non-symmetrized', color='b', lw=2, alpha=0.15 )
        plot( frame['sym-cp-y-'+k].f, frame['sym-cp-y-'+k][2,2][kind].fd_amp, label=k+'-symmetrized',color='b', alpha=1, lw=1 )
        
        title(y_raw.simname)
        ylabel(r'$|\tilde{h}_{22}|$')

        #
        yscale('log')
        xscale('log')
        legend()
        xlabel('$fM$')

        #
        kind = 'psi4'

        #
        ax2=subplot(2,1,2,sharex=ax1)
        
        k = 'td'
        plot( frame['cp-y-'+k].f, frame['cp-y-'+k][2,2][kind].fd_dphi, label=k+'-non-symmetrized', color='k', lw=2, alpha=0.15 )
        plot( frame['sym-cp-y-'+k].f, frame['sym-cp-y-'+k][2,2][kind].fd_dphi, label=k+'-symmetrized',color='k', alpha=1, lw=1 )
        
        k = 'fd'
        plot( frame['cp-y-'+k].f, frame['cp-y-'+k][2,2][kind].fd_dphi, label=k+'-non-symmetrized', color='b', lw=2, alpha=0.15 )
        plot( frame['sym-cp-y-'+k].f, frame['sym-cp-y-'+k][2,2][kind].fd_dphi, label=k+'-symmetrized',color='b', alpha=1, lw=1 )

        #
        xscale('log')
        xlim( 0.005, 0.3 )
        ylim( frame['cp-y-'+k][2,2]['strain'].fd_dphi[mask][-1]-250,frame['cp-y-'+k][2,2]['strain'].fd_dphi[mask][-1]+500 )
        xlabel('$fM$')
        ylabel(r'$\frac{d}{df}\arg(\tilde{h}_{22})$');
        
        #
        return fig,[ax1,ax2]
        
    #
    fig,ax = plot_amp_dphi(frame)
    file_path = data_dir+'%s.png'%frame['raw'].simname
    alert('Saving diagnostic plot to "%s"'%yellow(file_path))
    savefig( file_path )
    
    # Select and output amplitude and phase data
    f = frame['raw'].f
    fd_amp  = frame['sym-cp-y-fd'][2,2]['strain'].fd_amp
    fd_dphi = frame['sym-cp-y-fd'][2,2]['psi4'].fd_dphi
    td_amp  = frame['sym-cp-y-td'][2,2]['strain'].fd_amp
    td_dphi = frame['sym-cp-y-td'][2,2]['psi4'].fd_dphi
    data_array = array([ f, td_amp, fd_amp, td_dphi, fd_dphi ]).T

    #
    alert('Saving waveform data to "%s"'%yellow(txt_file_path))
    # pickle.dump( data_array, open( file_path, "wb" ) )
    savetxt( txt_file_path, data_array, header='[ f, td_amp, fd_amp, td_dphi, fd_dphi ], here td and fd refer to the frame type used; frequencies are positive, waveform info are symmetrized in the psi4 TD/FD coprecessing frame from NR simulation at %s'%frame['raw'].simdir )
    
alert('All done.')