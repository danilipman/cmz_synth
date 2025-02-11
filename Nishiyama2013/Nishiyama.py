import pylab as pl
import numpy as np
import matplotlib.colors as mc
import matplotlib as mpl
from matplotlib import rc
from matplotlib import cm
from matplotlib.patches import Polygon
from numpy import pi, sqrt, cos, sin, arctan2, linspace, log, exp, sinh, arcsinh
from astropy.io import fits
from scipy.spatial.distance import cdist
from vorbin.voronoi_2d_binning import voronoi_2d_binning
from scipy.spatial import Voronoi, voronoi_plot_2d
rc('text', usetex=True)

data = np.loadtxt('Star-Map_-10.5_wCC-ISAAC_J2000_comp.dat')

l       = data[:,0]
b       = data[:,1] 
lstar   = data[:,2]  
bstar   = data[:,3] 
N       = data[:,4] 
N_error = data[:,5] 
RA      = data[:,6] 
DE      = data[:,7]

Nl, Nb = 421, 181

l       = l.reshape(Nl,Nb)
b       = b.reshape(Nl,Nb)
lstar   = lstar.reshape(Nl,Nb)
bstar   = bstar.reshape(Nl,Nb)
N       = N.reshape(Nl,Nb)
N_error = N_error.reshape(Nl,Nb)
RA      = RA.reshape(Nl,Nb)
DE      = DE.reshape(Nl,Nb)

# plot number of stars in each bin
fig, ax = pl.subplots(figsize=(10,5))
levels  = np.linspace(0,120,121)
norm    = mc.BoundaryNorm(levels, 256)
cmap    = 'jet'
extent  = [l.max(),l.min(),b.max(),b.min()]
IM = ax.imshow(N.T,norm=norm,cmap=cmap,origin='l',extent=extent)
ax.set_xlabel(r'$l\, {\rm [deg]}$',fontsize=22)
ax.set_ylabel(r'$b\, {\rm [deg]}$',fontsize=22)
ax.set_title('Number of stars')
ax.set_aspect('equal')
ax.set_xlim(l.max(),l.min())
ax.set_ylim(b.min(),b.max())
ax.grid()
ax.tick_params(labelsize=18) 
cb = fig.colorbar(IM, ax=ax, ticks=np.arange(0,130,20))
cb.set_label(r'number of stars',fontsize=22)
fig.savefig('Nishiyama.pdf',bbox_inches='tight')
pl.close(fig)