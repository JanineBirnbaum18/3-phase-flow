# import packages
import numpy as np
import pandas as pd
import sys
sys.path.insert(0,"C:/Users/Janine/.local/share/emcee-master")
import emcee as mc
import corner

from trig_fund import *
from dambreak import dambreak
from syrup_prop import Visc


vid = sys.argv[1]

# load video data

meta = pd.read_hdf(vid + '.h5', key='meta')
dat = pd.read_hdf(vid + '.h5', key='df')
try:
    meta.K_post
    K_init = meta.K_post[0]
    tauy_init = meta.tauy_post[0]
    n_init = meta.n_post[0]
except AttributeError:
    K_init = meta.K_guess[0]
    tauy_init = meta.tauy_guess[0]
    n_init = meta.n_guess[0]

K_fluid = meta.K_fluid[0]
H = meta.H[0]
L = meta.L[0]
rho = meta.Rho[0]*1000
theta = meta.Slope[0]
T = meta.Temp[0]
phi_gas = meta.Phi_gas[0]
phi_solid = meta.Phi_solid[0]

t0 = 0 # Time of dam release [sec]
tloc = int(t0/(dat.Time[1] - dat.Time[0]))
t_pos = dat.Time.iloc[tloc:] - dat.Time.iloc[tloc]
x_pos = dat.X_pos.iloc[tloc:]

# Setup model for MCMC

def lnlike(mu, x, y, yerr):
    K, tauy, n, theta = mu
    h, t = dambreak(x_grid, h0, theta, rho, K, tauy, n, 9.81, t_pos.iloc[-1], 50*t_pos.size)
    X = np.zeros_like(t_pos)
    t_match = np.zeros_like(t_pos)

    threshold = 0.0001
    for j in np.arange(t_pos.shape[0]):
        X[j] = x_grid[np.max(np.nonzero(h[j*50,:]>threshold)) + 1] # find node with h>threshold
        t_match[j] = t[j*50]
    ind = np.nonzero(X[1:] - X[:-1])[0]
    model = np.interp(t_match,t_match[ind],X[ind])
    
    inv_sigma2 = 1.0/(yerr**2)
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(2*np.pi*yerr**2)))
    
# Initial guesses for MCMC
K_min = 0.5*K_init
K_max = 1.5*K_init

tauy_init += 0.1
tauy_min = 0
tauy_max = 50

n_min = 0.7
n_max = 2.5

theta_sig = 0.5
theta_range = np.linspace(-2,2,100)
def p_theta(theta):
    return 1/np.sqrt(2*np.pi*theta_sig**2)*np.exp(-theta**2/(2*theta_sig**2))

init_pos = [K_init, tauy_init, n_init, 0]

# load data to MCMC
x = t_pos
y = x_pos + L

# Define model resolution for MCMC
x_grid = np.linspace(0,y[-1]*1.2,52)
dx = x[1] - x[0]
yerr = 2*dx
h0 = np.zeros_like(x_grid)
h0[x_grid<0.2] = H

# Uniform prior distribution
def lnprior(mu):
    K, tauy, n, theta = mu
    if K_min < K < K_max and tauy_min < tauy < tauy_max and n_min < n < n_max:
        return np.log(p_theta(theta))
    return -np.inf

def lnprob(mu, x, y, yerr):
    lp = lnprior(mu)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(mu, x, y, yerr)
    
# Create walkers
ndim, nwalkers = 4, 16
mu_pos = [init_pos + [10, 1e-2, 1e-2, theta_sig]*np.random.randn(ndim) for i in range(nwalkers)]

# Run MCMC 
steps = 500
sampler = mc.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
sampler.run_mcmc(mu_pos, steps, progress=True);

s = 150 # ignore initial samples
samples = sampler.chain[:, s:, :].reshape((-1, ndim))
samples_full = sampler.chain[:, :, :].reshape((-1, ndim))

K_fit = samples[:,0]; tauy_fit = samples[:,1]; n_fit = samples[:,2]; theta_fit = samples[:,3]

# Find mode of each posterior distribution
N, bin_edges = np.histogram(K_fit, bins=20)
K_post = np.mean([bin_edges[np.argmax(N)], bin_edges[np.argmax(N)+1]]) 
N, bin_edges = np.histogram(tauy_fit, bins=20)
tauy_post = np.mean([bin_edges[np.argmax(N)], bin_edges[np.argmax(N)+1]]) 
N, bin_edges = np.histogram(n_fit, bins=20)
n_post = np.mean([bin_edges[np.argmax(N)], bin_edges[np.argmax(N)+1]]) 
N, bin_edges = np.histogram(theta_fit, bins=20)
theta_post = np.mean([bin_edges[np.argmax(N)], bin_edges[np.argmax(N)+1]]) 

# Find standard deviation
K_std = np.std(K_fit)
tauy_std = np.std(tauy_fit)
n_std = np.std(n_fit)
theta_std = np.std(theta_fit)

# Save samples
meta['K_post'] = K_post
meta['K_r'] = K_post/K_fluid
meta['tauy_post'] = tauy_post
meta['n_post'] = n_post
meta['K_std'] = K_std
meta['tauy_std'] = tauy_std
meta['n_std'] = n_std

meta.to_hdf(title_str + '.h5', key='meta', mode='a')

samples = pd.DataFrame({'K':K_fit, 'tauy':tauy_fit, 'n':n_fit, 'full':samples_full})
samples.to_hdf(title_str + '.h5', key='s', mode='a')