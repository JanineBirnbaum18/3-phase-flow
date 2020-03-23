# import packages
import numpy as np
import copy
import sys
import pandas as pd

from trig_fund import *
from dambreak import dambreak
from syrup_prop import Visc

vid = sys.argv[1]

path_to_file = r'/local/data/lava/janine/data/'
meta = pd.read_hdf(path_to_file + vid + '.h5', key='meta')
dat = pd.read_hdf(path_to_file + vid + '.h5', key='df')  

# Output: "A_EnKF"
def EnKF(A_EnKF, data, dataerr):
    # A matrix
    A = copy.deepcopy(A_EnKF)
    # Number of the measurement
    M = len(data)
    # Number of the ensemble
    N = A_EnKF.shape[1]
    # Number of the measurements + parameters
    Ndim = A_EnKF.shape[0]
    # Number of the parameters
    Np = Ndim - M

    # 1_N matrix
    OneN = np.ones([N, N]) / N

    # A bar matrix: A_bar = A OneN
    A_bar = A @ OneN
    # A' = A - A_bar
    A_prm = A - A_bar

    # H matrix
    H1 = np.zeros([Np, M])
    H2 = np.identity(M)
    H = np.vstack((H1, H2)).transpose()

    # Measurement Matrix
    d = np.kron(np.ones((N, 1)), data).transpose()

    # E-matrix: method I (using uncertainty)
    E = np.kron(np.ones((N, 1)), dataerr).transpose()

    # measurement + pertubation
    D = d + E

    # DHA = D - HA
    DHA = D - H @ A

    # HApE = HA' + E
    HApE = H @ A_prm + E

    # Singular value decomposition
    U, S, V = np.linalg.svd(HApE)
    SIGinv = (S @ S.T) ** (-1)
    
    # calculate the analysis ensemble
    X1 = SIGinv * U.transpose()
    X2 = X1 @ DHA
    X3 = U @ X2
    X4 = (H @ A_prm).transpose() @ X3
    Aa = A + A_prm @ X4

    return Aa
    
def dambreak_smooth(mu, res=50):
    K, tauy, n, theta = mu
    h, t = dambreak(x_grid, h0, theta, rho, K, tauy, n, 9.81, t_pos[-1], res*t_pos.size)
    X = np.zeros_like(t_pos)
    t_match = np.zeros_like(t_pos)

    threshold = 0.0001
    for j in np.arange(len(t_pos)):
        X[j] = min(x_grid[np.max(np.nonzero(h[j*res,:]>threshold)) + 1],len(x_grid)-1) # find node with h>threshold
        t_match[j] = t[j*res]
    ind = np.nonzero(X[1:] - X[:-1])[0]
    model = np.interp(t_match,t_match[ind],X[ind])
    
    return model

# load video data

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
t_pos = (dat.Time.iloc[tloc:] - dat.Time.iloc[tloc]).values
x_pos = (dat.X_pos.iloc[tloc:]).values

# number of ensemble
nens = 300
# number of iteration
nitr = 5

# number of the parameters
npar = 4
# Initial parameter guess
# point source location
# Consistency, K (Pas)
Kc0 = np.random.standard_normal(nens) * K_init/10 + K_init
# yield strength, tau_y (Pa)
tauyc0 = np.random.standard_normal(nens) * tauy_init/10 + tauy_init
# flow index, n
nc0 = np.random.standard_normal(nens) * n_init/20 + n_init
# slope (degrees)
thetac0 = np.random.standard_normal(nens)

# number of measurements
nobs = len(t_pos)
# allocate the A matrix
A = np.zeros([npar + nobs, nens])
# storing the parameters into the A matrix
A[0, :] = Kc0
A[1, :] = tauyc0
A[2, :] = nc0
A[3, :] = thetac0

x_grid=np.linspace(0,x_pos[-1]*1.2+0.2,52)
h0 = np.zeros_like(x_grid)
h0[x_grid<0.2]=H

# allocate the root mean square error
RMSE = np.zeros([nens])
error = np.zeros([nobs,nens,nitr])
# allocate the parameter estimation
param_est = np.zeros([npar, nens])
# storing the initial parameter
param_est[:, :] = A[0:npar, :]

# for all time steps
#for step in np.arange(1, len(t_pos)):    
    # get the measurements (Step 1)
Dat = x_pos + L
Err = 0.005*np.ones_like(x_pos)
    
# for all iterations
for iteration in np.arange(0, nitr):
    print('Iteration = {}'.format(iteration))
        
    # calculate forecast ensemble (Step 2)
    for i in np.arange(nens):
        model = dambreak_smooth(A[0:4, i]) 

        # store the data into A matrix
        A[npar:npar+nobs, i] = model
            
        # Root mean square error (RMSE)
        error[:,i,iteration] = model-Dat
        RMSE[i] = (np.sum((model - Dat) ** 2) / (len(Dat) - 1))**(1/2)
        
        #print('RMSE = {}'.format(RMSE[i].mean()))
        if (i%20 == 0):
            print('Ensemble = {}'.format(i))
      
    # EnKF analysis (Step 3 & 4)
    A = EnKF(A, Dat, Err)

# store the parameter estimation
param_est = A[0:npar, :]

EnKF_post = pd.DataFrame(param_est.transpose(), columns=['K', 'tauy', 'n', 'theta'])
RMS = pd.DataFrame((np.sum((error) ** 2, axis=0) / (len(Dat) - 1))**(1/2))

EnKF_post.to_hdf(path_to_file + vid + '.h5', key='EnKF_params', mode='a')     
RMS.to_hdf(path_to_file + vid + '.h5', key='EnKF_err', mode='a') 
    
print('All done')