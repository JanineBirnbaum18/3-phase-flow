import numpy as np
from trig_fund import *
from syrup_prop import *
import warnings

# filter warnings regularly ecountered in the model
warnings.filterwarnings("ignore", message="divide by zero encountered in true_divide")
warnings.filterwarnings("ignore", message="divide by zero encountered in power")
warnings.filterwarnings("ignore", message="invalid value encountered in multiply")
warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")
warnings.filterwarnings("ignore", message="invalid value encountered in less")
warnings.filterwarnings("ignore", message="overflow encountered in true_divide")

def dambreak(x, h0, theta, rho, K, tauy, n, g, ts, i):
    
    """Solve flow for a Herschel-Bulkley fluid with 2-step Runge-Kutta 
    in time, centered in space for on an inclined surface.
    
    :Input:
     - *x*     (ndarray(m+2))      - Evaluation points in x [m].
     - *h0*    (ndarray(m+2))      - Initial height at evaluation points in x [m].
     - *theta* (float)             - Slope [deg].
     - *rho*   (float)             - Fluid density [kg/m^3].
     - *K*     (float)             - Fluid consistency [Pa s].     
     - *tauy*  (float)             - Yield stress [Pa].
     - *n*     (float)             - Rheology power law exponent.
     - *g*     (float)             - Gravitational constant [m/s^2].
     - *ts*    (float)             - Final time for simulation [sec].
     - *i*     (int)               - Number of points for discretization in t.
     

    :Output:
     - (ndarray(m+2, i+2)) - Solution of ODE.
     - (ndarray(i+2)) - Evaluation times.
    """

    # Discretize domain 
    dx = x[1:] - x[:-1]
    dt = ts/(i+1)
    t = np.linspace(0, ts, i+2)

    # Initialize solution
    h = np.zeros((i+2,x.size))    
    h[0,:] = h0
    
    # Solve
    for k in np.arange(i+1):
        
        for j, s in enumerate([0.5, 1]):
            
            hr = (h[k+j,2:] + h[k+j,1:-1])/2
            hl = (h[k+j,1:-1] + h[k+j,:-2])/2
            hxr = sind(theta) - cosd(theta)*(h[k+j,2:] - h[k+j,1:-1])/dx[1:]
            hxl = sind(theta) - cosd(theta)*(h[k+j,1:-1] - h[k+j,:-2])/dx[:-1]
            Yr  = hr - tauy/(rho*g*np.abs(hxr))
            Yl  = hl - tauy/(rho*g*np.abs(hxl))
        
            Yr[Yr<0] = 0
            Yl[Yl<0] = 0
        
            Qr = n*((np.abs(hxr))**(1/n-1)*(Yr)**(1+1/n))/((n+1)*(2*n+1))*((1+2*n)*hr - n*Yr)*hxr
            Ql = n*((np.abs(hxl))**(1/n-1)*(Yl)**(1+1/n))/((n+1)*(2*n+1))*((1+2*n)*hl - n*Yl)*hxl
        
            Qr[hxr==0] = 0
            Ql[hxl==0] = 0
    
            h[k+1,1:-1] = h[k,1:-1] - s*dt/dx[:-1]*(rho*g/K)**(1/n)*(Qr - Ql)
            
            # Boundary conditions
            h[k+1,0] = h[k+1,1]
    
    return h, t

def dambreak_u(x, z, h0, theta, rho, K, tauy, n, g, ts, i):
    
    """Solve flow for a Herschel-Bulkley fluid with 2-step Runge-Kutta 
    in time, centered in space for on an inclined surface.
    
    :Input:
     - *x*     (ndarray(m+2))      - Evaluation points in x [m].
     - *z*     (ndarray(n+2))      - Evaluation points in z [n].
     - *h0*    (ndarray(m+2))      - Initial height at evaluation points in x [m].
     - *theta* (float)             - Slope [deg].
     - *rho*   (float)             - Fluid density [kg/m^3].
     - *K*     (float)             - Fluid consistency [Pa s].     
     - *tauy*  (float)             - Yield stress [Pa].
     - *n*     (float)             - Rheology power law exponent.
     - *g*     (float)             - Gravitational constant [m/s^2].
     - *ts*    (float)             - Final time for simulation [sec].
     - *i*     (int)               - Number of points for discretization in t.
     

    :Output:
     - (ndarray(m+2, i+2)) - Solution of ODE.
     - (ndarray(i+2)) - Evaluation times.
    """

    # Discretize domain 
    dx = x[1:] - x[:-1]
    dt = ts/(i+1)
    t = np.linspace(0, ts, i+2)

    # Initialize solution
    h = np.zeros((i+2,x.size))    
    h[0,:] = h0
    Y = np.zeros((i+2,x.size))
    u = np.zeros((i+2,x.size-2,z.size))
    
    # Solve
    for k in np.arange(i+1):
        
        for j, s in enumerate([0.5, 1]):
            
            hr = (h[k+j,2:] + h[k+j,1:-1])/2
            hl = (h[k+j,1:-1] + h[k+j,:-2])/2
            hxr = sind(theta) - cosd(theta)*(h[k+j,2:] - h[k+j,1:-1])/dx[1:]
            hxl = sind(theta) - cosd(theta)*(h[k+j,1:-1] - h[k+j,:-2])/dx[:-1]
            Yr  = hr - tauy/(rho*g*np.abs(hxr))
            Yl  = hl - tauy/(rho*g*np.abs(hxl))
        
            Yr[Yr<0] = 0
            Yl[Yl<0] = 0
        
            Qr = n*((np.abs(hxr))**(1/n-1)*(Yr)**(1+1/n))/((n+1)*(2*n+1))*((1+2*n)*hr - n*Yr)*hxr
            Ql = n*((np.abs(hxl))**(1/n-1)*(Yl)**(1+1/n))/((n+1)*(2*n+1))*((1+2*n)*hl - n*Yl)*hxl
        
            Qr[hxr==0] = 0
            Ql[hxl==0] = 0
    
            h[k+1,1:-1] = h[k,1:-1] - s*dt/dx[:-1]*(rho*g/K)**(1/n)*(Qr - Ql)
            
            # Boundary conditions
            h[k+1,0] = h[k+1,1]
        # Velocity
        h0r = np.rollaxis(np.tile(Yr,(z.size,1)),0,2) - np.tile(z,(x.size-2,1))
        h0r[h0r<0] = 0
        h0l = np.rollaxis(np.tile(Yl,(z.size,1)),0,2) - np.tile(z,(x.size-2,1))
        h0l[h0l<0] = 0
        ur = -n/(n+1)*((h0r)**(1/n+1)-np.rollaxis(np.tile(Yr,(z.size,1)),0,2)**(1/n+1))*(g*rho/K)**(1/n)*(np.rollaxis(np.tile(hxr,(z.size,1)),0,2)**1/n)
        ul = -n/(n+1)*((h0l)**(1/n+1)-np.rollaxis(np.tile(Yr,(z.size,1)),0,2)**(1/n+1))*(g*rho/K)**(1/n)*(np.rollaxis(np.tile(hxl,(z.size,1)),0,2)**1/n)
        u[k,:,:] = (ur + ul)/2
        u[k,np.tile(z,(x.size-2,1))>np.rollaxis(np.tile(h[k,1:-1],(z.size,1)),0,2)] = 0
        Y[k,1:-1] = (Yr + Yl)/2
    
    return h, t, Y, u