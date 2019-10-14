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

def dambreak_T(x, h0, T0, T_lft, theta, rho, K, tauy, n, g, k_T, cp, T_atm, h_atm, e, ts, i):
    
    """Solve flow for a Herschel-Bulkley fluid with 2-step Runge-Kutta 
    in time, centered in space for spatial and temporal variation in 
    density, consistency, yield stress, and power law exponent on an 
    inclined surface.
    
    :Input:
     - *x*     (ndarray(m+2))      - Evaluation points in x [m].
     - *h0*    (ndarray(m+2))      - Initial height at evaluation points in x [m].
     - *T0*    (ndarray(m+2))      - Initial temperature at evaluation points in x [deg C].
     - *T_lft* (ndarray(i+2))      - Temperature at left boundary [deg C].
     - *theta* (ndarray(m+2))      - Slope in degrees at evaluation points in x [deg].
     - *rho*   (ndarray(m+2, i+2)) - Fluid density [kg/m^3].
     - *K*     (ndarray(m+2, i+2)) - Fluid consistency [Pa s].     
     - *tauy*  (ndarray(m+2, i+2)) - Yield stress [Pa].
     - *n*     (ndarray(m+2, i+2)) - Rheology power law exponent.
     - *g*     (float)             - Gravitational constant [m/s^2].
     - *k_T*   (float)             - Thermal conductivity [W/(m K)].
     - *cp*    (float)             - Specific heat capacity [J/(kg K)].
     - *T_atm* (float)             - Air temperature [deg C].
     - *h_atm* (float)             - Air convective ehat transfer coefficient [W/(m^2 K)].
     - *e*     (float)             - Emissivity.
     - *ts*    (float)             - Final time for simulation [sec].
     - *i*     (int)               - Number of points for discretization in t.
     

    :Output:
     - (ndarray(m+2, i+2)) - Solution of ODE.
     - (ndarray(i+2)) - Evaluation times.
     - (ndarray(m+2, i+2)) - Yield Surface Y.
    """

    # Discretize domain 
    dx = x[1:] - x[:-1]
    dt = ts/(i+1)
    t = np.linspace(0, ts, i+2)

    # Initialize solution
    h = np.zeros((i+2,x.size))
    T = np.zeros((i+2,x.size))
    rhoh = np.zeros((i+2,x.size))
    
    # Initial conditions
    h[0,:] = h0
    T[0,:] = T0
    rhoh[0,:] = rho[0,:]*h0
    
    sigma = 5.67 * 10**(-8) # Stefan-Boltzmann constant in [W/(m^2 K^4)]
    
    # Solve
    for k in np.arange(i+1):
        
        for j, s in enumerate([0.5, 1]):  # Loops through twice for 2-step RK
                        
            hr = (h[k+j,2:] + h[k+j,1:-1])/2
            hl = (h[k+j,1:-1] + h[k+j,:-2])/2
            hxr = sind(theta[2:]) - (h[k+j,2:]*cosd(theta[2:]) - rho[k+j,1:-1]/rho[k+j,2:]*h[k+j,1:-1]*cosd(theta[1:-1]))/dx[1:]
            hxl = sind(theta[1:-1]) - (h[k+j,1:-1]*cosd(theta[1:-1]) - rho[k+j,:-2]/rho[k+j,1:-1]*h[k+j,:-2]*cosd(theta[:-2]))/dx[:-1]
            Yr  = hr - tauy[k+j,2:]/(rho[k+j,2:]*g*np.abs(hxr))
            Yl  = hl - tauy[k+j,1:-1]/(rho[k+j,1:-1]*g*np.abs(hxl))
        
            Yr[Yr<0] = 0
            Yl[Yl<0] = 0
            
            nr = n[k+j,2:]
            nl = n[k+j,1:-1]
        
            # Calculate mass flux
            Qr = rho[k+j,2:]*(rho[k+j,2:]*g/K[k+j,2:])**(1/nr)*nr*(np.abs(hxr)**(1/nr-1)*(Yr)**(1+1/nr))/((nr+1)*(2*nr+1))*((1+2*nr)*hr - nr*Yr)*hxr
            Ql = rho[k+j,1:-1]*(rho[k+j,1:-1]*g/K[k+j,1:-1])**(1/nl)*nl*(np.abs(hxl)**(1/nl-1)*(Yl)**(1+1/nl))/((nl+1)*(2*nl+1))*((1+2*nl)*hl - nl*Yl)*hxl
        
            Qr[hxr==0] = 0
            Ql[hxl==0] = 0
    
            # Update solution
            rhoh[k+1,1:-1] = rhoh[k,1:-1] - s*dt/dx[:-1]*(Qr - Ql)
            h[k+1,1:-1] =rhoh[k+1,1:-1]/rho[k+1,1:-1]
        
            # Heat   
            A = np.nan_to_num(dx[:-1]*np.sqrt(hxl**2 + 1)/(rhoh[k+j,1:-1]*cp)) # Surface area for conduction/radiation
            A = np.insert(A,0,1)
            A = np.append(A,0)
            A[h[k+j,:]==0] = 0
            
            # Calculate conductive/radiative flux
            Qrad = dt*A*(e*sigma*((T[k+j,:]+273)**4 - (T_atm+273)**4) + h_atm*(T[k+j,:] - T_atm))
            Qrad[np.abs(Qrad)>=np.abs(T[k+j,:] - T_atm)] = T[k+j,np.abs(Qrad)>=np.abs(T[k+j,:] - T_atm)] - T_atm # T cannot fall below T_atm

            Tr = (T[k+j,2:] + T[k+j,1:-1])/2
            Tl = (T[k+j,1:-1] + T[k+j,:-2])/2
            
            Txr = (T[k+j,2:] - T[k+j,1:-1])/dx[1:]
            Txl = (T[k+j,1:-1] - T[k+j,:-2])/dx[:-1]
            
            # Update solution
            T[k+1,1:-1] = np.nan_to_num((T[k,1:-1]*rhoh[k,1:-1] - s*dt/dx[:-1]*(T[k+j,1:-1]*Qr - T[k+j,:-2]*Ql))/rhoh[k+1,1:-1]) \
                          + 1/cp*s*dt/dx[:-1]*(k_T/rho[k+j,2:]*Txr - k_T/rho[k+j,1:-1]*Txl) \
                          - s*Qrad[1:-1]
            T[k+1,rhoh[k+1,:]==0] = T_atm # T = T_atm when flow has no mass
            
            # Left boudary conditions
            T[k+1,0] = T_lft[k+1]
            rhoh[k+1,0] = rhoh[k+1,1]
            h[k+1,0] = h[k+1,1]

    return h, t, T

def dambreak_T_K(x, h0, T0, T_lft, theta, tauy, n, g, k_T, cp, T_atm, h_atm, e, ts, i):
    
    """Solve flow for a Herschel-Bulkley fluid with 2-step Runge-Kutta 
    in time, centered in space for spatial and temporal variation in 
    density, consistency, yield stress, and power law exponent on an 
    inclined surface. Consistency is a function of temperature, requires
    function K = Visc(T) and rho = Density(T).
    
    :Input:
     - *x*     (ndarray(m+2))      - Evaluation points in x [m].
     - *h0*    (ndarray(m+2))      - Initial height at evaluation points in x [m].
     - *T0*    (ndarray(m+2))      - Initial temperature at evaluation points in x [deg C].
     - *T_lft* (ndarray(i+2))      - Temperature at left boundary [deg C].
     - *theta* (ndarray(m+2))      - Slope in degrees at evaluation points in x [deg].
     - *tauy*  (ndarray(m+2, i+2)) - Yield stress [Pa].
     - *n*     (ndarray(m+2, i+2)) - Rheology power law exponent.
     - *g*     (float)             - Gravitational constant [m/s^2].
     - *k_T*   (float)             - Thermal conductivity [W/(m K)].
     - *cp*    (float)             - Specific heat capacity [J/(kg K)].
     - *T_atm* (float)             - Air temperature [deg C].
     - *h_atm* (float)             - Air convective ehat transfer coefficient [W/(m^2 K)].
     - *e*     (float)             - Emissivity.
     - *ts*    (float)             - Final time for simulation [sec].
     - *i*     (int)               - Number of points for discretization in t.     

    :Output:
     - (ndarray(m+2, i+2)) - Solution of ODE.
     - (ndarray(i+2)) - Evaluation times.
     - (ndarray(m+2, i+2)) - Yield Surface Y.
    """

    # Discretize domain 
    dx = x[1:] - x[:-1]
    dt = ts/(i+1)
    t = np.linspace(0, ts, i+2)

    # Initialize solution
    h = np.zeros((i+2,x.size))
    rhoh = np.zeros((i+2,x.size))
    rho = np.zeros((i+2,x.size))
    T = np.zeros((i+2,x.size))
    K = np.zeros((i+2,x.size))
    
    # Initial Condition
    rho[0,:] = Density(T0)
    rhoh[0,:] = rho[0,:]*h0
    h[0,:] = h0
    T[0,:] = T0
    
    sigma = 5.67 * 10**(-8) # Stefan-Boltzmann constant in [W/(m^2 K^4)]
    
    # Solve
    for k in np.arange(i+1):
        
        for j, s in enumerate([0.5, 1]): # Loops through twice for two-step RK
            
            # Temperature-dependent viscosity
            K[k+j] = Visc(T[k+j])
            
            hr = (h[k+j,2:] + h[k+j,1:-1])/2
            hl = (h[k+j,1:-1] + h[k+j,:-2])/2
            hxr = sind(theta[2:]) - (h[k+j,2:]*cosd(theta[2:]) - rho[k+j,1:-1]/rho[k+j,2:]*h[k+j,1:-1]*cosd(theta[1:-1]))/dx[1:]
            hxl = sind(theta[1:-1]) - (h[k+j,1:-1]*cosd(theta[1:-1]) - rho[k+j,:-2]/rho[k+j,1:-1]*h[k+j,:-2]*cosd(theta[:-2]))/dx[:-1]
            Yr  = hr - tauy[k+j,2:]/(rho[k+j,2:]*g*np.abs(hxr))
            Yl  = hl - tauy[k+j,1:-1]/(rho[k+j,1:-1]*g*np.abs(hxl))
        
            # Yield Surface
            Yr[Yr<0] = 0
            Yl[Yl<0] = 0
            
            nr = n[k+j,2:]
            nl = n[k+j,1:-1]
        
            # Calculate mass flux
            Qr = rho[k+j,2:]*(rho[k+j,2:]*g/K[k+j,2:])**(1/nr)*nr*(np.abs(hxr)**(1/nr-1)*(Yr)**(1+1/nr))/((nr+1)*(2*nr+1))*((1+2*nr)*hr - nr*Yr)*hxr
            Ql = rho[k+j,1:-1]*(rho[k+j,1:-1]*g/K[k+j,1:-1])**(1/nl)*nl*(np.abs(hxl)**(1/nl-1)*(Yl)**(1+1/nl))/((nl+1)*(2*nl+1))*((1+2*nl)*hl - nl*Yl)*hxl
        
            Qr[hxr==0] = 0
            Ql[hxl==0] = 0
            
            # Update solution
            rhoh[k+1,1:-1] = rhoh[k,1:-1] - s*dt/dx[:-1]*(Qr - Ql)
        
            # Heat   
            A = np.nan_to_num(dx[:-1]*np.sqrt(hxl**2 + 1)/(rhoh[k+j,1:-1]*cp)) # Scaled surface area for conduction/radiation
            A = np.insert(A,0,1)
            A = np.append(A,0)
            A[h[k+j,:]==0] = 0
            
            # Heat flux through conduction and radiation
            Qrad = dt*A*(e*sigma*((T[k+j,:]+273)**4 - (T_atm+273)**4) + h_atm*(T[k+j,:] - T_atm))
            Qrad[np.abs(Qrad)>=np.abs(T[k+j,:] - T_atm)] = T[k+j,np.abs(Qrad)>=np.abs(T[k+j,:] - T_atm)] - T_atm

            Tr = (T[k+j,2:] + T[k+j,1:-1])/2
            Tl = (T[k+j,1:-1] + T[k+j,:-2])/2
            
            Txr = (T[k+j,2:] - T[k+j,1:-1])/dx[1:]
            Txl = (T[k+j,1:-1] - T[k+j,:-2])/dx[:-1]
            
            # Update solution
            T[k+1,1:-1] = np.nan_to_num((T[k,1:-1]*rhoh[k,1:-1] - s*dt/dx[:-1]*(T[k+j,1:-1]*Qr - T[k+j,:-2]*Ql))/rhoh[k+1,1:-1]) \
                          + 1/cp*s*dt/dx[:-1]*(k_T/rho[k+j,2:]*Txr - k_T/rho[k+j,1:-1]*Txl) \
                          - s*Qrad[1:-1]
            T[k+1,rhoh[k+1,:]==0] = T_atm
            T[k+1,0] = T_lft[k+1] #Left boundary condition
            
            # Find height with temperature-dependent density
            rho[k+1,:] = Density(T[k+1,:])
            h[k+1,1:-1] =rhoh[k+1,1:-1]/rho[k+1,1:-1]
            
            # Left boundary condition
            rhoh[k+1,0] = rhoh[k+1,1]
            h[k+1,0] = h[k+1,1]

    return h, t, T