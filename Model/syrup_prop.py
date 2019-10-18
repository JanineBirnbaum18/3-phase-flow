# Corn syrup temperature-dependent fluid properties:

def Visc(T): 
    """Calculate viscosity of 77% sugar corn syrup at specified Temperature based on 
    best-fit curve calibrated between 23.5 - 46 deg C.
    
    :Input:
     - *T* (array_like)      - Evaluation temperatures in [deg C].
     
    :Output:
     - (array_like)          - Viscosity at temperatures T in [Pa s].

    """
    return 15603*T**(-2.44)

def Density(T): 
    """Calculate density of 77% sugar corn syrup at specified Temperature based on 
    best-fit curve calibrated between 23.5 - 46 deg C.
    
    :Input:
     - *T* (array_like)      - Evaluation temperatures in [deg C].
     
    :Output:
     - (array_like)          - Density at temperatures T in [Kg/m^3].

    """
    return 1404.6 - 0.5012*(T-20)