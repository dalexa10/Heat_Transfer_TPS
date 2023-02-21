
def compute_viscosity(T):
    """
    Computes the AIR absolute viscosity as function of absolute temperature in Kelvin (Sutherland Law).
    """
    # mu_ref = 1.789e-5       # Expressions from Anderson - Hypersonic and High temperature gas dyn.
    # T_ref = 288
    # S = 110
    # mu = mu_ref * (T/T_ref)**(3/2) * (T_ref + S)/(T + S)
    mu = 1.46e-6 * (T**(3/2) / (T + 111))  # Expressions from Heisser
    return mu

def compute_conductivity(T):
    """ Computes the Air conductivity as function of the absolute temperature in Kelvin. Eqn 2.5b Heisser """
    k = 1.99e-3 * (T**(3/2) / (T + 112))
    return k

def compute_Reynolds_number(rho, U, mu, L=1):
    """
    Computes the Reynolds number. Default L is set one to compute Reynolds number per unit length
    """
    Re = rho * U * L / mu
    return Re

def compute_total_temperature(T, Ma, gam):
    """ Computes the total temperature from static temperature, Mach and gamma """
    Tt = (1 + (gam - 1) * Ma**2 / 2) * T
    return Tt

def compute_total_pressure(Ma, P, gam):
    """ Computes the total pressure from static pressure, Mach and gamma """
    Pt = ((1 + (gam - 1) * Ma**2 / 2)**(gam / (gam - 1))) * P
    return Pt

def compute_density_ideal_gas(P, T, R):
    """  Computes the gas density using ideal gas law """
    rho = P / (R * T)
    return rho
