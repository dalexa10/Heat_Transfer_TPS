from scipy.linalg import solve
import numpy as np
import math
from constructors.materials import Multi_Layer_TPS
from aux_functions.thermodynamic_functions import compute_viscosity, \
    compute_conductivity, compute_density_ideal_gas, compute_Reynolds_number

def compute_fourier_number(alfa, dx, dt):
    """ Eqn 5-44 Cengel - Heat and Mass Transfer"""
    tau = alfa * dt / dx**2
    return tau

def compute_Prandl_number(mu, Cp, k):
    """ Eqn 9-10 Heisser - Hypersonic Airbreathing Engines """
    Pr = mu * Cp / k
    return Pr

def compute_recovery_factor(Pr, reg='lam'):
    """ Eqn 9-5 Heisser - Hypersonic Airbreathing Engines"""
    try:
        if reg == 'turb':
            return Pr ** (1/3)
        elif reg == 'lam':
            return Pr ** (1/2)
    except ValueError:
        print('Invalid flow regime')
        return None

def compute_reference_temperature_White(Te, Me, Tw):
    """
    Compute a reference temperature as function of the temperature and Mach number at the edge of BL.
    Tw corresponds to the wall temperature. Viscous Fluid Flow - White (Eqn 7.57)
    """
    T_star = Te * (0.5 + 0.039 * Me**2 + 0.5 * (Tw/Te))
    return T_star

def compute_reference_temperature_MS(T_e, Ma_e, Tw, Pr, gam, reg='lam'):
    """ Meador - Smart reference temperature method
    Compute a reference temperature as function of the temperature and Mach number at the edge of BL.
    Tw corresponds to the wall temperature.
    Anderson - Hypersonic and High Temperature Gas Dynamics Sect. 6.9.1
    """
    try:
        r = compute_recovery_factor(Pr=Pr, reg=reg)
        if reg == 'lam':
            T_star = T_e * (0.45 + 0.55 * (Tw / T_e) + 0.16 * r * ((gam - 1) /2) * Ma_e**2)
            return T_star
        elif reg == 'turb':
            T_star = T_e * (0.5 * (1 + (Tw / T_e)) + 0.16 * r * ((gam - 1) / 2) * Ma_e ** 2)
            return T_star
    except ValueError:
        print ('Invalid flow regime')
        return None

def compute_adiabatic_wall_temperature(r, Tte, Te):
    """ Eqn 16.49 Anderson - Fundamentals of Aerodynamics """
    Taw = r * (Tte - Te) + Te
    return Taw

def compute_friction_coefficient(Re, T_star, Te, regime='lam'):
    """ Eqn 9.9 Heisser - Hypersonic Airbreathing Engines"""
    if regime == 'lam':
        Cf = (0.664 * Te) / (Re**0.5 * T_star)
    else:
        Cf = (0.0574 * Te) / (Re**(1/5) * T_star)
    return Cf

def compute_Chapman_Rubesin_parameter(T_star, Te):
    """ Eqn. 7-55  Viscous Fluid Flow - White (Eqn 7.57) """
    CR = (T_star/Te) ** (-1/3)
    return CR

def compute_Stanton_number(Pr, Re, reg='lam'):
    """ Eqn 9.6 Hisser - Hypersonic Airbreathing Engines """
    try:
        if reg == 'lam':
            St = 0.332 / (Pr**(2/3) * Re**(1/2))
            return St
        elif reg == 'turb':
            St = 0.0287 / (Pr**(2/5) * Re**(1/5))
            return St
    except ValueError:
        print ('Invalid flow regime')
        return None

def compute_convective_heat(St, rho, Ue, Cp, Taw, Tw):
    qw = St * rho * Ue * Cp * (Taw - Tw)
    return qw

class Convective_Boundary_Condition:
    """
    Class to define boundary conditions of a supersonic flow over a flat plate
    for a heat transfer analysis
    It is assumed that edge conditions are the same as freestream conditions
    Inputs:
        - flow (object): thermodynamic state of a
        - Tw (float): wall temperature [K]
    """

    def __init__(self, flow, Tw, method='Meador', reg='turb'):
        if method == 'Meador':
            self.T_star = compute_reference_temperature_MS(T_e=flow.T, Ma_e=flow.Ma, Tw=300, Pr=0.71, gam=flow.gam,
                                                           reg=reg)
        else:  # White (should include exception)
            self.T_star = compute_reference_temperature_White(flow.T, flow.Ma, Tw)

        self.Cp_star = flow.Cp
        self.mu_star = compute_viscosity(self.T_star)
        self.k_star = compute_conductivity(self.T_star)
        self.rho_star = compute_density_ideal_gas(flow.P, self.T_star, flow.R)  # This assummes constant pressure over a flat plate (see Anderson - Example 6.5)
        self.Pr_star = compute_Prandl_number(self.mu_star, self.Cp_star, self.k_star)
        self.Re_star = compute_Reynolds_number(self.rho_star, flow.V, self.mu_star)
        self.St_star = compute_Stanton_number(self.Pr_star, self.Re_star, reg=reg)
        self.r = compute_recovery_factor(self.Pr_star, reg=reg)
        self.Taw = compute_adiabatic_wall_temperature(self.r, flow.Tt, flow.T)
        self.q_conv = compute_convective_heat(self.St_star, self.rho_star, flow.V, self.Cp_star, self.Taw, Tw)
        self.k_conv = self.St_star * self.rho_star * flow.V * self.Cp_star


def global_heat_transfer_matrix(TPS, flow, T_v, dt, reg='lam'):  # T_v stands for Temperature vector
    A = np.zeros([TPS.NoN, TPS.NoN])
    B = np.zeros([TPS.NoN, 1])

    # Convective Heat Object instantiation
    q_conv = Convective_Boundary_Condition(flow, T_v[0, 0], reg=reg)

    # # Convective boundary condition (Explicit time)
    # A[0, 0] = (2 * TPS.layer[0].tau + 1)
    # A[0, 1] = - 2 * TPS.layer[0].tau
    # B[0, 0] = T_v[0, 0] + q_conv * 2 * TPS.layer[0].tau * (TPS.layer[0].geom.mesh.dx / TPS.layer[0].mat.k)

    # Convective boundary condition (Implicit time)
    A[0, 0] = 2 * TPS.layer[0].tau * (1 + q_conv.k_conv * (TPS.layer[0].geom.mesh.dx / TPS.layer[0].mat.k)) + 1
    A[0, 1] = - 2 * TPS.layer[0].tau
    B[0, 0] = T_v[0, 0] + 2 * TPS.layer[0].tau * q_conv.k_conv * (TPS.layer[0].geom.mesh.dx / TPS.layer[0].mat.k) * q_conv.Taw

    # Insulated boundary condition
    A[-1, -2] = (2 * TPS.layer[-1].tau)
    A[-1, -1] = - (2 * TPS.layer[-1].tau + 1)
    B[-1, 0] = - T_v[-1, 0]

    # Conduction across material
    nd_count_ls = TPS.nd_count_ls
    nd_count = int(1)
    lay_i = int(0)

    for j in range(1, TPS.NoN - 1):
        if (nd_count == nd_count_ls[lay_i]):  # Layers interface
            Cp_rho_avg = (TPS.layer[lay_i].mat.rho * TPS.layer[lay_i].mat.Cp + TPS.layer[lay_i + 1].mat.rho *
                          TPS.layer[lay_i + 1].mat.rho) / 2
            dx_avg = (TPS.layer[lay_i].geom.mesh.dx + TPS.layer[lay_i + 1].geom.mesh.dx) / 2
            A[j, j - 1] = TPS.layer[lay_i].mat.k / TPS.layer[lay_i].geom.mesh.dx
            A[j, j] = - ((TPS.layer[lay_i].mat.k / TPS.layer[lay_i].geom.mesh.dx) +
                         (TPS.layer[lay_i + 1].mat.k / TPS.layer[lay_i + 1].geom.mesh.dx) + (Cp_rho_avg * dx_avg / dt))
            A[j, j + 1] = TPS.layer[lay_i + 1].mat.k / TPS.layer[lay_i + 1].geom.mesh.dx
            B[j, 0] = - (Cp_rho_avg * dx_avg * T_v[j, 0]) / dt
            lay_i += 1
            nd_count += 1
        else:
            A[j, j - 1] = TPS.layer[lay_i].tau
            A[j, j] = - (2 * TPS.layer[lay_i].tau + 1)
            A[j, j + 1] = TPS.layer[lay_i].tau
            B[j, 0] = - T_v[j, 0]
            nd_count += 1
    return A, B, q_conv.q_conv


def transient_heat_transfer(TPS, flow, Tw0, tf, dt, reg='lam'):
    for layer in TPS.layer:
        layer.tau = compute_fourier_number(layer.mat.alfa, layer.geom.mesh.dx, dt)
        if math.isinf(layer.tau) or math.isnan(layer.tau):
            raise Warning('Fourier number was not computed properly')

    t = 0.
    T_v = Tw0 * np.ones([TPS.NoN, 1])
    T_v_ls = []
    q_conv_ls = []
    t_ls = []
    T_v_ls.append(T_v)

    while t < tf + dt:
        A, B, q_conv = global_heat_transfer_matrix(TPS, flow, T_v, dt, reg=reg)
        if (np.isnan(B).any()):
            print('Simulation blew up at {:.3f}'.format(t))
        q_conv_ls.append(q_conv)
        T_v = solve(A, B)
        T_v_ls.append(T_v)
        t_ls.append(t)
        t += dt

    return T_v_ls, q_conv_ls, t_ls, t

if __name__ == '__main__':

    from constructors.data.TPS_layouts_data import TPS_dict_B

    TPS_instance = Multi_Layer_TPS(TPS_dict_B)
    print('The nodes in the interface of the layers are:')
    print(TPS_instance.nd_count_ls)

