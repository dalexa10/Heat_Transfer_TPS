from aux_functions.thermodynamic_functions import *

class Fluid_State_CPG:
    """
    Class that contains the thermodynamic state of the fluid at a point in
    space-time in Calorically Perfect Formulation
    """
    def __init__(self, **kwargs):
        self.V = kwargs.get('V', None)          # Flow speed [m/s]
        self.Ma = kwargs.get('Ma', None)        # Flow Mach number
        self.a = kwargs.get('a', None)          # Speed of sound [m/s]
        self.T = kwargs.get('T', None)          # Static temperature [K]
        self.P = kwargs.get('P', None)          # Static pressure [Pa]
        self.R = kwargs.get('R', None)          # Gas constant
        self.Cp = kwargs.get('Cp', None)        # Gas specific heat at constant pressure
        self.Cv = kwargs.get('Cv', None)        # Gas specific heat at constant volume
        self.gam = kwargs.get('gam', None)      # Specific heats ratio
        self.rho = kwargs.get('rho', None)      # Gas density

        # Simple hack to avoid code breaking if other (weird) kwargs arguments are given
        # in class instance
        att_ls = ['V', 'Ma', 'a', 'T', 'P', 'rho', 'Sa', 'Cp', 'gam', 'R']
        [setattr(self, att_i, kwargs[att_i]) for att_i in kwargs if att_i not in att_ls]

    @property
    def mu(self):
        return compute_viscosity(self.T)
    @property
    def k(self):  # Gas conductivity
        return compute_conductivity(self.T)
    @property
    def Tt(self):
        return compute_total_temperature(self.T, self.Ma, self.gam)
    @property
    def Pt(self):
        return compute_total_pressure(self.Ma, self.P, self.gam)
    @property
    def rho(self):
        if self._rho is None:
            return compute_density_ideal_gas(self.P, self.T, self.R)
        else:
            return self._rho
    @rho.setter
    def rho(self, val):
        self._rho = val
