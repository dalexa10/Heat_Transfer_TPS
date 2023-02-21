# --------------------------------------------------------
#              Thermal protection system A
# --------------------------------------------------------
TPS_dict_A = {0: {'material': {'name': 'Fiber_Carbon',
                                'k': 6.6,  # Average axial and radial
                                'Cp': 1380,
                                'rho': 120},
                  'geometry': {'thick': 12e-3}}}

# ---------------------------------------------------------
#               Thermal protection system B
# ---------------------------------------------------------
TPS_dict_B = {0: {'material': {'name': 'P2000',
                               'k': 27.7,
                               'Cp': 770.3,
                               'rho': 7116.7},
                  'geometry': {'thick': 2e-3}},
              1: {'material': {'name': 'SiO2',
                               'k': 0.033,
                               'Cp': 753.6,
                               'rho': 96.1},
                  'geometry': {'thick': 2e-3}},
              2: {'material': {'name': 'Ti',
                               'k': 18,
                               'Cp': 942,
                               'rho': 4437.1},
                  'geometry': {'thick': 2e-3}}}
