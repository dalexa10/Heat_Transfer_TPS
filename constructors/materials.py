from constructors.data.materials_data import materials_db

class Material_Properties:
    """
    Class that contains solid material properties for a given material name
    """
    def __init__(self, mat_name):
        self.name = str(mat_name)
        if self.name in materials_db.keys():
            self.k = materials_db[self.name]['k']
            self.Cp = materials_db[self.name]['Cp']
            self.rho = materials_db[self.name]['rho']
        else:  # Not in local database
            try:
                self.k = mat_name['k']
                self.Cp = mat_name['Cp']
                self.rho = mat_name['rho']
            except NameError:
                print('Something went wrong. Check input data from material.')

    @property
    def alfa(self):
        """ Compute thermal diffusivity """
        return self.k / (self.rho * self.Cp)


class Mesh:
    """
    Class that contain information of layer discretization.
    """
    def __init__(self):
        self.NoN = None
        self.dx = None


class Layer_Geometry:
    """
    Contains primary layer's geometrical design variables
    Add other attributes if necessary
    """
    def __init__(self, geom_dict):
        self.thick = geom_dict['thick']
        self.mesh = Mesh()


class Single_Layer:
    """
    Single layer thermal protection system class.
    Instantiate from a set of keyword arguments
    """
    def __init__(self, threshold=5e-5, **kwargs):
        self.mat = Material_Properties(kwargs['material'])
        self.geom = Layer_Geometry(kwargs['geometry'])
        self._discretize_layer(threshold)
        self.tau = None

    def _discretize_layer(self, threshold):
        NoN = 5  # Number of nodes in layer
        dx = 1  # Length from node to node
        while dx >= threshold:
            NoN *= 2
            dx = self.geom.thick / NoN
        self.geom.mesh.dx = dx
        self.geom.mesh.NoN = NoN

class Multi_Layer_TPS:
    """
    Multi layer thermal protection system class.
    It is a composition of the Single layer TPS class
    """
    def __init__(self, TPS_dict):
        self.n_layers = 0
        self.layer = []
        self.NoN = 1  # Initialize number of nodes in discretization
        for k, properties in TPS_dict.items():
            layer_i = Single_Layer(**properties)
            self.layer.append(layer_i)
            self.n_layers += 1
            self.NoN += layer_i.geom.mesh.NoN
        self.nd_count_ls = self._interface_position_calculator()

    def _interface_position_calculator(self):
        nd_count_ls = []
        nd_count = 1
        for j in range(self.n_layers):
            nd_count += self.layer[j].geom.mesh.NoN
            nd_count_ls.append(nd_count)
        return nd_count_ls


if __name__ == '__main__':

    # ------------------------------------------------------------------------
    # Example of instantiation of a mulilayer thermal protection system object
    # Taken from: Thermal Management in a Scramjet Powered Hypersonic Cruise Vehicle by Marley
    # Table 3.1
    # ------------------------------------------------------------------------
    # Load dictionary data
    from data.TPS_layouts_data import TPS_dict_B

    TPS_inst = Multi_Layer_TPS(TPS_dict_B)
    print('Number of nodes in the TPS is {:d}'.format(TPS_inst.NoN))
