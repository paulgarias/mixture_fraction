import numpy as np

class mf(object):
    """
    Mixture fraction class.
    Uses cantera object gas = ct.Solution in order to work.
    The constructure requires gas (or ct.Solution), mass fraction of oxidizer stream composition Y_o,
    and mass fraction of fuel stream composition Y_f.

    Example:
    h = mf(gas,Y_o,Y_f) 
    where Y_o and Y_f are mass fraction arrays of the same size as gas.n_species. 

    Attributes:
    aa: elemental composition matrix
    gammaf: weighing factor for Bilger mixture fraction formulation
    beta_o: oxidizer stream coupling function for  mixture fraction
    beta_f: fuel stream coupling function for  mixture fraction
    Y_o: oxidizer stream mass fraction 
    Y_f: fuel stream mass fraction
    stoich: stoichiometric mixture fraction value
    """
    def __init__(self,gas,Y_o,Y_f):
        self.gas = gas ; self.Y_o = Y_o ; self.Y_f = Y_f
        self.aa = self.getelemcompmatrix()
        self.gammaf = self.getgammaf()
        self.beta_o = self.computebeta(Y_o)
        self.beta_f = self.computebeta(Y_f)
        self.stoich = -self.beta_o/(self.beta_f - self.beta_o)

    def getelemcompmatrix(self):
        aa = 0.0*np.ones((self.gas.n_elements,self.gas.n_species)) # elemental composition matrix
        for i in range(self.gas.n_elements):
            for j in range(self.gas.n_species):
                aa[i,j] = self.gas.n_atoms(self.gas.species_name(j),self.gas.element_name(i))

        return aa

    def getgammaf(self):
        gammaf = 0.0 * np.ones(self.gas.n_elements) # weighing factor
        for i in range(self.gas.n_elements):
            if self.gas.element_name(i) == 'C':
                gammaf[i] = 2.0/self.gas.atomic_weights[i]
            elif self.gas.element_name(i) == 'H':
                gammaf[i] = 1.0/(2.0*self.gas.atomic_weights[i])
            elif self.gas.element_name(i) == 'O':
                gammaf[i] = -1.0/self.gas.atomic_weights[i]
    
        return gammaf

    def getelemmassfrac(self,Y):
        elem_mass_frac = 0.0*np.ones(self.gas.n_elements)
        for i in range(self.gas.n_elements):
            elem_mass_frac[i] = 0.0
            for j in range(self.gas.n_species):
                elem_mass_frac[i] = elem_mass_frac[i] + self.aa[i,j]*self.gas.atomic_weights[i]*Y[j]/self.gas.molecular_weights[j]
        
        return elem_mass_frac
                
    def computebeta(self,Y):
        emf = self.getelemmassfrac(Y)
        beta = 0.0
        for i in range(self.gas.n_elements):
            beta = beta + self.gammaf[i]*emf[i]

        return beta
     
    def spec2mf(self,Y):
        """
        Computes the mixture fraction value from the composition Y.
        Inputs:
            Y: gas composition mass fraction
        Outputs:
            mixfr: mixture fraction for the given gas composition
        Example:
            Z = h.spec2mf(Y)
        """
        beta = self.computebeta(Y)
        mixfr = (beta - self.beta_o)/(self.beta_f - self.beta_o)
        
        return mixfr

    def mf2spec(self,mf):
        """
        Computes the gas composition for a given mixture fraction value.
        Inputs:
            mixfr: mixture fraction value
        Outputs:
            Y: gas composition mass fraction from mixture fraction value
        Example:
            Y = h.mf2spec(mixfr)
        """
        Y = 0.0*np.ones(self.gas.n_species)
        for i in range(self.gas.n_species):
            Y[i] = mf * self.Y_f[i] + (1.0-mf) * self.Y_o[i]

        return Y
