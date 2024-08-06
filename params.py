import numpy as np

kB = 1.38065e-23 # in kg*m^2/s^2/K
hbar = 1.05457e-34 # in kg*m^2/s
mu0 = 1.25664e-6 # in kg*m/s^2/A^2
ge = 2.00232
gammae = ge * 9.27401e-24 / hbar # in Hz/T
g = 9.8 # in m/s^2

class Setup:
    
    def __init__(self, default = None, **kwargs):
        if default == 'existing':
            self.R = 2e-5 # in m
            self.M = 7e5 # in A/m
            self.rho = 7400 # in kg/m^3

            self.T = 4. # in K
            self.gamma = 1e-2 # in Hz
            self.Vphi = 1e-14 # in J

            self.jn = 1.

            self.Np = 25
            self.Rp = 1e-3 # in m
            self.ap = 1e-4 # in m
            self.Ls = 8e-11 # in H
            self.Li = 1.8e-6 # in H
            self.k = 0.85
            self.kappa = np.full(2, 1000 * hbar) # in kg*m^2/s
            self.etatilde_phi = 5e-9

        if default == 'future':
            self.R = 2e-3 # in m
            self.M = 7e5 # in A/m
            self.rho = 7400 # in kg/m^3

            self.T = 0.05 # in K
            self.gamma = 2e-6 # in Hz
            self.Vratio = 1e-3

            self.jn = 1.

            self.Np = np.full(2, 6)
            self.Rp = np.full(2, 8e-3) # in m
            self.ap = np.full(2, 1e-4) # in m
            self.Ls = np.full(2, 8e-11) # in H
            self.Li = np.full(2, 1.8e-6) # in H
            self.k = np.full(2, 0.85)
            self.kappa = np.full(2, hbar) # in kg*m^2/s

        if default == 'freefall':
            self.R = 0.02 # in m
            self.M = 7e5 # in A/m
            self.rho = 7400 # in kg/m^3

            self.T = 300 # in K
            self.gamma = 1e-10 # in Hz
            self.Vtheta = 1e-10 # in J
            self.Vphi = 1e-10 # in J

            self.jn = 1.

            self.etatilde = np.full(2, 1e-5) # in sqrt(J)
            self.kappa = np.full(2, hbar) # in kg*m^2/s

        if 'R' in kwargs: self.R = kwargs['R'] # in m
        if 'M' in kwargs: self.M = kwargs['M'] # in A/m
        if 'rho' in kwargs: self.rho = kwargs['rho'] # in kg/m^3

        if 'T' in kwargs: self.T = kwargs['T'] # in K
        if 'gamma' in kwargs: self.gamma = kwargs['gamma'] # in Hz
        if 'Vtheta' in kwargs: self.Vtheta = kwargs['Vtheta'] # in J
        if 'Vphi' in kwargs: self.Vphi = kwargs['Vphi'] # in J

        if 'jn' in kwargs: self.jn = kwargs['jn']

        if 'Np' in kwargs: self.Np = kwargs['Np']
        if 'Rp' in kwargs: self.Rp = kwargs['Rp'] # in m
        if 'ap' in kwargs: self.ap = kwargs['ap'] # in m
        if 'Ls' in kwargs: self.Ls = kwargs['Ls'] # in H
        if 'Li' in kwargs: self.Li = kwargs['Li'] # in H
        if 'k' in kwargs: self.k = kwargs['k']
        if 'kappa' in kwargs: self.kappa = kwargs['kappa'] # in kg*m^2/s

        self.vol = 4 * np.pi / 3 * self.R ** 3 # in m^3
        self.mass = self.rho * self.vol # in kg
        self.mu = self.M * self.vol # in A*m^2
        self.I = 0.4 * self.mass * self.R ** 2 # in kg*m^2

        if default != 'freefall':
            self.z0 = (mu0 * self.M ** 2 * self.R ** 3 / (16 * self.rho * g)) ** 0.25 # in m
            self.Vtheta = mu0 * self.mu ** 2 / (32 * np.pi * self.z0 ** 3) # in J
            if default == 'future':
                self.Vphi = self.Vratio * self.Vtheta # in J

            self.eta = self.Np * mu0 * self.mu / (2 * self.Rp) # in Wb
            self.Mi = self.k * np.sqrt(self.Ls * self.Li) # in H
            self.Lp = self.Np ** 2 * mu0 * self.Rp * (np.log(8 * self.Rp / self.ap) - 2) # in H
            self.Ltot = self.Li + self.Lp # in H
            self.etatilde = (self.Mi / self.Ltot) * self.eta / np.sqrt(self.Ls) # in sqrt(J)
            if default == 'existing':
                self.etatilde = np.array([self.etatilde, self.etatilde_phi])
                
        self.omegaI = self.mu / gammae / self.I # in Hz
        self.vtheta = self.Vtheta * gammae / self.mu # in Hz
        self.vphi = self.Vphi * gammae / self.mu # in Hz

        self.etatilde_res = np.sqrt(4 * kB * self.I * self.gamma * self.T / self.kappa) # in sqrt(J)
        self.etatilde_broad = np.sqrt(np.array([self.Vtheta, self.Vphi])) # in sqrt(J)

    def print_freqs(self):
        print(f'omega_I: 2pi * {self.omegaI / (2 * np.pi):.1e} Hz')
        print(f'v_tt: 2pi * {self.vtheta / (2 * np.pi):.1e} Hz')
        print(f'v_pp: 2pi * {self.vphi / (2 * np.pi):.1e} Hz')

    def print_etas(self):
        for i in range(2):
            alpha = ['t', 'p']
            print(f'~eta_{alpha[i]}: {self.etatilde[i]:.1e} sqrt(J)')
            print(f'~eta^res_{alpha[i]}: {self.etatilde_res[i]:.1e} sqrt(J)')
            print(f'~eta^broad_{alpha[i]}: {self.etatilde_broad[i]:.1e} sqrt(J)')

    def Sth(self):
        return 4 * kB * self.I * self.gamma * self.T / self.mu ** 2 * np.identity(2) # in T^2/Hz

    def Sback(self):
        return np.diag(self.kappa * self.etatilde ** 2) / self.mu ** 2 # in T^2/Hz

    def chiInv(self, omega):
        return self.I * np.array([[-omega ** 2 + self.omegaI * self.vtheta, -1j * self.jn * self.omegaI * omega],
                                  [1j * self.jn * self.omegaI * omega, -omega ** 2 + self.omegaI * self.vphi]]) # in J

    def Simp(self, omega):
        return np.einsum('ijf, j, j, j, jkf -> fik', self.chiInv(omega), 1 / self.etatilde, self.kappa, 1 / self.etatilde, self.chiInv(omega)) / self.mu ** 2 # in T^2/Hz

    def Stot(self, omega):
        return self.Sth() + self.Sback() + self.Simp(omega) # in T^2/Hz
