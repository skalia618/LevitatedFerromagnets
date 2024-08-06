from matplotlib import ticker
import matplotlib.pyplot as plt
import numpy as np
from params import *

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'font.size': 14})

fig, ax = plt.subplots(figsize = (6., 6.))

class CustomTicker(ticker.LogFormatterSciNotation): 
    def __call__(self, x, pos = None): 
        if x not in np.concatenate((0.1 * np.arange(1, 10), np.arange(1, 10), 10 * np.arange(1, 10))): 
            return ticker.LogFormatterSciNotation.__call__(self, x, pos = None) 
        else:
            return "{x:g}".format(x = x)

c = 2.99792e8 # in m/s
qe = 1.60218e-19 # in C
Hz_to_eV = 6.58212e-16 # converts from omega [in Hz] to mass [in eV/c^2]
rhoDM = 4.8e-5 # in J/m^3
vDM = 1e-3 * c # in m/s

SNR = 3
Tint = 365 * 86400 # in s

existing = Setup(default = 'existing')
future = Setup(default = 'future')
freefall = Setup(default = 'freefall')

xlim1 = 2 * np.pi * Hz_to_eV * 1e-3
xlim2 = 2 * np.pi * Hz_to_eV * 1e3
ylim1 = 1e-17
ylim2 = 3e-6
ax.set_xlim(xlim1, xlim2)
ax.set_ylim(ylim1, ylim2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$m_a$\,[eV]')
ax.set_ylabel(r'$g_{ae}$')

ax.tick_params(which = 'both', direction = 'in')
secxax = ax.secondary_xaxis('top', functions = (lambda x: x / (2 * np.pi * Hz_to_eV), lambda x: 2 * np.pi * Hz_to_eV * x), zorder = 1)
secxax.tick_params(which = 'both', direction = 'in')
secxax.xaxis.set_major_formatter(CustomTicker())
secxax.set_xlabel(r'$f_a$\,[Hz]')
secyax = ax.secondary_yaxis('right')
secyax.tick_params(which = 'both', direction = 'in')
plt.setp(secyax.get_yticklabels(), visible = False)

masses = np.logspace(np.log10(xlim1), np.log10(xlim2), 1000)
ax.fill_between(masses, 1.6e-13, 1, color = '0.96', zorder = 0.7)
ax.axhline(1.6e-13, color = '0.81', zorder = 0.7)
ax.text(2e-17, 3e-13, 'Red Giants', ha = 'center', va = 'center')
ax.fill_between(masses, 1.9e-12, 1, color = '0.87', zorder = 0.8)
ax.axhline(1.9e-12, color = '0.72', zorder = 0.8)
ax.text(2e-17, 3.5e-12, 'XENONnT', ha = 'center', va = 'center')
ax.fill_between(masses, 7.5e-9, 1, color = '0.84', zorder = 0.9)
ax.axhline(7.5e-9, color = '0.69', zorder = 0.9)
ax.text(1.4e-17, 1.4e-8, 'Torsion', ha = 'center', va = 'center')
comag = np.loadtxt('axelectron_constraints/comag.txt').T
ax.fill_between(comag[0], comag[1], 1, color = '0.81', zorder = 1)
ax.plot(comag[0], comag[1], color = '0.66', zorder = 1)
ax.text(3e-15, 8e-8, r'Comagnetometers', ha = 'center', va = 'center')

colors = [(0.317647, 0.654902, 0.752941), (1., 0.721569, 0.219608), (0.921569, 0.494118, 0.431373), (0.705882, 0.494118, 0.545098)]
ax.fill_between(masses, 3.5e-10 * np.sqrt(masses), 3.5e-9 * np.sqrt(masses), color = colors[3], alpha = 0.5)
ax.text(6e-15, 8e-17, 'ALP Cogenesis', rotation = 180 / np.pi * np.arctan(0.5 * np.log(xlim2 / xlim1) / np.log(ylim2 / ylim1)), ha = 'center', va = 'center')

BDM = 2 * np.sqrt(hbar * c) / ge / qe * vDM / c ** 2 * np.sqrt(2 * rhoDM)
Tcoh = lambda mass: np.minimum(Tint, c ** 2 / vDM ** 2 * 2 * np.pi * Hz_to_eV / mass) # in s
ax.plot(masses, np.sqrt(6 * SNR) / BDM / (np.sum(np.real(np.linalg.eigvals(existing.Stot(masses / Hz_to_eV))) ** -2, axis = 1) * Tcoh(masses) * Tint) ** 0.25, color = colors[0], label = r'Existing')
ax.plot(masses, np.sqrt(6 * SNR) / BDM / (np.sum(np.real(np.linalg.eigvals(future.Stot(masses / Hz_to_eV))) ** -2, axis = 1) * Tcoh(masses) * Tint) ** 0.25, color = colors[1], label = r'Future')
ax.plot(masses, np.sqrt(6 * SNR) / BDM / (np.sum(np.real(np.linalg.eigvals(freefall.Stot(masses / Hz_to_eV))) ** -2, axis = 1) * Tcoh(masses) * Tint) ** 0.25, color = colors[2], label = r'Freefall')

ax.text(3.5e-13, 4e-17, 'Axion-Electron', ha = 'center', va = 'center', fontsize = 16, bbox = dict(boxstyle = 'round', facecolor = 'white', alpha = 0.5))
#ax.legend(loc = 'lower right')

fig.tight_layout()
#fig.show()
fig.savefig('axelectron_projection.pdf')
