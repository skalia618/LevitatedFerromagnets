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
Hz_to_eV = 6.58212e-16 # converts from omega [in Hz] to mass [in eV/c^2]
GeV_to_J = 1.60218e-10 # converts from GeV to J
rhoDM = 4.8e-5 # in J/m^3
vDM = 1e-3 * c # in m/s

SNR = 3
const = 0.1
L = 0.1 # in m
Tint = 365 * 86400 # in s

existing = Setup(default = 'existing')
future = Setup(default = 'future')
freefall = Setup(default = 'freefall')

xlim1 = 2 * np.pi * Hz_to_eV * 1e-3
xlim2 = 2 * np.pi * Hz_to_eV * 1e3
ylim1 = 6e-18
ylim2 = 3e2
ax.set_xlim(xlim1, xlim2)
ax.set_ylim(ylim1, ylim2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$m_a$\,[eV]')
ax.set_ylabel(r'$g_{a\gamma}\,[\mathrm{GeV}^{-1}]$')
ax.yaxis.set_major_formatter(CustomTicker())

ax.tick_params(which = 'both', direction = 'in')
secxax = ax.secondary_xaxis('top', functions = (lambda x: x / (2 * np.pi * Hz_to_eV), lambda x: 2 * np.pi * Hz_to_eV * x), zorder = 1)
secxax.tick_params(which = 'both', direction = 'in')
secxax.xaxis.set_major_formatter(CustomTicker())
secxax.set_xlabel(r'$f_a$\,[Hz]')
secyax = ax.secondary_yaxis('right')
secyax.tick_params(which = 'both', direction = 'in')
plt.setp(secyax.get_yticklabels(), visible = False)

masses = np.logspace(np.log10(xlim1), np.log10(xlim2), 1000)
chandra = np.loadtxt('axphoton_constraints/chandra.txt').T
ax.fill_between(chandra[0], chandra[1], ylim2, color = '0.98', zorder = 0.7)
ax.plot(chandra[0], chandra[1], color = '0.83', zorder = 0.7)
ax.text(8e-14, 1.5e-12, r'Chandra', ha = 'center', va = 'center')
ax.fill_between(masses, 5.2304e-12, ylim2, color = '0.94', zorder = 0.8)
ax.axhline(5.2304e-12, color = '0.79', zorder = 0.8)
ax.text(4e-13, 1.7e-11, r'SN1987A', ha = 'center', va = 'center')
ax.fill_between(masses, 6.6316e-11, ylim2, color = '0.86', zorder = 0.9)
ax.axhline(6.6316e-11, color = '0.71', zorder = 0.9)
ax.text(1.8e-12, 1.8e-10, r'CAST', ha = 'center', va = 'center')
snipehunt = np.loadtxt('axphoton_constraints/snipehunt_axion.txt').T
ax.fill_between(2 * np.pi * Hz_to_eV * snipehunt[0], snipehunt[1], ylim2, color = '0.81', zorder = 1)
ax.plot(2 * np.pi * Hz_to_eV * snipehunt[0], snipehunt[1], color = '0.66', zorder = 1)
ax.text(7e-15, 5e-5, 'SNIPE\nHunt', ha = 'center', va = 'center')
supermag = np.loadtxt('axphoton_constraints/supermag_axion.txt').T
ax.fill_between(2 * np.pi * Hz_to_eV * supermag[0], supermag[1], ylim2, color = '0.81', zorder = 1)
ax.plot(2 * np.pi * Hz_to_eV * supermag[0], supermag[1], color = '0.66', zorder = 1)
ax.text(1.7e-17, 5e-5, 'SuperMAG', ha = 'center', va = 'center')

colors = [(0.317647, 0.654902, 0.752941), (1., 0.721569, 0.219608), (0.921569, 0.494118, 0.431373), (0.705882, 0.494118, 0.545098)]
ax.fill_between(masses, 8e-10 * np.sqrt(masses), 8e-9 * np.sqrt(masses), color = colors[3], alpha = 0.5)
ax.text(6e-15, 1.5e-16, 'ALP Cogenesis', rotation = 180 / np.pi * np.arctan(0.5 * np.log(xlim2 / xlim1) / np.log(ylim2 / ylim1)), ha = 'center', va = 'center')

BDM = lambda mu: const * np.sqrt(2 * hbar * c * rhoDM) * mu0 * mu / L ** 2 / GeV_to_J # in T * GeV
Tcoh = lambda mass: np.minimum(Tint, c ** 2 / vDM ** 2 * 2 * np.pi * Hz_to_eV / mass) # in s
ax.plot(masses, np.sqrt(6 * SNR) / BDM(existing.mu) / (np.sum(np.real(np.linalg.eigvals(existing.Stot(masses / Hz_to_eV))) ** -2, axis = 1) * Tcoh(masses) * Tint) ** 0.25, color = colors[0], label = r'Existing')
ax.plot(masses, np.sqrt(6 * SNR) / BDM(future.mu) / (np.sum(np.real(np.linalg.eigvals(future.Stot(masses / Hz_to_eV))) ** -2, axis = 1) * Tcoh(masses) * Tint) ** 0.25, color = colors[1], label = r'Future')
ax.plot(masses, np.sqrt(6 * SNR) / BDM(freefall.mu) / (np.sum(np.real(np.linalg.eigvals(freefall.Stot(masses / Hz_to_eV))) ** -2, axis = 1) * Tcoh(masses) * Tint) ** 0.25, color = colors[2], label = r'Freefall')

ax.text(4e-13, 7e-17, 'Axion-Photon', ha = 'center', va = 'center', fontsize = 16, bbox = dict(boxstyle = 'round', facecolor = 'white', alpha = 0.5))
#ax.legend(loc = 'lower left')

fig.tight_layout()
#fig.show()
fig.savefig('axphoton_projection.pdf')
