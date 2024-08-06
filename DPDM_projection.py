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
eV_to_kg = 1.78266e-36 # converts from eV/c^2 to kg
rhoDM = 4.8e-5 # in J/m^3
vDM = 1e-3 * c # in m/s

SNR = 3
L = 0.1 # in m
Tint = 365 * 86400 # in s

existing = Setup(default = 'existing')
future = Setup(default = 'future')
freefall = Setup(default = 'freefall')

xlim1 = 2 * np.pi * Hz_to_eV * 1e-3
xlim2 = 2 * np.pi * Hz_to_eV * 1e3
ylim1 = 1e-14
ylim2 = 1
ax.set_xlim(xlim1, xlim2)
ax.set_ylim(ylim1, ylim2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r"$m_{A'}$\,[eV]")
ax.set_ylabel(r'$\varepsilon$')
ax.yaxis.set_major_formatter(CustomTicker())

ax.tick_params(which = 'both', direction = 'in')
secxax = ax.secondary_xaxis('top', functions = (lambda x: x / (2 * np.pi * Hz_to_eV), lambda x: 2 * np.pi * Hz_to_eV * x), zorder = 1)
secxax.tick_params(which = 'both', direction = 'in')
secxax.xaxis.set_major_formatter(CustomTicker())
secxax.set_xlabel(r"$f_{A'}$\,[Hz]")
secyax = ax.secondary_yaxis('right')
secyax.tick_params(which = 'both', direction = 'in')
plt.setp(secyax.get_yticklabels(), visible = False)

darkages = np.loadtxt('DPDM_constraints/darkages.txt').T
ax.fill_between(darkages[0], darkages[1], 1, color = '0.98', zorder = 0.7)
ax.plot(darkages[0], darkages[1], color = '0.83', zorder = 0.7)
ax.text(2.5e-13, 2e-12, r"Resonant $A'\rightarrow\gamma$", ha = 'center', va = 'center')
leoT = np.loadtxt('DPDM_constraints/leoT.txt').T
ax.fill_between(leoT[0], leoT[1], 1, color = '0.95', zorder = 0.8)
ax.plot(leoT[0], leoT[1], color = '0.8', zorder = 0.8)
ax.text(1e-12, 2e-9, r'Leo T', ha = 'center', va = 'center')
firas = np.loadtxt('DPDM_constraints/firas.txt').T
ax.fill_between(firas[0], firas[1], 1, color = '0.92', zorder = 0.9)
ax.plot(firas[0], firas[1], color = '0.77', zorder = 0.9)
ax.text(2e-13, 4e-6, r"FIRAS $\gamma\rightarrow A'$", ha = 'center', va = 'center')
sqsn = np.loadtxt('DPDM_constraints/amails.txt').T
ax.fill_between(np.concatenate((sqsn[0,::1000], [sqsn[0,-1]], [sqsn[0,-1]])), np.concatenate((sqsn[1,::1000], [sqsn[1,-1]], [1])), 1, color = '0.86', zorder = 0.9)
ax.plot(np.concatenate((sqsn[0,::1000], [sqsn[0,-1]], [sqsn[0,-1]])), np.concatenate((sqsn[1,::1000], [sqsn[1,-1]], [1])), color = '0.71', zorder = 0.9)
ax.text(2e-13, 2e-2, r'AMAILS', ha = 'center', va = 'center')
snipehunt = np.loadtxt('DPDM_constraints/snipehunt_DPDM.txt').T
ax.fill_between(2 * np.pi * Hz_to_eV * snipehunt[0], snipehunt[1], 1, color = '0.81', zorder = 1)
ax.plot(2 * np.pi * Hz_to_eV * snipehunt[0], snipehunt[1], color = '0.66', zorder = 1)
ax.text(7e-15, 8e-5, 'SNIPE\nHunt', ha = 'center', va = 'center')
supermag = np.loadtxt('DPDM_constraints/supermag_DPDM.txt').T
ax.fill_between(2 * np.pi * Hz_to_eV * supermag[0], supermag[1], 1, color = '0.81', zorder = 1)
ax.plot(2 * np.pi * Hz_to_eV * supermag[0], supermag[1], color = '0.66', zorder = 1)
ax.text(1.7e-17, 4e-3, 'SuperMAG', ha = 'center', va = 'center')

colors = [(0.317647, 0.654902, 0.752941), (1., 0.721569, 0.219608), (0.921569, 0.494118, 0.431373)]
masses = np.logspace(np.log10(xlim1), np.log10(xlim2), 1000)
BDM = lambda mass: np.sqrt(2 * mu0 * rhoDM) * c / hbar * mass * eV_to_kg * L # in T
Tcoh = lambda mass: np.minimum(Tint, c ** 2 / vDM ** 2 * 2 * np.pi * Hz_to_eV / mass) # in s
ax.plot(masses, np.sqrt(6 * SNR) / BDM(masses) / (np.sum(np.real(np.linalg.eigvals(existing.Stot(masses / Hz_to_eV))) ** -2, axis = 1) * Tcoh(masses) * Tint) ** 0.25, color = colors[0], label = r'Existing')
ax.plot(masses, np.sqrt(6 * SNR) / BDM(masses) / (np.sum(np.real(np.linalg.eigvals(future.Stot(masses / Hz_to_eV))) ** -2, axis = 1) * Tcoh(masses) * Tint) ** 0.25, color = colors[1], label = r'Future')
ax.plot(masses, np.sqrt(6 * SNR) / BDM(masses) / (np.sum(np.real(np.linalg.eigvals(freefall.Stot(masses / Hz_to_eV))) ** -2, axis = 1) * Tcoh(masses) * Tint) ** 0.25, color = colors[2], label = r'Freefall')

ax.text(3.6e-17, 6e-14, 'Dark Photon', ha = 'center', va = 'center', fontsize = 16, bbox = dict(boxstyle = 'round', facecolor = 'white', alpha = 0.5))
ax.legend(loc = (0.02, 0.2))

fig.tight_layout()
#fig.show()
fig.savefig('DPDM_projection.pdf')
