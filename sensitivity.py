from matplotlib import ticker
import matplotlib.pyplot as plt
import numpy as np
from params import Setup

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'font.size': 14})

fig, ax = plt.subplots(figsize = (6., 6.))

class CustomTicker(ticker.LogFormatterSciNotation): 
    def __call__(self, x, pos = None): 
        if x not in np.concatenate((0.1 * np.arange(1, 10), np.arange(1, 10), 10 * np.arange(1, 10))): 
            return ticker.LogFormatterSciNotation.__call__(self, x, pos = None) 
        else:
            return "{x:g}".format(x = x)

default = 'existing'
xlim1 = 1
xlim2 = 3000
ylim1 = 1e-32
ylim2 = 1e-24
legend = 'upper left'
xtxt = 2.5
ytxt = 2.5e-32

##default = 'future'
##xlim1 = 1e-3
##xlim2 = 1000
##ylim1 = 1e-50
##ylim2 = 1e-28
##legend = 'upper left'
##xtxt = 4e-3
##ytxt = 1.5e-49

##default = 'freefall'
##xlim1 = 1e-5
##xlim2 = 10
##ylim1 = 1e-50
##ylim2 = 1e-28
##legend = 'upper left'
##xtxt = 2.5
##ytxt = 1.5e-49

ax.set_xlim(xlim1, xlim2)
ax.set_ylim(ylim1, ylim2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$f\,[\mathrm{Hz}]$')
ax.set_ylabel(r'$S_{BB}\,[\mathrm{T}^2/\mathrm{Hz}]$')
ax.xaxis.set_major_formatter(CustomTicker())
ax.yaxis.set_minor_locator(ticker.NullLocator())

ax.tick_params(which = 'both', direction = 'in')
secxax = ax.secondary_xaxis('top')
secxax.tick_params(which = 'both', direction = 'in')
plt.setp(secxax.get_xticklabels(), visible = False)
secxax.xaxis.set_minor_formatter(ticker.NullFormatter())
secyax = ax.secondary_yaxis('right')
secyax.tick_params(which = 'both', direction = 'in')
plt.setp(secyax.get_yticklabels(), visible = False)
secyax.yaxis.set_minor_locator(ticker.NullLocator())

colors = [(0.921569, 0.494118, 0.431373), (1., 0.721569, 0.219608), (0.317647, 0.654902, 0.752941)]
freqs = np.logspace(np.log10(xlim1), np.log10(xlim2), 1000)
setup = Setup(default = default)
ax.axhline(setup.Sth()[0, 0], color = colors[0], label = r'Thermal')
ax.axhline(setup.Sback()[0, 0], color = colors[1], label = r'Back-action')
ax.plot(freqs, np.amin(np.abs(np.linalg.eigvals(setup.Simp(2 * np.pi * freqs))), axis = 1), color = colors[2], label = r'Imprecision (min.)')
ax.plot(freqs, np.amax(np.abs(np.linalg.eigvals(setup.Simp(2 * np.pi * freqs))), axis = 1), linestyle = '--', color = colors[2], label = r'Imprecision (max.)')
ax.plot(freqs, 1 / np.sqrt(np.trace(np.linalg.matrix_power(setup.Stot(2 * np.pi * freqs), -2), axis1 = 1, axis2 = 2).real), color = 'k', label = r'Total')
if legend != None: ax.legend(loc = legend)

ax.text(xtxt, ytxt, default.capitalize(), ha = 'center', va = 'center', fontsize = 16, bbox = dict(boxstyle = 'round', facecolor = 'white', alpha = 0.5))

fig.tight_layout()
#fig.show()
fig.savefig(f'sensitivity_{default}.pdf')
