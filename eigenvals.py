from matplotlib import ticker
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'font.size': 14})

fig, ax = plt.subplots(figsize = (6., 6.))

jn = 1.
# Partially trapped parameters
omegaI = 2 * np.pi * 1
vtheta = 2 * np.pi * 1e4
vphi = 2 * np.pi * 1e-4
xlim1 = 1e-6
xlim2 = 1e6
ylim1 = 1e-14
ylim2 = 1e6
# Gyroscope parameters
#omegaI = 2 * np.pi * 1e4
#vtheta = 2 * np.pi * 1
#vphi = 2 * np.pi * 1e-4
#xlim1 = 1e-6
#xlim2 = 1e6
#ylim1 = 1e-14
#ylim2 = 1e2

chiInv = lambda omega: np.array([[-omega ** 2 + omegaI * vtheta, -1j * jn * omegaI * omega],
                                 [1j * jn * omegaI * omega, -omega ** 2 + omegaI * vphi]])

ax.set_xlim(xlim1, xlim2)
ax.set_ylim(ylim1, ylim2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$f$\,[Hz]')
ax.set_ylabel(r'$I|\lambda|\,[\mathrm{Hz}^{-2}]$')

xcoor = np.concatenate((np.logspace(np.log10(xlim1), -2.1, 400),
                        np.linspace(10 ** -2.1, 10 ** -1.9, 100),
                        np.logspace(-1.9, 1.9, 200),
                        np.linspace(10 ** 1.9, 10 ** 2.1, 100),
                        np.logspace(2.1, np.log10(xlim2), 400)))
colors = [(0.317647, 0.654902, 0.752941), (1., 0.721569, 0.219608)]
lambdas = 1 / np.abs(np.linalg.eigvals(np.moveaxis(chiInv(2 * np.pi * xcoor), 2, 0)))
ax.plot(xcoor, np.amax(lambdas, axis = 1), color = colors[0], label = r'$I|\lambda_\mathrm{max}|$', zorder = 2.1)
ax.plot(xcoor, np.amin(lambdas, axis = 1), color = colors[1], label = r'$I|\lambda_\mathrm{min}|$', zorder = 2)
ax.legend()

class CustomTicker(ticker.LogFormatterSciNotation): 
    def __call__(self, x, pos = None): 
        if x not in [0.1, 1, 10]: 
            return ticker.LogFormatterSciNotation.__call__(self, x, pos = None) 
        else: 
            return "{x:g}".format(x = x)

ax.tick_params(which = 'both', direction = 'in')
ax.xaxis.set_major_formatter(CustomTicker())
ax.yaxis.set_major_formatter(CustomTicker())
secxax = ax.secondary_xaxis('top', zorder = 1)
secxax.tick_params(which = 'both', direction = 'in')
plt.setp(secxax.get_xticklabels(), visible = False)
secyax = ax.secondary_yaxis('right', zorder = 1)
secyax.tick_params(which = 'both', direction = 'in')
plt.setp(secyax.get_yticklabels(), visible = False)

fig.tight_layout()
#fig.show()
fig.savefig('eigenvals_trapped.pdf')
