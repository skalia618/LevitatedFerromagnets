from matplotlib import ticker
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'font.size': 14})

fig, ax = plt.subplots(figsize = (6., 6.))

jn = 1.
# Partially trapped parameters
##omegaI = 2 * np.pi * 1
##vtheta = 2 * np.pi * 1e4
##vphi = 2 * np.pi * 1e-4
##xlim1 = 1e-8
##xlim2 = 1e8
##ylim1 = 1e-14
##ylim2 = 1e6
# Gyroscope parameters
omegaI = 2 * np.pi * 1e4
vtheta = 2 * np.pi * 1
vphi = 2 * np.pi * 1e-4
xlim1 = 1e-8
xlim2 = 1e8
ylim1 = 1e-14
ylim2 = 1e2

chiInv = lambda omega: np.array([[-omega ** 2 + omegaI * vtheta, -1j * jn * omegaI * omega],
                                 [1j * jn * omegaI * omega, -omega ** 2 + omegaI * vphi]])

ax.set_xlim(xlim1, xlim2)
ax.set_ylim(ylim1, ylim2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$f$\,[Hz]')
ax.set_ylabel(r'$I|\chi_{\alpha\beta}|\,[\mathrm{Hz}^{-2}]$')

# Partially trapped shading
##ax.fill_betweenx([ylim1, ylim2], [vphi / (2 * np.pi * jn), vphi / (2 * np.pi * jn)], [jn * omegaI / (2 * np.pi), jn * omegaI / (2 * np.pi)], color = 'limegreen', alpha = 0.1)
##ax.text(np.sqrt(xlim1 * vphi / (2 * np.pi * jn)), 3e-14, 'Libration', color = 'black', ha = 'center', va = 'center', size = 'medium')
##ax.text(np.sqrt(vphi * omegaI) / (2 * np.pi), 8e-14, 'Libration/\nPrecession', color = 'black', ha = 'center', va = 'center', size = 'medium')
##ax.text(np.sqrt(xlim2 * jn * omegaI / (2 * np.pi)), 3e-14, 'Libration', color = 'black', ha = 'center', va = 'center', size = 'medium')
# Gyroscope shading
ax.fill_betweenx([ylim1, ylim2], [vphi / (2 * np.pi * jn), vphi / (2 * np.pi * jn)], [vtheta / (2 * np.pi * jn), vtheta / (2 * np.pi * jn)], color = 'limegreen', alpha = 0.1)
ax.fill_betweenx([ylim1, ylim2], [vtheta / (2 * np.pi * jn), vtheta / (2 * np.pi * jn)], [jn * omegaI / (2 * np.pi), jn * omegaI / (2 * np.pi)], color = 'limegreen', alpha = 0.2)
ax.text(np.sqrt(xlim1 * vphi / (2 * np.pi * jn)), 3e-14, 'Libration', color = 'black', ha = 'center', va = 'center', size = 'medium')
ax.text(np.sqrt(vphi * vtheta) / (2 * np.pi * jn), 8e-14, 'Libration/\nPrecession', color = 'black', ha = 'center', va = 'center', size = 'medium')
ax.text(np.sqrt(vtheta * omegaI) / (2 * np.pi), 3e-14, 'Precession', color = 'black', ha = 'center', va = 'center', size = 'medium')
ax.text(np.sqrt(xlim2 * jn * omegaI / (2 * np.pi)), 3e-14, 'Libration', color = 'black', ha = 'center', va = 'center', size = 'medium')

xcoor = np.concatenate((np.logspace(np.log10(xlim1), -2.1, 400),
                        np.linspace(10 ** -2.1, 10 ** -1.9, 100),
                        np.logspace(-1.9, -0.1, 200),
                        np.linspace(10 ** -0.1, 10 ** 0.1, 100),
                        np.logspace(0.1, 1.9, 200),
                        np.linspace(10 ** 1.9, 10 ** 2.1, 100),
                        np.logspace(2.1, np.log10(xlim2), 400)))
colors = [(0.317647, 0.654902, 0.752941), (1., 0.721569, 0.219608), (0.921569, 0.494118, 0.431373)]
chi = np.abs(np.linalg.inv(np.moveaxis(chiInv(2 * np.pi * xcoor), 2, 0)))
ax.plot(xcoor, chi[:,0,0], color = colors[0], label = r'$I|\chi_{\theta\theta}|$')
ax.plot(xcoor, chi[:,0,1], color = colors[1], label = r'$I|\chi_{\theta\phi}|$')
ax.plot(xcoor, chi[:,1,1], color = colors[2], label = r'$I|\chi_{\phi\phi}|$')
ax.legend()

class CustomTicker(ticker.LogFormatterSciNotation): 
    def __call__(self, x, pos = None): 
        if x not in [0.1, 1, 10]: 
            return ticker.LogFormatterSciNotation.__call__(self, x, pos = None) 
        else: 
            return "{x:g}".format(x = x)

ax.tick_params(which = 'both', direction = 'in')
ax.set_xticks([1e-8, 1e-4, 1, 1e4, 1e8])
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
fig.savefig('elements_gyro.pdf')
