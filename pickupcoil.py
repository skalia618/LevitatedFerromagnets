from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
import numpy as np

fig = plt.figure(figsize = (6., 6.))
ax = fig.add_subplot(projection = '3d', computed_zorder = False)
ax.view_init(15., -60.)

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'font.size': 14, 'text.latex.preamble': r'\usepackage{bm}'})

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        return np.min(zs)
    
xs = np.linspace(-1., 1., 1000)
ys = np.linspace(-1., 1., 1000)
XS, YS = np.meshgrid(xs, ys)
ax.plot_surface(XS, YS, -1. + 0 * XS, color = '0.9')

ax.add_artist(Arrow3D([0., 0.], [0., 0.], [0., -1.], mutation_scale = 10., lw = 0.8, arrowstyle='<|-|>', color = 'black'))
ax.text(-0.09, 0., -0.5, r'$z_0$', color = 'black', ha = 'center', va = 'center', size = 'medium')
ax.add_artist(Arrow3D([0., -1. / np.sqrt(2)], [0., -1. / np.sqrt(2)], [0., 0.], mutation_scale = 10., lw = 0.8, arrowstyle='<|-|>', color = 'black'))
ax.text(-0.4 / np.sqrt(2), -0.4 / np.sqrt(2), 0.06, r'$R_{p,\theta}$', color = 'black', ha = 'center', va = 'center', size = 'medium')
ax.add_artist(Arrow3D([0., 0.7 / np.sqrt(2)], [0., 0.], [0., -0.7 / np.sqrt(2)], mutation_scale = 10., lw = 0.8, arrowstyle='<|-|>', color = 'black'))
ax.text(0.28, 0., -0.43, r'$R_{p,\phi}$', color = 'black', ha = 'center', va = 'center', size = 'medium')

rad = 0.2
thetas = np.linspace(0., np.pi, 100)
phis = np.linspace(0., 2 * np.pi, 100)
THETAS, PHIS = np.meshgrid(thetas, phis)
ax.plot_surface(rad * np.sin(THETAS) * np.cos(PHIS), rad * np.sin(THETAS) * np.sin(PHIS), rad * np.cos(THETAS), color = '0.8', alpha = 0.7)

ts = np.linspace(0., 1., 1000)
ax.plot(np.cos(2 * np.pi * ts), np.sin(2 * np.pi * ts), 0.02 + 0 * ts, color = 'black')
ax.plot(np.cos(2 * np.pi * ts), np.sin(2 * np.pi * ts), 0 * ts, color = 'black')
ax.plot(np.cos(2 * np.pi * ts), np.sin(2 * np.pi * ts), -0.02 + 0 * ts, color = 'black')
ax.plot(0.7 * np.cos(2 * np.pi * ts), 0.05 + 0 * ts, 0.7 * np.sin(2 * np.pi * ts), color = 'black')
ax.plot(0.7 * np.cos(2 * np.pi * ts), 0 * ts, 0.7 * np.sin(2 * np.pi * ts), color = 'black')
ax.plot(0.7 * np.cos(2 * np.pi * ts), -0.05 + 0 * ts, 0.7 * np.sin(2 * np.pi * ts), color = 'black')

ax.add_artist(Arrow3D([0., 0.5], [0., 0.], [0., 0.], mutation_scale = 20., lw = 1.5, arrowstyle='-|>', color = 'crimson'))
ax.text(0.45, 0., 0.09, r'$\bm{\hat n}$', color = 'crimson', ha = 'center', va = 'center', size = 'large')
ax.add_artist(Arrow3D([0., 0.], [0., 0.], [0., 0.5], mutation_scale = 20., lw = 1.5, arrowstyle='-|>', color = 'black'))
ax.text(0., 0., 0.55, r'$\bm{\hat\theta}$', color = 'black', ha = 'center', va = 'center', size = 'large')
ax.add_artist(Arrow3D([0., 0.], [0., 0.5], [0., 0.], mutation_scale = 20., lw = 1.5, arrowstyle='-|>', color = 'black'))
ax.text(0., 0.6, 0., r'$\bm{\hat\phi}$', color = 'black', ha = 'center', va = 'center', size = 'large')

ax.set_xlim(-0.8, 0.8)
ax.set_ylim(-0.8, 0.8)
ax.set_zlim(-1., 0.6)
ax.set_box_aspect((1., 1., 1.))
ax.axis('off')

fig.tight_layout()
#fig.show()
fig.savefig('pickupcoil.png')
