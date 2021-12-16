"""
pyalignit.draw
--------------
Utilities for drawing pharmacophores

"""
import numpy as np

from .cpyalignit import *

HAS_MATPLOTLIB = True

try:
    import matplotlib.pyplot as plt
except ImportError:
    HAS_MATPLOTLIB = False


PHARM_COLORS = {
    FuncGroup.AROM: '#ffa500',  # orange
    FuncGroup.HDON: '#ff00ff',  # magenta
    FuncGroup.HACC: '#00ff00',  # green
    FuncGroup.LIPO: '#00ffff',  # cyan
    FuncGroup.POSC: '#0000ff',  # blue
    FuncGroup.NEGC: '#ff0000',  # red
    FuncGroup.HYBH: '#964b00',  # brown
    FuncGroup.HYBL: '#cbd123',  # yellow
    FuncGroup.EXCL: '#808080',  # grey
    FuncGroup.UNDEF: '#0f0f0f'  # black
}


def draw_pharmacophore(pharmacophore, ax=None):
    if not HAS_MATPLOTLIB:
        raise ValueError('Drawing requires matplotlib')
    if ax is None:
        ax = plt.gca(projection='3d')
    else:
        if ax.name != '3d':
            raise ValueError('Axes must be 3-dimensional')
    for point in pharmacophore:
        _draw_pharmacophore_point(point, ax)
    _set_axes_equal(ax)
    return ax


def _draw_pharmacophore_point(p, ax):
    center = [p.point.x, p.point.y, p.point.z]
    sph = _wireframe_sphere(center, p.alpha)
    color = PHARM_COLORS[p.func]
    ax.plot_wireframe(*sph, color=color, alpha=0.5)


def _wireframe_sphere(center, radius, n_meridians=20, n_latitude=None):
    if n_latitude is None:
        n_latitude = max(n_meridians/2, 4)
    u, v = np.mgrid[0:2*np.pi:n_meridians*1j, 0:np.pi:n_latitude*1j]
    sphere_x = center[0] + radius * np.cos(u) * np.sin(v)
    sphere_y = center[1] + radius * np.sin(u) * np.sin(v)
    sphere_z = center[2] + radius * np.cos(v)
    return sphere_x, sphere_y, sphere_z


def _set_axes_equal(ax):
    x_lim = ax.get_xlim3d()
    y_lim = ax.get_ylim3d()
    z_lim = ax.get_zlim3d()
    x_range = abs(x_lim[1] - x_lim[0])
    x_middle = np.mean(x_lim)
    y_range = abs(y_lim[1] - y_lim[0])
    y_middle = np.mean(y_lim)
    z_range = abs(z_lim[1] - z_lim[0])
    z_middle = np.mean(z_lim)
    plot_radius = 0.5 * max([x_range, y_range, z_range])
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
