import os
import numpy as np
from scipy.io import savemat
import matplotlib.pyplot as plt

from . import make_global as g


def save_results(dirname, tool_idx, sdm):
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    filename = f"MCA_aed_{round(g.AE/g.D*1000)}_tool_{tool_idx}"
    savemat(os.path.join(dirname, f"{filename}.mat"), {
        'sso': sdm.sso,
        'dc': sdm.dc,
        'ei': sdm.ei,
        'forcedHinf': sdm.forcedHinf
    })

    X = sdm.sso*60/2/np.pi
    Y = sdm.dc*1000

    fig, ax = plt.subplots()
    cs = ax.contour(X, Y, 20*np.log10(sdm.ei), levels=[0])
    ax.clabel(cs, inline=True)
    ax.set_xlabel('Spindle speed (rpm)')
    ax.set_ylabel('Axial depth of cut (mm)')
    fig.savefig(os.path.join(dirname, f'SLD_{filename}.png'))
    plt.close(fig)

    fig, ax = plt.subplots()
    cs = ax.contourf(X, Y, 20*np.log10(sdm.ei))
    fig.colorbar(cs)
    ax.contour(X, Y, 20*np.log10(sdm.ei), levels=[0], colors='r')
    ax.set_xlabel('Spindle speed (rpm)')
    ax.set_ylabel('Axial depth of cut (mm)')
    fig.savefig(os.path.join(dirname, f'conv_SLD_{filename}.png'))
    plt.close(fig)

    fig, ax = plt.subplots()
    cs = ax.contourf(X, Y, 20*np.log10(sdm.forcedHinf))
    fig.colorbar(cs)
    ax.set_xlabel('Spindle speed (rpm)')
    ax.set_ylabel('Axial depth of cut (mm)')
    fig.savefig(os.path.join(dirname, f'Forced_{filename}.png'))
    plt.close(fig)
