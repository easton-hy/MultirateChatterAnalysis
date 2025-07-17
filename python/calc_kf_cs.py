import numpy as np
from . import make_global as g


def calc_kf(tool, sdm, work):
    k = sdm.k
    intk = sdm.intk
    flutes = tool['flutes']
    fist = tool['fist']
    fiex = tool['fiex']
    Kt = work.Kt
    Kn = work.Kn

    Kf_mat = np.zeros((2, 2, k))
    phi = np.zeros(intk)

    for i in range(k):
        dtr = 2*np.pi/flutes/k
        for j in range(flutes):
            garr = np.zeros(intk)
            for h in range(intk):
                phi[h] = (dtr*(i + (h+1)/intk) + j*2*np.pi/flutes) % (2*np.pi)
                if fist <= phi[h] <= fiex:
                    garr[h] = 1
            Kf_mat[0,0,i] += np.sum(garr*(Kt*np.cos(phi)+Kn*np.sin(phi))*np.sin(phi))/intk
            Kf_mat[0,1,i] += np.sum(garr*(Kt*np.cos(phi)+Kn*np.sin(phi))*np.cos(phi))/intk
            Kf_mat[1,0,i] += np.sum(garr*(-Kt*np.sin(phi)+Kn*np.cos(phi))*np.sin(phi))/intk
            Kf_mat[1,1,i] += np.sum(garr*(-Kt*np.sin(phi)+Kn*np.cos(phi))*np.cos(phi))/intk
    return Kf_mat
