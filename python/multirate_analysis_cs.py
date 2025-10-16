import numpy as np
from scipy.linalg import eigvals, norm
from . import make_global as g
from .make_sdm import SDM
from .make_work import Work
from .make_tool import build_tool
from .calc_kf_cs import calc_kf
from .calc_f_g_h import calc_fgh


def multirate_analysis(tool_idx):
    tool = build_tool(tool_idx)
    sdm = SDM()
    work = Work()
    work.Kf = calc_kf(tool, sdm, work)

    sdm.sso = np.zeros((sdm.stx+1, sdm.sty+1))
    sdm.dc = np.zeros_like(sdm.sso)
    sdm.ei = np.zeros_like(sdm.sso)
    sdm.forcedHinf = np.zeros_like(sdm.sso)

    for ix in range(sdm.stx + 1):
        o = g.O_ST + ix*(g.O_FI - g.O_ST)/sdm.stx
        tool = build_tool(tool_idx)  # update due to o? same as MATLAB call
        for iy in range(sdm.sty + 1):
            w = g.W_ST + iy*(g.W_FI - g.W_ST)/sdm.sty
            F, G, H = calc_fgh(sdm, tool, work, o, w)
            sdm.sso[ix, iy] = o
            sdm.dc[ix, iy] = w
            eig_max = max(abs(eigvals(F)))
            sdm.ei[ix, iy] = eig_max**(o/g.O_NORM)
            inv_if = np.linalg.inv(np.eye(F.shape[0]) - F)
            sdm.forcedHinf[ix, iy] = norm(H @ inv_if @ G)
    return sdm
