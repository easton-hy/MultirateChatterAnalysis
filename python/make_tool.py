import numpy as np
from scipy import signal
from . import make_global as g


def build_tool(tool_idx):
    s = signal.lti([], [], [])  # placeholder for lti representation
    if tool_idx == 1:
        mass = 0.4
        omegan = 1435*2*np.pi
        zeta = 0.012
        M = mass
        K = M*omegan**2
        B = 2*zeta/omegan*K
        # second order system 1/(Ms^2 + Bs + K)
        num = [1]
        den = [M, B, K]
    elif tool_idx == 2:
        K0 = 1e8
        omegani = [1925*2*np.pi, 1200*2*np.pi]
        zetai = [0.01, 0.006]
        num = [1/K0]
        den = [1]
        for on, z in zip(omegani, zetai):
            num = np.polymul(num, [on**2])
            den = np.polymul(den, [1, 2*z*on, on**2])
    elif tool_idx == 3:
        K = 0.25*1e9
        on = np.array([1052, 1648])*2*np.pi
        zeta = [0.005, 0.005]
        num = [1/K]
        den = [1]
        for o_i, z in zip(on, zeta):
            num = np.polymul(num, [o_i**2])
            den = np.polymul(den, [1, 2*z*o_i, o_i**2])
    elif tool_idx == 4:
        K = 0.58*1e9
        on = 1190*2*np.pi
        zeta = 0.0019
        num = [1/K]
        den = [1, 2*zeta*on, on**2]
    elif tool_idx == 5:
        K = 0.25*1e9
        on = np.array([1470, 1510, 1650, 1900, 2100])*2*np.pi
        zeta = [0.005]*5
        num = [1/K]
        den = [1]
        for o_i, z in zip(on, zeta):
            num = np.polymul(num, [o_i**2])
            den = np.polymul(den, [1, 2*z*o_i, o_i**2])
    else:
        raise ValueError("Unsupported tool")

    # convert to state-space, replicate for 2 DOF
    sys = signal.StateSpace(*signal.tf2ss(num, den))
    A, B, C, D = sys.A, sys.B, sys.C, sys.D
    # 2x2 plant with identical dynamics
    A = np.kron(np.eye(2), A)
    B = np.kron(np.eye(2), B)
    C = np.kron(np.eye(2), C)
    D = np.kron(np.eye(2), D)
    return {
        'plant': {'A': A, 'B': B, 'C': C, 'D': D},
        'flutes': g.FLUTES,
        'fist': angle_in_out(g.UP_OR_DOWN, g.AE/g.D)[0],
        'fiex': angle_in_out(g.UP_OR_DOWN, g.AE/g.D)[1]
    }


def angle_in_out(up_or_down, aD):
    if up_or_down == 1:
        fist = 0
        fiex = np.arccos(1 - 2*aD)
    else:
        fist = np.arccos(2*aD - 1)
        fiex = np.pi
    return fist, fiex
