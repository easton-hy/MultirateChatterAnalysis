import numpy as np
from . import make_global as g

# Semi-Discretization Method parameters equivalent to data/make_SDM.m
class SDM:
    stx = 50  # steps of spindle speed
    sty = 50  # steps of depth of cut
    k = round(2*np.pi/g.O_ST*8*1e3/g.FLUTES)  # tau = k dt
    intk = 20  # tau_Kf = intk * dt
    wa = 0.5
    wb = 0.5
