import numpy as np

# Global parameters equivalent to data/make_Global.m
AE = 3/1000  # radial depth of cut (m)
D = 10/1000  # diameter of tool (m)
FLUTES = 4   # number of teeth
UP_OR_DOWN = -1  # down milling (-1)
W_ST = 0/1000  # starting depth of cut (m)
W_FI = 10/1000  # final depth of cut (m)
O_ST = (2*np.pi/60) * 4000  # starting spindle speed (rad/s)
O_FI = (2*np.pi/60) * 8000  # final spindle speed (rad/s)
TAU_NORM = 10e-3  # nominal period for SLD (s)

# derived nominal speed (rpm)
O_NORM = 60/TAU_NORM/FLUTES
