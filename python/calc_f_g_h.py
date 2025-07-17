import numpy as np
from scipy.linalg import expm


def calc_fgh(sdm, tool, work, o, w):
    k = sdm.k
    wa = sdm.wa
    wb = sdm.wb
    flutes = tool['flutes']
    Kf = work.Kf

    tau = 2*np.pi/o/flutes
    dt = tau/k

    A = tool['plant']['A']
    B = tool['plant']['B']
    C = tool['plant']['C']

    n = A.shape[0]

    Dmat = np.zeros((n*(k+1), n*(k+1)))
    Dmat[n: , n: ] = np.eye(n*k)

    F = np.eye(n*(k+1))
    G = np.zeros((n*(k+1), k))
    E = np.zeros((n*(k+1), 1))
    H = np.zeros((2, (k+1)*n))

    for i in range(k):
        Bc = w*(B@Kf[:,:,i])
        Ac2 = Bc@C
        Ac1 = A - Ac2
        P = expm(Ac1*dt)
        R = (P - np.eye(n)) @ np.linalg.inv(Ac1) @ Ac2
        Q = (P - np.eye(n)) @ np.linalg.inv(Ac1) @ Bc
        Dmat[0:n, 0:n] = P
        Dmat[0:n, n*(k+1)-2*n:n*(k+1)-n] = wa*R
        Dmat[0:n, n*(k+1)-n:n*(k+1)] = wb*R
        E[0:n,0] = Q[:,0]

        Fi0 = F.copy()
        F[0:n,:] = P@Fi0[0:n,:] + Dmat[0:n, n*(k+1)-2*n:n*(k+1)] @ Fi0[n*(k+1)-2*n:n*(k+1),:]
        F[n:n*(k+1),:] = Fi0[0:n*k,:]

        G0 = G.copy()
        G[:,i] = E[:,0]
        G[0:n,0:i+1] = P@G[0:n,0:i+1] + Dmat[0:n, n*(k+1)-2*n:n*(k+1)] @ G0[n*(k+1)-2*n:n*(k+1), 0:i+1]
        G[n:n*(k+1),0:i+1] = G0[0:n*k,0:i+1]

        Dmat.fill(0)

    for i in range(k+1):
        H[0:2, n*i:n*(i+1)] = C

    return F, G, H
