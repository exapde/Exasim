from numpy import *
import Preprocessing, Mesh

def cubesphere(order, R0, R1, n, m):
    pc, tc = Mesh.cubemesh(n, n, m, 1)
    pc[2, :] = loginc(pc[2, :], 3)

    pc[0, :] = -pi / 4 + pi * pc[0, :] / 2
    pc[1, :] = -pi / 4 + pi * pc[1, :] / 2
    pc[2, :] = R0 + (R1 - R0) * pc[2, :]

    dgnod = Preprocessing.createdgnodes(pc, tc, 0, [], [], order)
    npe, nd, ne = dgnod.shape
    dgnod = dgnod.transpose(1, 0, 2)
    dgnod = dgnod.reshape((nd, npe * ne), order='F') 

    offe = tc.shape[1]
    offp = pc.shape[1]
    offd = npe * ne

    ph = zeros((3, 6 * offp))
    th = zeros((8, 6 * offe), dtype=int)
    dg = zeros((3, 6 * npe * ne))

    p = transform(pc)
    g = transform(dgnod)

    ph[:, 0:offp] = p
    dg[:, 0:offd] = g
    th[:, 0:offe] = tc

    R = rotmat(array([0, 0, 1]))
    for i in range(1, 4):
        p = R @ p
        g = R @ g
        ph[:, i * offp:(i + 1) * offp] = p
        dg[:, i * offd:(i + 1) * offd] = g
        th[:, i * offe:(i + 1) * offe] = tc + i * offp

    R = rotmat(array([0, 1, 0]))
    p5 = R.T @ ph[:, :offp]
    g5 = R.T @ dg[:, :offd]
    ph[:, 4 * offp:5 * offp] = p5
    dg[:, 4 * offd:5 * offd] = g5
    th[:, 4 * offe:5 * offe] = tc + 4 * offp

    p6 = R @ ph[:, :offp]
    g6 = R @ dg[:, :offd]
    ph[:, 5 * offp:6 * offp] = p6
    dg[:, 5 * offd:6 * offd] = g6
    th[:, 5 * offe:6 * offe] = tc + 5 * offp

    ph, th = fixmesh(ph.T, th.T)
    p = ph.T
    t = th.T

    dg = dg.reshape((nd, npe, 6 * ne), order='F') 
    dgnodes = dg.transpose(1, 0, 2)

    # dgnodes = dg.reshape(3, npe, 6 * ne).transpose(1, 0, 2)
    return p, t, dgnodes


def rotmat(u):
    W = array([[0, -u[2], u[1]],
               [u[2], 0, -u[0]],
               [-u[1], u[0], 0]])
    return eye(3) + W + W @ W


def fixmesh(p, t):
    snap = 100 * max(max(p, axis=0) - min(p, axis=0)) * 1024 * finfo(float).eps
    rounded = round(p / snap) * snap
    _, ix, jx = unique(rounded, axis=0, return_index=True, return_inverse=True)
    p = p[ix]
    t = jx[t]
    if t.ndim == 1:
        t = t[newaxis, :]

    pix, ix, jx = unique(t, return_index=True, return_inverse=True)
    t = jx.reshape(t.shape)
    p = p[pix]
    
    return p, t


def transform(p):
    tanth = tan(p[0, :])
    tanph = tan(p[1, :])
    x = p[2, :] / sqrt(1 + tanth**2 + tanph**2)
    return vstack((x, x * tanth, x * tanph))


def loginc(x, alpha):
    a = min(x)
    b = max(x)
    return a + (b - a) * (exp(alpha * (x - a) / (b - a)) - 1) / (exp(alpha) - 1)

