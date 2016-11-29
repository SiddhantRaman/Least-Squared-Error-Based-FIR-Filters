import numpy as np
import scipy.linalg as ln


def lpfls(N,wp,ws,W):
    M = (N-1)/2
    nq = np.arange(0,2*M+1)
    nb = np.arange(0,M+1)
    q = (wp/np.pi)*np.sinc((wp/np.pi)*nq) - W*(ws/np.pi)*np.sinc((ws/np.pi)*nq)
    b = (wp/np.pi)*np.sinc((wp/np.pi)*nb)
    q[0] = wp/np.pi + W*(1-ws/np.pi) # since sin(pi*n)/pi*n = 1, not 0
    b = b.transpose()

    Q1 = ln.toeplitz(q[0:M+1])
    Q2 = ln.hankel(q[0:M+1],q[M:])
    Q = Q1+Q2

    a = ln.solve(Q,b)
    h = list(nq)
    for i in nb:
        h[i] = 0.5*a[M-i]
        h[N-1-i] = h[i]
    h[M] = 2*h[M]
    hmax = max(np.absolute(h))
    for i in nq:
        h[i] = (8191/hmax)*h[i]
    return h

def lpfls2notch(N,wp,ws,wn1,wn2,W):
    M = (N-1)/2
    nq = np.arange(0,2*M+1)
    nb = np.arange(0,M+1)
    q = (wp/np.pi)*np.sinc((wp/np.pi)*nq) - W*(ws/np.pi)*np.sinc((ws/np.pi)*nq)
    b = (wp/np.pi)*np.sinc((wp/np.pi)*nb)
    q[0] = wp/np.pi + W*(1-ws/np.pi) # since sin(pi*n)/pi*n = 1, not 0
    b = np.asmatrix(b)
    b = b.transpose()

    Q1 = ln.toeplitz(q[0:M+1])
    Q2 = ln.hankel(q[0:M+1],q[M:])
    Q = Q1+Q2

    G1 = np.cos(wn1*nb)
    G2 = np.cos(wn2*nb)
    G = np.matrix([G1,G2])

    d = np.array([0,0])
    d = np.asmatrix(d)
    d = d.transpose()

    c = np.asmatrix(ln.solve(Q,b))

    mu = ln.solve(G*ln.inv(Q)*G.transpose(),G*c - d)

    a = c - ln.solve(Q,G.transpose()*mu)
    h = np.zeros(N)
    for i in nb:
        h[i] = 0.5*a[M-i]
        h[N-1-i] = h[i]
    h[M] = 2*h[M]
    hmax = max(np.absolute(h))
    for i in nq:
        h[i] = (8191/hmax)*h[i]
    return h

def lpfls1notch(N,wp,ws,wn1,W):
    M = (N-1)/2
    nq = np.arange(0,2*M+1)
    nb = np.arange(0,M+1)
    q = (wp/np.pi)*np.sinc((wp/np.pi)*nq) - W*(ws/np.pi)*np.sinc((ws/np.pi)*nq)
    b = (wp/np.pi)*np.sinc((wp/np.pi)*nb)
    q[0] = wp/np.pi + W*(1-ws/np.pi) # since sin(pi*n)/pi*n = 1, not 0
    b = np.asmatrix(b)
    b = b.transpose()

    Q1 = ln.toeplitz(q[0:M+1])
    Q2 = ln.hankel(q[0:M+1],q[M:])
    Q = Q1+Q2

    G1 = np.cos(wn1*nb)
    G = np.matrix([G1])

    d = np.array([0])
    d = np.asmatrix(d)

    c = np.asmatrix(ln.solve(Q,b))

    mu = ln.solve(G*ln.inv(Q)*G.transpose(),G*c - d)

    a = c - ln.solve(Q,G.transpose()*mu)
    h = np.zeros(N)
    for i in nb:
        h[i] = 0.5*a[M-i]
        h[N-1-i] = h[i]
    h[M] = 2*h[M]
    hmax = max(np.absolute(h))
    for i in nq:
        h[i] = (8191/hmax)*h[i]
    return h

def bpfls(N,ws1,wp1,wp2,ws2,W):
    h = np.zeros(N)
    return h

def bpfls1notch(N,ws1,wp1,wp2,ws2,wn1,W):
    h = np.zeros(N)
    return h
