import numpy as np
import scipy.linalg as ln


def lpfls(N,wp,ws,W):
    M = (N-1)/2
    nq = np.arange(0,2*M+1)
    nb = np.arange(0,M+1)
    q = (wp/np.pi)*np.sinc((wp/np.pi)*nq) - W*(ws/np.pi)*np.sinc((ws/np.pi)*nq)
    b = (wp/np.pi)*np.sinc((wp/np.pi)*nb)
    b[0] = wp/np.pi
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
    M = (N-1)/2
    nq = np.arange(0,2*M+1)
    nb = np.arange(0,M+1)
    q = W*np.sinc(nq) - (W*ws2/np.pi) * np.sinc(nq* (ws2/np.pi)) + (wp2/np.pi) * np.sinc(nq*(wp2/np.pi)) - (wp1/np.pi) * np.sinc(nq*(wp1/np.pi)) + (W*ws1/np.pi) * np.sinc(nq*(ws1/np.pi))
    b = (wp2/np.pi)*np.sinc((wp2/np.pi)*nb) - (wp1/np.pi)*np.sinc((wp1/np.pi)*nb)
    b[0] = wp2/np.pi - wp1/np.pi
    q[0] = W - W*ws2/np.pi + wp2/np.pi - wp1/np.pi + W*ws1/np.pi # since sin(pi*n)/pi*n = 1, not 0
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

def bpfls1notch(N,ws1,wp1,wp2,ws2,wn1,W):
    h = np.zeros(N)
    return h

def bpfls2notch(N,ws1,wp1,wp2,ws2,wn1,wn2,W):
    h = np.zeros(N)
    return h

def hpfls(N,ws,wp,W):
    M = (N-1)/2
    nq = np.arange(0,2*M+1)
    nb = np.arange(0,M+1)
    b = 1 - (wp/np.pi)* np.sinc(nb * wp/np.pi)
    b[0] = 1- wp/np.pi
    q = 1 - (wp/np.pi)* np.sinc(nq * wp/np.pi) + W * (ws/np.pi) * np.sinc(nq * ws/np.pi) # since sin(pi*n)/pi*n = 1, not 0
    q[0] = b[0] + W* ws/np.pi
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
