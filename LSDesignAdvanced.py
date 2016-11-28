import numpy as np
import scipy.linalg as ln
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

Fpass = 11.0     # MHz Passband end frequency
Fstop = 15.0     #MHz Stopband start frequency
Fsamp = 50.0     # MHz Sampling Frequency
Weight = 100    # Weight of stop band error
Taps = 71      # FIR Filter Taps
Fnotch1 = 60   #Notch out freq Fnotch1
Fnotch2 = 60

n_global = np.arange(0,Taps)
h_global = np.zeros(Taps)

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
    

fig, myFilter = plt.subplots()
plt.subplots_adjust(left=0.25,bottom=0.25)
myFilter.set_xlabel('Freq(MHz)')
myFilter.set_ylabel('Magnitude(dB)')
 
h = lpfls(Taps,2*np.pi*(Fpass/Fsamp),2*np.pi*(Fstop/Fsamp),Weight)
h_quant = np.round(h)
h_global = h_quant
H = np.fft.fft(h,1024)
H_quant = np.fft.fft(h_quant,1024) 
Mag = 20*np.log10(abs(H[0:512]))
Mag_quant = 20*np.log10(abs(H_quant[0:512]))
freq = np.zeros(512)
for i in range(0,512) :
    freq[i] = Fsamp/2 * i/512
l1, = plt.plot(freq,Mag, lw=2, color='blue')
plt.axis([0, Fsamp/2, -50, 100])
l2, = plt.plot(freq,Mag_quant, lw=2, color='red')

axcolor = 'lightgoldenrodyellow'
axTaps = plt.axes([0.6, 0.125, 0.25, 0.02], axisbg=axcolor)
axFsamp = plt.axes([0.6, 0.1, 0.25,0.02], axisbg=axcolor)
axFpass = plt.axes([0.6, 0.075, 0.25, 0.02],axisbg=axcolor)
axFstop = plt.axes([0.6, 0.05,0.25, 0.02],axisbg=axcolor)
axW = plt.axes([0.1, 0.125, 0.25, 0.02],axisbg=axcolor)   
axFnotch1 = plt.axes([0.1, 0.1, 0.25, 0.02],axisbg=axcolor)  
axFnotch2 = plt.axes([0.1, 0.075, 0.25, 0.02],axisbg=axcolor)

sTaps = Slider(axTaps,'Taps',15,200,valinit=Taps)
sFsamp = Slider(axFsamp,'Samp Freq',1.0,100,valinit=Fsamp)
sFpass = Slider(axFpass,'Passband',1.0,40.0,valinit=Fpass)
sFstop = Slider(axFstop,'Stopband',1.0,45.0,valinit=Fstop)
sW = Slider(axW,'Weight',1.0,1000,valinit=Weight)
sFnotch1 = Slider(axFnotch1,'Notch 1',0.0,60,valinit=Fnotch1)
sFnotch2 = Slider(axFnotch2,'Notch 2',0.0,60,valinit=Fnotch2)

def update(val):
    Taps = int(sTaps.val)
    if(Taps%2 == 0) :
        Taps = Taps +1 
    Fsamp = sFsamp.val
    Fpass = sFpass.val
    Fstop = sFstop.val
    W = sW.val
    Fnotch1 = sFnotch1.val
    Fnotch2 = sFnotch2.val
    
    wp = 2*np.pi*(Fpass/Fsamp)
    ws = 2*np.pi*(Fstop/Fsamp)
    wn1 = 2*np.pi*(Fnotch1/Fsamp)
    wn2 = 2*np.pi*(Fnotch2/Fsamp)
    #h = np.zeros(0,Taps)
    if((Fnotch1 <= Fsamp/2) and (Fnotch2 <= Fsamp/2)) :
        h = lpfls2notch(Taps,wp,ws,wn1,wn2,W)
    elif((Fnotch1 > Fsamp/2) and (Fnotch2 <= Fsamp/2)):
        h = lpfls1notch(Taps,wp,ws,wn2,W)
    elif((Fnotch1 <= Fsamp/2) and (Fnotch2 > Fsamp/2)):
        h = lpfls1notch(Taps,wp,ws,wn1,W)
    else:
        h = lpfls(Taps,wp,ws,W)
    h_quant = np.round(h)
    global h_global
    h_global = h_quant
    H = np.fft.fft(h,1024)
    H_quant = np.fft.fft(h_quant,1024)
    Mag = 20*np.log10(abs(H[0:512]))
    Mag_quant = 20*np.log10(abs(H_quant[0:512]))
    
    freq = np.zeros(512)
    for i in range(0,512) :
        freq[i] = Fsamp/2 * i/512  
    l1.set_ydata(Mag)
    l1.set_xdata(freq)
    l2.set_ydata(Mag_quant)
    l2.set_xdata(freq)
    fig.canvas.draw_idle()
sTaps.on_changed(update)
sFsamp.on_changed(update)
sFpass.on_changed(update)
sFstop.on_changed(update)
sW.on_changed(update)
sFnotch1.on_changed(update)
sFnotch2.on_changed(update)

    
resetax = plt.axes([0.1, 0.025, 0.05, 0.04])
button1 = Button(resetax, 'Reset', color=axcolor, hovercolor='green')
def reset(event):
    sTaps.reset()
    sFsamp.reset()
    sFpass.reset()
    sFstop.reset()
    sW.reset()
    sFnotch1.reset()
    sFnotch2.reset()
button1.on_clicked(reset)
        
generateax = plt.axes([0.3, 0.025, 0.05, 0.04])
button2 = Button(generateax, 'Generate', color=axcolor, hovercolor='green')
def generate(event):
    global h_global
    print h_global[0:(len(h_global)+1)/2]
button2.on_clicked(generate)                        
                        
plt.show()
    
    
