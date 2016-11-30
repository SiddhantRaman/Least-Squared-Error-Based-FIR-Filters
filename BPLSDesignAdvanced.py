import numpy as np
import LSFIR
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

Fpass1 = 11.0     # MHz Passband end frequency
Fstop1 = 5.0      #MHz Stopband start frequency
Fpass2 = 30.0
Fstop2 = 35.0
Fsamp = 90.0      # MHz Sampling Frequency
Weight = 100      # Weight of stop band error
Taps = 71         # FIR Filter Taps
Fnotch1 = 60      #Notch out freq Fnotch1
Fnotch2 = 60

h_global = np.zeros(Taps)

fig, myFilter = plt.subplots()
plt.subplots_adjust(left=0.25,bottom=0.25)
myFilter.set_xlabel('Freq(MHz)')
myFilter.set_ylabel('Magnitude(dB)')
 
h = LSFIR.bpfls(Taps,2*np.pi*(Fstop1/Fsamp),2*np.pi*(Fpass1/Fsamp), 2*np.pi*(Fpass2/Fsamp),2*np.pi*(Fstop2/Fsamp), Weight)
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
plt.axis([0, 50.0, -50, 100])
l2, = plt.plot(freq,Mag_quant, lw=2, color='red')

axcolor = 'lightgoldenrodyellow'
axTaps = plt.axes([0.6, 0.125, 0.25, 0.02], axisbg=axcolor)
axFsamp = plt.axes([0.6, 0.1, 0.25,0.02], axisbg=axcolor)
axFpass2 = plt.axes([0.6, 0.075, 0.25, 0.02],axisbg=axcolor)
axFstop2 = plt.axes([0.6, 0.05,0.25, 0.02],axisbg=axcolor)
axW = plt.axes([0.1, 0.125, 0.25, 0.02],axisbg=axcolor)   
axFpass1 = plt.axes([0.1, 0.1, 0.25, 0.02],axisbg=axcolor)  
axFstop1 = plt.axes([0.1, 0.075, 0.25, 0.02],axisbg=axcolor)

sTaps = Slider(axTaps,'Taps',15,200,valinit=Taps)
sFsamp = Slider(axFsamp,'Samp Freq',1.0,100,valinit=Fsamp)
sFpass2 = Slider(axFpass2,'PassEnd',1.0,40.0,valinit=Fpass2)
sFstop2 = Slider(axFstop2,'StopEnd',1.0,45.0,valinit=Fstop2)
sW = Slider(axW,'Weight',1.0,1000,valinit=Weight)
sFpass1 = Slider(axFpass1,'PassBegin',0.0,60,valinit=Fpass1)
sFstop1 = Slider(axFstop1,'StopBegin',0.0,60,valinit=Fstop1)

def update(val):
    Taps = int(sTaps.val)
    if(Taps%2 == 0) :
        Taps = Taps +1 
    Fsamp = sFsamp.val
    Fpass2 = sFpass2.val
    Fstop2 = sFstop2.val
    W = sW.val
    Fpass11 = sFpass1.val
    Fstop1 = sFstop1.val
    
    wp2 = 2*np.pi*(Fpass2/Fsamp)
    ws2 = 2*np.pi*(Fstop2/Fsamp)
    wp1 = 2*np.pi*(Fpass1/Fsamp)
    ws1 = 2*np.pi*(Fstop1/Fsamp)
    #h = np.zeros(0,Taps)
    h = LSFIR.bpfls(Taps,ws1,wp1,wp2,ws2,W)
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
sFpass2.on_changed(update)
sFstop2.on_changed(update)
sW.on_changed(update)
sFpass1.on_changed(update)
sFstop1.on_changed(update)

    
resetax = plt.axes([0.1, 0.025, 0.05, 0.04])
button1 = Button(resetax, 'Reset', color=axcolor, hovercolor='green')
def reset(event):
    sTaps.reset()
    sFsamp.reset()
    sFpass1.reset()
    sFstop1.reset()
    sW.reset()
    sFpass2.reset()
    sFstop2.reset()
button1.on_clicked(reset)
        
generateax = plt.axes([0.3, 0.025, 0.05, 0.04])
button2 = Button(generateax, 'Generate', color=axcolor, hovercolor='green')
def generate(event):
    global h_global
    print h_global[0:(len(h_global)+1)/2]
button2.on_clicked(generate)                        
                        
plt.show()
    
    
