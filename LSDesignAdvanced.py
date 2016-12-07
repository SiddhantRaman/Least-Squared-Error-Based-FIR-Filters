import numpy as np
import LSFIR
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

Fpass2 = 11.0     # MHz Passband end frequency
Fstop2 = 15.0     #MHz Stopband start frequency
Fstop1 = 5.0
Fpass1 = 4.0
Fsamp = 50.0     # MHz Sampling Frequency
Weight = 100    # Weight of stop band error
Taps = 71      # FIR Filter Taps
Fnotch1 = 60   #Notch out freq Fnotch1
Fnotch2 = 60

h_global = np.zeros(Taps)
design = 0

fig, myFilter = plt.subplots(figsize=(18, 9))
plt.subplots_adjust(left=0.25,bottom=0.3)
myFilter.set_xlabel('Freq(MHz)')
myFilter.set_ylabel('Magnitude(dB)')
 
h = LSFIR.lpfls(Taps,2*np.pi*(Fpass2/Fsamp),2*np.pi*(Fstop2/Fsamp),Weight)
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
axTaps = plt.axes([0.6, 0.175, 0.25, 0.02], axisbg=axcolor)
axFsamp = plt.axes([0.6, 0.150, 0.25,0.02], axisbg=axcolor)
axFpass2 = plt.axes([0.6, 0.125, 0.25, 0.02],axisbg=axcolor)
axFstop2 = plt.axes([0.6, 0.100,0.25, 0.02],axisbg=axcolor)
axFstop1 = plt.axes([0.6, 0.075,0.25, 0.02],axisbg=axcolor)
axFpass1 = plt.axes([0.6, 0.050,0.25, 0.02],axisbg=axcolor)
axW = plt.axes([0.1, 0.175, 0.25, 0.02],axisbg=axcolor)   
axFnotch1 = plt.axes([0.1, 0.150, 0.25, 0.02],axisbg=axcolor)  
axFnotch2 = plt.axes([0.1, 0.125, 0.25, 0.02],axisbg=axcolor)

sTaps = Slider(axTaps,'Taps',15,200,valinit=Taps)
sFsamp = Slider(axFsamp,'Samp Freq',1.0,100,valinit=Fsamp)
sFpass1 = Slider(axFpass1,'PassFreq2',1.0,40.0,valinit=Fpass1)
sFstop1 = Slider(axFstop1,'StopFreq2',1.0,45.0,valinit=Fstop1)
sFstop2 = Slider(axFstop2, 'StopFreq1',1.0,45.0,valinit=Fstop2)
sFpass2 = Slider(axFpass2,'PassFreq1',1.0,45.0,valinit=Fpass2)
sW = Slider(axW,'Weight',1.0,1000,valinit=Weight)
sFnotch1 = Slider(axFnotch1,'Notch 1',0.0,60,valinit=Fnotch1)
sFnotch2 = Slider(axFnotch2,'Notch 2',0.0,60,valinit=Fnotch2)

def update(val):
    global design
    Taps = int(sTaps.val)
    if(Taps%2 == 0) :
        Taps = Taps +1 
    Fsamp = sFsamp.val
    Fpass1 = sFpass1.val
    Fstop1 = sFstop1.val
    Fstop2 = sFstop2.val
    Fpass2 = sFpass2.val
    W = sW.val
    Fnotch1 = sFnotch1.val
    Fnotch2 = sFnotch2.val
    
    wp = 2*np.pi*(Fpass2/Fsamp)
    ws = 2*np.pi*(Fstop2/Fsamp)
    wp2 = 2*np.pi*(Fpass1/Fsamp)
    ws2 = 2*np.pi*(Fstop1/Fsamp)
    wn1 = 2*np.pi*(Fnotch1/Fsamp)
    wn2 = 2*np.pi*(Fnotch2/Fsamp)

    if(design == 0):
        if((Fnotch1 <= Fsamp/2) and (Fnotch2 <= Fsamp/2)) :
            h = LSFIR.lpfls2notch(Taps,wp,ws,wn1,wn2,W)
        elif((Fnotch1 > Fsamp/2) and (Fnotch2 <= Fsamp/2)):
            h = LSFIR.lpfls1notch(Taps,wp,ws,wn2,W)
        elif((Fnotch1 <= Fsamp/2) and (Fnotch2 > Fsamp/2)):
            h = LSFIR.lpfls1notch(Taps,wp,ws,wn1,W)
        else:
            h = LSFIR.lpfls(Taps,wp,ws,W)
    elif(design == 1):
        h = LSFIR.bpfls(Taps,ws2,wp2,wp,ws,W)
    elif(design == 2):
        h = LSFIR.hpfls(Taps, ws2, wp2, W)
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
sFpass1.on_changed(update)
sFstop1.on_changed(update)
sFpass2.on_changed(update)
sFstop2.on_changed(update)
sW.on_changed(update)
sFnotch1.on_changed(update)
sFnotch2.on_changed(update)

    
resetax = plt.axes([0.1, 0.025, 0.05, 0.04])
button1 = Button(resetax, 'Reset', color=axcolor, hovercolor='green')
def reset(event):
    sTaps.reset()
    sFsamp.reset()
    sFpass1.reset()
    sFstop1.reset()
    sFstop2.reset()
    sFpass2.reset()
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

lowpassax = plt.axes([0.025, 0.5, 0.05, 0.04])
button3 = Button(lowpassax, 'LowPass', color=axcolor, hovercolor='green')
def LowPassDesignSet(event):
    global design
    design = 0

button3.on_clicked(LowPassDesignSet)
button3.on_clicked(reset)
button3.on_clicked(update)

highpassax = plt.axes([0.025, 0.4, 0.05, 0.04])
button4 = Button(highpassax, 'HighPass', color=axcolor, hovercolor='green')
def HighPassDesignSet(event):
    global design
    design = 2

button4.on_clicked(HighPassDesignSet)
button4.on_clicked(reset)
button4.on_clicked(update)

bandpassax = plt.axes([0.025, 0.3, 0.05, 0.04])
button5 = Button(bandpassax, 'BandPass', color=axcolor, hovercolor='green')
def BandPassDesignSet(event):
    global design
    design = 1

button5.on_clicked(BandPassDesignSet)
button5.on_clicked(reset)
button5.on_clicked(update)
     
                        
plt.show()
    
    
