import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob
import math


plt.ion()

#####THIS IS A FUNCTION FOR CORRELATION
def xcorr(a,b):
    aft=np.fft.rfft(a)
    bft=np.fft.rfft(b)
    return np.fft.irfft(aft*np.conj(bft))


########tHIS IS A FUNCTION TO SMOOTH VECTORS WITH BOXCART
def smooth_vec(a,width):
    aft=np.fft.rfft(a)
    vec=np.zeros(len(a))
    vec[:width]=1
    vec[-width+1:]=1
    vec=vec/np.sum(vec)
    vecft=np.fft.rfft(vec)
    return np.fft.irfft(vecft*aft,len(a))



######READ template FUNCTIONS
def read_template(filename):
    dataFile=h5py.File(filename,'r')
    template=dataFile['template']
    th=template[0]
    tl=template[1]
    return th,tl
def read_file(filename):
    dataFile=h5py.File(filename,'r')
    dqInfo = dataFile['quality']['simple']
    qmask=dqInfo['DQmask'][...]

    meta=dataFile['meta']
    #gpsStart=meta['GPSstart'].value
    gpsStart=meta['GPSstart'][()]
    #print meta.keys()
    #utc=meta['UTCstart'].value
    utc=meta['UTCstart'][()]
    #duration=meta['Duration'].value
    duration=meta['Duration'][()]
    #strain=dataFile['strain']['Strain'].value
    strain=dataFile['strain']['Strain'][()]
    dt=(1.0*duration)/len(strain)

    dataFile.close()
    return strain,dt,utc



########I upload all files
strain1,dt1,utc1=read_file("H-H1_LOSC_4_V2-1126259446-32.hdf5")
strain2,dt2,utc2=read_file("L-L1_LOSC_4_V2-1126259446-32.hdf5")
strain3,dt3,utc3=read_file("H-H1_LOSC_4_V2-1128678884-32.hdf5")
strain4,dt4,utc4=read_file("L-L1_LOSC_4_V2-1128678884-32.hdf5")
strain5,dt5,utc5=read_file("H-H1_LOSC_4_V2-1135136334-32.hdf5")
strain6,dt6,utc6=read_file("L-L1_LOSC_4_V2-1135136334-32.hdf5")
strain7,dt7,utc7=read_file("H-H1_LOSC_4_V1-1167559920-32.hdf5")
strain8,dt8,utc8=read_file("L-L1_LOSC_4_V1-1167559920-32.hdf5")

th1,tl2=read_template("GW150914_4_template.hdf5")
th3,tl4=read_template("LVT151012_4_template.hdf5")
th5,tl6=read_template("GW151226_4_template.hdf5")
th7,tl8=read_template("GW170104_4_template.hdf5")


######################################################################
############################ Problem A)###############################
######################################################################
####This is a function for the  match filter

def matchedfilter(strain,template):
    x=np.linspace(-1,1,len(strain))*np.pi
    
    #For window, I used a cos, but I set it flat in the middle. 
    win=0.5+0.5*np.cos(x)
    i=0
    
    for e in win:
        
        if e>0.6:
            win[i]=0.6
        i+=1     

    
    strain_windowed=win*strain
    
    
    strain_ft=np.fft.rfft(strain_windowed)
    Aft=np.fft.rfft(win*template)
    
    dft=strain_ft
    #NOISE MODEL HERE
    #I use smooth vec to get a smoother version of power spectrum
    #I have np.maximum for cases when the noise is too underweighted
    N=np.abs(strain_ft)**2
    N2=np.abs(smooth_vec(N,5))
    N=np.maximum(N,N2)
    
    ######Calculation for sigma#####
    Ninv=1/N
    det=np.abs(xcorr(Ninv*Aft,Aft))
    
    sigmaFT= np.sqrt(1/det)
    sigma = np.fft.irfft(sigmaFT)
    
    ################################
    mf_ft=np.conj(Aft)*(dft/N)
    mf=np.fft.irfft(mf_ft)
    sigma = np.fft.irfft(sigmaFT, len(mf))
    return mf, sigma


mf1, sigma1 = matchedfilter(strain1,th1)
mf2, sigma2 = matchedfilter(strain2,tl2)
mf3, sigma3 = matchedfilter(strain3,th3)
mf4, sigma4 = matchedfilter(strain4,tl4)
mf5, sigma5 = matchedfilter(strain5,th5)
mf6, sigma6 = matchedfilter(strain6,tl6)
mf7, sigma7 = matchedfilter(strain7,th7)
mf8, sigma8 = matchedfilter(strain8,tl8)


print("Problem a): see noise model inside matchfilter")
print(" ")
######################################################################
############################ Problem B)###############################
######################################################################
fig, axs = plt.subplots(2)
fig.suptitle('1st event')
axs[0].plot(np.abs(mf1))
axs[1].plot(np.abs(mf2))
plt.savefig("event1.png")

fig, axs = plt.subplots(2)
fig.suptitle('2nd event')
axs[0].plot(np.abs(mf3))
axs[1].plot(np.abs(mf4))
plt.savefig("event2.png")

fig, axs = plt.subplots(2)
fig.suptitle('3rd event')
axs[0].plot(np.abs(mf5))
axs[1].plot(np.abs(mf6))
plt.savefig("event3.png")

fig, axs = plt.subplots(2)
fig.suptitle('4th event')
axs[0].plot(np.abs(mf7))
axs[1].plot(np.abs(mf8))
plt.savefig("event4.png")

print("problem b): see plots ")
print(" ")
######################################################################
############################ Problem C)###############################
######################################################################

SNR1=np.average(np.abs(mf1/sigma1))
SNR2=np.average(np.abs(mf2/sigma2))
SNR3=np.average(np.abs(mf3/sigma3))
SNR4=np.average(np.abs(mf4/sigma4))
SNR5=np.average(np.abs(mf5/sigma5))
SNR6=np.average(np.abs(mf6/sigma6))
SNR7=np.average(np.abs(mf7/sigma7))
SNR8=np.average(np.abs(mf8/sigma8))


SNRCombined1= np.sqrt(SNR1**2+SNR2**2)
SNRCombined2= np.sqrt(SNR3**2+SNR4**2)
SNRCombined3= np.sqrt(SNR5**2+SNR6**2)
SNRCombined4= np.sqrt(SNR7**2+SNR8**2)



print("Problem c): The 8 SNR Are")
print("SNR1: ", SNR1)
print("SNR2: ", SNR2)
print("SNR3: ", SNR3)
print("SNR4: ", SNR4)
print("SNR5: ", SNR5)
print("SNR6: ", SNR6)
print("SNR7: ", SNR7)
print("SNR8: ", SNR8)
print("The 4 combined SNRS are")
print("SNRCombined1:", SNRCombined1)
print("SNRCombined2:", SNRCombined2)
print("SNRCombined3:", SNRCombined3)
print("SNRCombined4:", SNRCombined4)
print(" ")
######################################################################
############################ Problem D)###############################
######################################################################
def analyticSNR(strain, template):
    #Noise model
    N=np.abs(strain)**2
    Ninv = 1/N
    A=template
    d=strain
    #Below, I am calculating m/sigma which is the SNR
    # Where m = A^T*N^-1*d/A^T*N^-1*A
    # sigma = 1/sqrt(A^T*N^-1*A)
    SNR = xcorr(Ninv*A,d)/np.sqrt(np.abs(xcorr(Ninv*A,A)))
    return SNR
    
ASNR1=np.average(np.abs(analyticSNR(strain1,th1)))
ASNR2=np.average(np.abs(analyticSNR(strain2,tl2)))
ASNR3=np.average(np.abs(analyticSNR(strain3,th3)))
ASNR4=np.average(np.abs(analyticSNR(strain4,tl4)))
ASNR5=np.average(np.abs(analyticSNR(strain5,th5)))
ASNR6=np.average(np.abs(analyticSNR(strain6,tl6)))
ASNR7=np.average(np.abs(analyticSNR(strain7,th7)))
ASNR8=np.average(np.abs(analyticSNR(strain8,tl8)))


print("Problem D): The 8 analytical SNRS are")
print("ASNR1: ", ASNR1)
print("ASNR2: ", ASNR2)
print("ASNR3: ", ASNR3)
print("ASNR4: ", ASNR4)
print("ASNR5: ", ASNR5)
print("ASNR6: ", ASNR6)
print("ASNR7: ", ASNR7)
print("ASNR8: ", ASNR8)
print(" ")



print("The SNR disagree, the analytical ones are bigger. The reason is probably that the non-analytica values have smaller noise")
print(" ")
######################################################################
############################ Problem E)###############################
######################################################################
def problemE(strain,template,event):
    x=np.linspace(-1,1,len(strain))*np.pi
    win=0.5+0.5*np.cos(x)
    i=0
    for e in win:       
       if e>0.6:
           win[i]=0.6
       i+=1     
       
    strain_windowed=win*strain   
    strain_ft=np.fft.rfft(strain_windowed)
    N=np.abs(strain_ft)**2
    N2=smooth_vec(N,5)
    N=np.maximum(N,N2)
    
    powerSpectrum= np.abs(np.fft.rfft(win*template)/np.sqrt(N))**2

    xfreq=np.fft.fftfreq(len(powerSpectrum),1/4096)
    y=np.cumsum(powerSpectrum)/np.sum(powerSpectrum)
    
    yy = [x - 0.5 for x in y]
    
    yy=np.abs(yy)
    
    index=np.argmin(yy)
    chosenFreq=xfreq[index]
    print("Chosen freq for event ", event," is ", chosenFreq)

   
    
    plt.plot(xfreq,y)

print("Problem E): ")
print(" ")

problemE(strain1,th1,1)
problemE(strain2,tl2,2)
problemE(strain3,th3,3)
problemE(strain4,tl4,4)
problemE(strain5,th5,5)
problemE(strain6,tl6,6)
problemE(strain7,th7,7)
problemE(strain8,tl8,8)
   





######################################################################
############################ Problem F)###############################
######################################################################
print("Problem F): please see section Problem F) for answer")
'''
To localize the time of arrival, I attempted an FWHM 
of the event, it gave me 1.76

'''

print("The estimation for FWHM gives me ", 2*np.sqrt(2*(math.log(2)))*np.std(mf1[1796:1815]))

b =np.where(mf2==np.amax(mf2))
a =np.where(mf1==np.amax(mf1))

#print("one number is", a," and the other is ", b)
'''
From this info, one arrives at 1790, other at 1768
This gives a difference of 22
with a step time of 4096/s, the time difference is
0.00537s
The uncertainty is half of the 1/(step time) so 0.0001s
From the equation 

L/c = D*sin(theta)/c

L/c = A = distance difference of arrival/speed of light =0.00537s
c = speed of light =299792km/s
D =distance between interferometers = 2000km
and small angle approximation sin(theta)=theta

We have A*c/D = theta

for uncertainty 0.0001
The uncertainty in theta is
0.001 *c/D = 0.015 roughly
'''





















#spec,nu=measure_ps(strain,do_win=True,dt=dt,osamp=16)
#strain_white=noise_filter(strain,numpy.sqrt(spec),nu,nu_max=1600.,taper=5000)

#th_white=noise_filter(th,numpy.sqrt(spec),nu,nu_max=1600.,taper=5000)
#tl_white=noise_filter(tl,numpy.sqrt(spec),nu,nu_max=1600.,taper=5000)


#matched_filt_h=numpy.fft.irfft(numpy.fft.rfft(strain_white)*numpy.conj(numpy.fft.rfft(th_white)))
#matched_filt_l=numpy.fft.irfft(numpy.fft.rfft(strain_white)*numpy.conj(numpy.fft.rfft(tl_white)))




#copied from bash from class
# strain2=np.append(strain,np.flipud(strain[1:-1]))
# tobs=len(strain)*dt
# k_true=np.arange(len(myft))*dnu