import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm


#function to use
def R(x,y):
    return np.sqrt(x**2 + y**2)

def potential(x,y):

   
    
    r = R(x,y)
    V = np.log(r)
    return V

def green(x, y, ):
    r = R(x,y)
    r[r<0.1] = 0.1

    g = 1/(4*np.pi*r)
    return g

def convolution(f, g):
    F = np.fft.fft(f)
    G = np.fft.fft(g)
    ans = np.fft.ifft(F*G)/len(f)

    return ans

################################################################################
################################################################################

#PART.A)

#
#Setting up a 2d plane
N = 100

l = np.linspace(-N/2, N/2, N)

x, y = np.meshgrid(l, l)


V = potential(x, y)


#Calculate v[o,o] assuming v[0,1] is avergae 
o = int(N/2)

V[o, o] = 8*V[o+1, o] - V[o+2, o] - V[o+1, o+1] - V[o+1, o-1]-V[o+2,o-1]-V[o+2,o+1]+V[o,o-1]+V[o,o+1]

#set rho
rho = V - 0.25*(np.roll(V,1,axis=0) + np.roll(V,-1,axis=0) + np.roll(V,1,axis=1) + np.roll(V,-1,axis=1))



#Rescale to set rho[0,0]=1
shif=1-rho[o,o]
rho=rho+shif

#Shift V to set V[0,0] = 1
shift = 1 - V[o, o]
V = V + shift


#V=V*rho
print(V[o,o])

print('V[1,0] = ', V[1,0])
print('V[2,0] = ', V[2,0])
print('V[5,0] = ', V[5,0])

################################################################################
################################################################################

#PART B)
green = green(x,y)


V=convolution(green,rho)

print("V from convolution gives center", V[o,o])
