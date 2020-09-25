

'''
Assignment 2, problem 1:
    
I didn't add any comments between lines 83 and 107 because that's just the starter 
code of the variable step size available in the Github.
When running this file, it will output the neval of the class function and my 
function for the three different functions.
My function is between lines 28 and 75

'''
import numpy as np


sig=0.1

def fun(x):
    return 1.0/(1.0+x**2)

def fun2(x):
    return 1.0+np.exp(-0.5*x**2/(sig**2))
'''
For this function, I added three parameters: The values of the left and right
edge, and the side; whether its right or left. 

'''
def simple_integrate(fun,a,b,tol,leftEdge,rightEdge,side):
    
    if side==None:
        x=np.linspace(a,b,5)
    elif side=="left":
        #Here I only create a linspace with 4 items to later recycle the first evaluation
        x=np.linspace(a+((b-a)/4),b,4)
    elif side=="right":
        #Here I only create a linspace with 4 items to later recycle the last evaluation
        x=np.linspace(a,b-((b-a)/4),4)
    
    y=fun(x)
    neval=len(x) 
    
    if side==None:
        #Like in class
        f1=(b-a)*(y[0]+4*y[2]+y[4])/6.0
        f2=(b-a)*( y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/12.0
        leftEdge=y[0]
        rightEdge=y[4]
    elif side=="left":
        #Here I recycle the firt evaluation
        f1=(b-a)*(leftEdge+4*y[1]+y[3])/6.0
        f2=(b-a)*(leftEdge+4*y[0]+2*y[1]+4*y[2]+y[3])/12.0
        leftEdge=leftEdge
        rightEdge=y[3]
    elif side=="right":
        #Here I recycle the last evaluation
        f1=(b-a)*(y[0]+4*y[2]+rightEdge)/6.0
        f2=(b-a)*(y[0]+4*y[1]+2*y[2]+4*y[3]+rightEdge)/12.0
        leftEdge=y[0]
        rightEdge=rightEdge

    myerr=np.abs(f2-f1)
    #print([a,b,f1,f2])
    if (myerr<tol):
        #return (f2)/1.0,myerr,neval
        return (16.0*f2-f1)/15.0,myerr,neval
    else:
        mid=0.5*(b+a)
        #Calling the function with parameter "left"
        f_left,err_left,neval_left=simple_integrate(fun,a,mid,tol/2.0,leftEdge,rightEdge,"left")
        #Calling the function with parameter "right"
        f_right,err_right,neval_right=simple_integrate(fun,mid,b,tol/2.0,leftEdge,rightEdge,"right")
        neval=neval+neval_left+neval_right
        f=f_left+f_right
        err=err_left+err_right
        return f,err,neval
    
f,err,neval1=simple_integrate(np.exp,-1,1,10e-4,0,0,None);pred=np.exp(1)-np.exp(-1)
f,err,neval2=simple_integrate(fun,-1,1,1e-4,0,0,None);pred=np.arctan(1)-np.arctan(-1)
a=-5;b=5;f,err,neval3=simple_integrate(fun2,a,b,1e-4,0,0,None);pred=(b-a)+np.sqrt(2*np.pi)*sig



def simple_integrate(fun,a,b,tol):
    x=np.linspace(a,b,5)
    
    #np.median(np.diff(x))
    y=fun(x)
    neval=len(x) #let's keep track of function evaluations
    f1=(y[0]+4*y[2]+y[4])/6.0*(b-a)
    f2=(y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/12.0*(b-a)
    myerr=np.abs(f2-f1)
    
    if (myerr<tol):
        #return (f2)/1.0,myerr,neval
        return (16.0*f2-f1)/15.0,myerr,neval
    else:
        mid=0.5*(b+a)
        f_left,err_left,neval_left=simple_integrate(fun,a,mid,tol/2.0)
        f_right,err_right,neval_right=simple_integrate(fun,mid,b,tol/2.0)
        neval=neval+neval_left+neval_right
        f=f_left+f_right
        err=err_left+err_right
        return f,err,neval
        
f,err,neval4=simple_integrate(np.exp,-1,1,10e-4);pred=np.exp(1)-np.exp(-1)
f,err,neval5=simple_integrate(fun,-1,1,1e-4);pred=np.arctan(1)-np.arctan(-1)
a=-5;b=5;f,err,neval6=simple_integrate(fun2,a,b,1e-4);pred=(b-a)+np.sqrt(2*np.pi)*sig


#Here I print the nevals for each function for each code
print("function1:")
print("In class neval",neval1, "home neval",neval4)
print("function2:")
print("In class neval",neval2, "home neval",neval5)
print("function3:")
print("In class neval",neval3, "home neval",neval6)


