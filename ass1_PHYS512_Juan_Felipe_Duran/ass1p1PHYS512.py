# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 12:21:15 2020

Assignment 1 PHYS 512
"""

'''
Problem 1
a) As seen in lecture 1, Our estimate was [f(x+dx)-f(x-dx)]/2dx 
since this cancels the second order term in the taylor expansion.
In order to cancel the third term, I used [f(x-2dx)-8*f(x-dx)+8*f(x+dx)-f(x+2*dx)]/(12*dx).
b)  With this new estimation, dx should be roughly (10^-16)^(1/4)~10^-4.
Below is a routine that shows for exp(x) and exp(0.01x) that dx works best at 10^-4.
'''


import numpy as np

expvals=np.linspace(-10,-2,17)
x0=0.01
xsmall = 0.01
truth=np.exp(x0)
truth2=np.exp(xsmall)
for myexp in expvals:
    dx=10**myexp
    f1=np.exp(x0+2*dx)
    f2=np.exp(x0+dx)
    f3=np.exp(x0-dx)   
    f4=np.exp(x0-2*dx)
    
    f1small=np.exp(xsmall+2*dx)
    f2small=np.exp(xsmall+dx)
    f3small=np.exp(xsmall-dx)   
    f4small=np.exp(xsmall-2*dx)
    
    deriv=(f4-8*f3+8*f2-f1)/(12*dx)  
    derivSmall=(f4small-8*f3small+8*f2small-f1small)/(12*dx)
    
    print(myexp,np.abs(deriv-truth2),np.abs(derivSmall-truth))
    

