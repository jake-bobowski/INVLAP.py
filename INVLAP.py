# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 19:18:36 2022

@author: Jake Bobowski
"""

# This Python function is used to numerically evaluate the inverse Laplace 
# transform of F(s) to determine f(t).  The method used to find the inverse
# Laplace transfrom was developed by Juraj Valsa and Lubomir Brancik.  The
# method is described in a 1998 paper entitled "Approximate Formulae for
# Numerical Inversion of Laplace Transforms" in the International Journal of
# Numerical Modelling: Electronic Newtorks, Devices and Fields.  The journal
# reference is Int. J. Numer. Model. 11, 153-166 (1998) - 
# DOI: https://doi.org/10.1002/(SICI)1099-1204(199805/06)11:3%3C153::AID-JNM299%3E3.0.CO;2-C

# Valsa and Brancik also coded a MATLAB implementation of their 
# numerical inverse Laplace transform method called INVLAP.m.  The code
# was posted on the MathWorks File Exchange:
# https://www.mathworks.com/matlabcentral/fileexchange/32824-numerical-inversion-of-laplace-transforms-in-matlab

# INVLAP.py is a Python implentation of Valsa and Brancik's method for
# numerical inverse Laplace transforms.  It is essentially an translation
# of the MATLAB code containted in INVLAP.m.  The method can be used to quickly 
# find the inverse Laplace transform f(t) of complex-valued functions F(s).

#############################################################################
# Here is the orignal INVLAP.m wirtten by Valsa and Brancik.                #
#############################################################################

# % INVLAP – Numerical Inversion of Laplace Transforms 
# function [radt,ft]=INVLAP(Fs,tini,tend,nnt,a,ns,nd);
# % Fs is formula for F(s) as a string
# % tini, tend are limits of the solution interval
# % nnt is total number of time instants
# % a, ns, nd are parameters of the method 
# % if not given, the method uses implicit values a=6, ns=20, nd=19
# % it is recommended to preserve a=6
# % increasing ns and nd leads to lower error
# % an example of function calling  
# % [t,ft]=INVLAP('s/(s^2+4*pi^2)',0,10,1001);
# % to plot the graph of results write plot(t,ft), grid on, zoom on
# FF=strrep(strrep(strrep(Fs,'*','.*'),'/','./'),'^','.^');
# if nargin==4
#   a=6; ns=20; nd=19;  end;    % implicit parameters
# radt=linspace(tini,tend,nnt); % time vector
# if tini==0  radt=radt(2:1:nnt); nnt=nnt-1;  end;  % t=0 is not allowed
# % Note that Jake Bobowski added the nnt=nnt-1 in the line above.
# % nnt needs to be decremented to make INVLAP.m work when tini=0 is used.
# tic					% measure the CPU time
# for n=1:ns+1+nd               % prepare necessary coefficients
#    alfa(n)=a+(n-1)*pi*j;
#    beta(n)=-exp(a)*(-1)^n;
# end;
# n=1:nd;
# bdif=fliplr(cumsum(gamma(nd+1)./gamma(nd+2-n)./gamma(n)))./2^nd;
# beta(ns+2:ns+1+nd)=beta(ns+2:ns+1+nd).*bdif;
# beta(1)=beta(1)/2;
# for kt=1:nnt                  % cycle for time t
#    tt=radt(kt);
#    s=alfa/tt;                 % complex frequency s
#    bt=beta/tt;
#    btF=bt.*eval(FF);          % functional value F(s)
#    ft(kt)=sum(real(btF));     % original f(tt)
# end;
# toc


#############################################################################
# Here is the Python implementation of INVLAP.                              #
#############################################################################

# When calling the INVLAP.py function, a string representation of F(s) is
# passed as one of the arguments.  The string expression is generated from
# a SymPy expression.  In the actaul INVLAP function, the string F(s)
# expression is then evaluted as a NumPy expression.  These conversions will
# be simplest if the script calling the INVLAP function uses:
#   from scipy import *
# and the function itself uses:
#   from numpy import *
#from numpy import *
#import math
#import time

# Here's the start of the INVLAP function.  The first 4 arguments are required
# The last three have default values set when alternative values are not
# specified by the user.
#def INVLAP(Fs, tini, tend, nnt, a = 6, ns = 20, nd = 19):
#    radt = linspace(tini, tend, nnt)
# An initial time tini - 0 is not allowed.  If the user specifies tini = 0,
# make a minor adjustment so that the function will still execute without
# an error.    
#    if tini == 0:
#        radt = radt[1:]
#        nnt = nnt - 1
#    t = time.time()
#    alfa = []
#    beta = []
#    for n in range(ns + nd + 1):
#        alfa = alfa + [a + n*pi*1j]
#        beta = beta + [-exp(a)*(-1)**(n + 1)]
#    beta = array(beta)
#    bdifSub = []
#    for n in arange(1, nd + 1):
#       bdifSub = bdifSub + [math.gamma(nd + 1)/math.gamma(nd + 2 - n)/math.gamma(n)]
#    bdif = flip(cumsum(bdifSub))/2**nd
#    beta[ns + 1:ns + nd + 1] = beta[ns + 1:ns + nd + 1]*bdif
#    beta[0]=beta[0]/2
#    ft = []
#    for kt in range(nnt):
#        tt = radt[kt]
#        s = alfa/tt
#        bt = beta/tt
#        btF = bt*eval(Fs) 
#        ft = ft + [sum(real(btF))]
#    elapsed = time.time() - t
# Print the elapsed time.
#    print("Elasped time", "{:.3f}".format(elapsed), "seconds.")
#    return radt, ft

## Example use of the INVLAP() function in INVLAP.py.  
## Calculate the transient response of a length of transmission line.  
## A voltage step is applied at one end of a trasmission line with 
## characteristic impedance Z0 = 50 ohms.  The voltage source has output
## impedance Rg = 1000 ohms.  The opposite end of the transmission line is
## open.  The voltage is measured at the node between Rg and the input of the
## transmission line.  For more details, see Am. J. Phys. 89, 96 (2021) - 
## DOI: https://doi.org/10.1119/10.0001896 or:
## https://arxiv.org/abs/2006.14381
## This problem has an analytic solution and does not require numerical
## methods.  However, we use it as an example to test the code.
##
## Import the required modules.  Remember, that calling the INVLAP function
## will be the most straightforward if one imports SymPy using 
## "from sympy import *".
##     
# import matplotlib.pyplot as plt
# import INVLAP
# from sympy import *
#
## Define the symbol s to use in the symbolic expression for F(s). 
# s = Symbol('s')
#
## Enter some parameter values.
# Z0 = 50 # ohms
# v0 = 0.6795*3e8 # m/s
# ell = 6.32 # m
# V0 = 1 # volts
# Rg = 1000 # ohms
# 
## Enter the symbolic expression for F(s).
# Fs = (V0/s)*(1/((Rg/Z0)*tanh(s*ell/v0) + 1))
## Convert the symbolic expression to a string.
# Fs = str(Fs) # F(s) as a string
## In this example, the F(s) expression does not have complex quantities with
## j or I = sqrt(-1).  However, if F(s) does have complex quantities, the string
## representation of the function will contain "I".  All of the I's need to be
## replaced by "1j" to work in the INVLAP function. 
#Fs = Fs.replace("I", "1j")
## Call the INVAP function contained in INVLAP.py.
#t, fcn = INVLAP.INVLAP(Fs, 1e-10, 0.25e-5, 10000, 6, 8000, 160)
#
## Plot the results.
#plt.plot(t, fcn)










## An optimized version of the INVLAP function written with the assistance of
## CHATgpt and the original INVLAP script written by Juraj Valsa and Lubomir
## Brancik.  
## February 19, 2025  


from numpy import *
import time
import scipy.special as sp

def INVLAP(Fs, tini, tend, nnt, a = 6, ns = 20, nd = 19):
    # Generate time array
    radt = linspace(tini, tend, nnt)
    
    # If tini == 0, adjust the time vector and nnt
    if tini == 0:
        radt = radt[1:]
        nnt -= 1

    # Start timing
    t = time.time()

    # Precompute alpha and beta arrays
    alfa = a + (arange(ns + nd + 1) * pi * 1j)
    beta = (-exp(a) * (-1) ** (arange(ns + nd + 1) + 1))

    # Adjust beta values for the additional terms
    bdifSub = sp.gamma(nd + 1) / (sp.gamma(nd + 2 - arange(1, nd + 1)) * sp.gamma(arange(1, nd + 1)))
    bdif = flip(cumsum(bdifSub)) / 2 ** nd
    beta[ns + 1:ns + nd + 1] *= bdif
    beta[0] /= 2

    # Initialize the ft array
    ft = zeros(nnt)

    # Loop over each time step instead of using broadcasting
    for kt in range(nnt):
        tt = radt[kt]
        
        # Broadcast for each time step and calculate btF for each time step
        s = alfa / tt
        bt = beta / tt
        btF = bt * Fs(s)
        
        # Sum and store the result
        ft[kt] = real(sum(btF))
    
    # Elapsed time
    elapsed = time.time() - t
    print(f"Elapsed time: {elapsed:.3f} seconds.")
    
    return radt, ft

def F(s, Z0, v0, ell, Rg, V0):
    # Lower the clamping threshold further
    max_value = 300  # Reduced max value for better control
    s_clamped = where(abs(s * ell / v0) > max_value, sign(s) * max_value, s * ell / v0)
    return (V0 / s) * (1 / ((Rg / Z0) * tanh(s_clamped) + 1))




# # Example use of the INVLAP() function in INVLAP.py.
# # Calculate the transient response of a length of transmission line.
# # A voltage step is applied at one end of a trasmission line with
# # characteristic impedance Z0 = 50 ohms.  The voltage source has output
# # impedance Rg = 1000 ohms.  The opposite end of the transmission line is
# # open.  The voltage is measured at the node between Rg and the input of the
# # transmission line. For more details, see Am. J. Phys. 89, 96 (2021) -
# # DOI: https://doi.org/10.1119/10.0001896 or:
# # https://arxiv.org/abs/2006.14381

# # Import the function
# import INVLAP
# import matplotlib.pyplot as plt

# # Define your parameters
# Z0 = 50
# v0 = 0.6795 * 3e8
# ell = 6.32
# V0 = 1
# Rg = 1000

# # Use F directly in the INVLAP function
# Fs = lambda s: INVLAP.F(s, Z0, v0, ell, Rg, V0)

# # Reduce nnt for testing to see if it helps with performance
# nnt = 1000  # Start with a smaller number of time steps
# t, fcn = INVLAP.INVLAP(Fs, 1e-10, 0.25e-5, nnt, 6, 8000, 160)

# plt.plot(t, fcn)
# plt.show()