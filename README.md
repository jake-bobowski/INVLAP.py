# INVLAP.py
Python function INVLAP() used to numerically evaluate inverse Laplace transforms.

This Python function is used to numerically evaluate the inverse Laplace transform of F(s) to determine f(t).  The method used to find the inverse Laplace transfrom was developed by Juraj Valsa and Lubomir Brancik.  The method is described in a 1998 paper entitled "Approximate Formulae for Numerical Inversion of Laplace Transforms" in the International Journal of Numerical Modelling: Electronic Newtorks, Devices and Fields.  The journal reference is *Int. J. Numer. Model.* **11**, 153-166 (1998) - DOI: https://doi.org/10.1002/(SICI)1099-1204(199805/06)11:3%3C153::AID-JNM299%3E3.0.CO;2-C

Valsa and Brancik also coded a MATLAB implementation of their numerical inverse Laplace transform method called INVLAP.m.  The code was posted on the MathWorks File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/32824-numerical-inversion-of-laplace-transforms-in-matlab

INVLAP.py is a Python implentation of Valsa and Brancik's method for numerical inverse Laplace transforms.  It is essentially an translation of the MATLAB code containted in INVLAP.m.  The method can be used to quickly find the inverse Laplace transform f(t) of complex-valued functions F(s).
