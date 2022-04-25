# INVLAP.py
Python function INVLAP() used to numerically evaluate inverse Laplace transforms.

This Python function is used to numerically evaluate the inverse Laplace transform of F(s) to determine f(t).  The method used to find the inverse Laplace transfrom was developed by Juraj Valsa and Lubomir Brancik.  The method is described in a 1998 paper entitled "Approximate Formulae for Numerical Inversion of Laplace Transforms" in the International Journal of Numerical Modelling: Electronic Newtorks, Devices and Fields.  The journal reference is *Int. J. Numer. Model.* **11**, 153-166 (1998) - DOI: https://doi.org/10.1002/(SICI)1099-1204(199805/06)11:3%3C153::AID-JNM299%3E3.0.CO;2-C

Valsa and Brancik also coded a MATLAB implementation of their numerical inverse Laplace transform method called INVLAP.m.  The code was posted on the MathWorks File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/32824-numerical-inversion-of-laplace-transforms-in-matlab

INVLAP.py is a Python implentation of Valsa and Brancik's method for numerical inverse Laplace transforms.  It is essentially an translation of the MATLAB code containted in INVLAP.m.  The method can be used to quickly find the inverse Laplace transform f(t) of complex-valued functions F(s).

```python
#################################################################
## Example use of the INVLAP() function in INVLAP.py.          ##
#################################################################

# Calculate the transient response of a length of transmission line.  
# A voltage step is applied at one end of a trasmission line with 
# characteristic impedance Z0 = 50 ohms.  The voltage source has output
# impedance Rg = 1000 ohms.  The opposite end of the transmission line is
# open.  The voltage is measured at the node between Rg and the input of the
# transmission line.  For more details, see Am. J. Phys. 89, 96 (2021) - 
# DOI: https://doi.org/10.1119/10.0001896 or:
# https://arxiv.org/abs/2006.14381
# This problem has an analytic solution and does not require numerical
# methods.  However, we use it as an example to test the code.
#
# Import the required modules.  Calling the INVLAP function
# will be the most straightforward if one imports SymPy using 
# "from sympy import *".
#     
import matplotlib.pyplot as plt
import INVLAP
from sympy import *

# Define the symbol s to use in the symbolic expression for F(s). 
s = Symbol('s')

# Enter some parameter values.
Z0 = 50 # ohms
v0 = 0.6795*3e8 # m/s
ell = 6.32 # m
V0 = 1 # volts
Rg = 1000 # ohms
 
# Enter the symbolic expression for F(s).
Fs = (V0/s)*(1/((Rg/Z0)*tanh(s*ell/v0) + 1))
# Convert the symbolic expression to a string.
Fs = str(Fs) # F(s) as a string
# In this example, the F(s) expression does not have complex quantities with
# j or I = sqrt(-1).  However, if F(s) does have complex quantities, the string
# representation of the function will contain "I".  All of the I's need to be
# replaced by "1j" to work in the INVLAP function. 
Fs = Fs.replace("I", "1j")
# Call the INVAP function contained in INVLAP.py.
t, fcn = INVLAP.INVLAP(Fs, 1e-10, 0.25e-5, 10000, 6, 8000, 160)

# Plot the results.
plt.plot(t, fcn)
```
