from math import *

tau = 0.25
fac = pow(tau/pi,1.5)
hi = ceil(sqrt(2.3*16/tau))


x = 0.5
y = 0.5
z = 0.5

s = 0.0 # verify norm of gaussian is computed as 1
u = 0.0 # verify sum of screened potential
for X in range (-hi,hi+1):
    for Y in range(-hi,hi+1):
        for Z in range(-hi,hi+1):
            rR = sqrt((x+X)**2 + (y+Y)**2 + (z+Z)**2)
            s += fac*exp(-tau*rR**2)
            u += erfc(sqrt(tau)*rR)/rR

print("tau", tau, "fac", fac, "hi", hi)
print("norm", s, 1.0)
print("potn", u, pi/tau)


            
            
