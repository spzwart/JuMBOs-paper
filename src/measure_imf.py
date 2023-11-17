from amuse.lab import *
import numpy

N = 10000
Njumbos = N*600/2000
print("aim at n=", Njumbos)

x = -2.0
x = -1.5
x = -0.99
mmin = 1 | units.MJupiter
dmmin = 0.5 |units.MJupiter
while True:
    print(mmin.value_in(units.MJupiter))
    m = new_salpeter_mass_distribution(N,
                                       mmin,
                                       15|units.MJupiter, alpha=x)

    m = numpy.sort(m.value_in(units.MJupiter))
    n = len(m[m>1])
    print(mmin, n)
    if n==Njumbos:
        break
    elif n<Njumbos:
        mmin += dmmin
    else:
        mmin -= dmmin
    dmmin = dmmin/2
print(mmin)
mmin = 0.3|units.MJupiter
print("mmin=", mmin)

m = new_salpeter_mass_distribution(N,
                                   mmin,
                                   15|units.MJupiter, alpha=x)
m = numpy.sort(m.value_in(units.MJupiter))
n = len(m[m>1])
print("Nstars=", len(m), "Nj=", n)
