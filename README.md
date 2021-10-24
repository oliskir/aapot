# pyaa
Computes the alpha-alpha scattering length for a given parametrisation of the nuclear potential.

Solves the 1-dimensional Schrodinger equation by placing the system in a large box, as described in the [scattering primer](docs/scattering_primer.pdf).

In addition to computing the alpha-alpha scattering length, the [pyaa.py](pyaa.py) script computes the capture rate to a Gaussian final state. 

Running the script should produce the below figure.

![Alpha-alpha scattering](output.png)


## Notes

Raquel's EPJA 2007 paper mentions test performed 
with three different a-a potentials, which all 
reproduce the experimental phase shifts, but result 
in somewhat different resonance energies. For the 
Be(2+) resonance,

 * E_R=3.1MeV, G=1.5MeV ("AB(a')")
 * E_R=3.2MeV, G=1.4MeV ("attr")
 * E_R=3.0MeV, G=1.4MeV ("AB(a'd')")

Note that the energies and widths differ by as 
much as 100-200 keV. Also, values are only quoted 
to a precision of 100 keV indicating that Raquel et al.
where not too concerned about the precise values of 
E_R and G (which are actually important for determining 
the Dalitz distributions of the 3a breakup).

It is also worth noting that the calculations were 
performed with a Gaussian three-body potential of the 
form S * exp(-rho^2 / b^2) with b fixed to 6 fm (corresponding 
to three touching alpha particles) and S adjusted to 
reproduce the observed 12C resonance energies. However, 
the effect of the choice of b on the Dalitz distributions
was not investigated. 

It is interesting to note that according to Bhattacharya 
and Adelberger 2002, the R-matrix energy and width of the 
Be(2+) resonance are constrained to a precision of 23keV 
and 37keV respectively, by the experimental phase shift data.

Therefore, it would appear that the R-matrix energy and width 
are more strongly constrained by the data than the resonance 
energy and width inferred from the type of a-a potentials 
considered by Ali and Bodmer.

It is also worth nothing that Ali and Bodmer use the same 
attractive potential for all three partial saves (s,d,f).
As far as I can tell, this restrictive choice is not 
dictated by laws of physics. Relaxing this assumption, 
one could likely produce a greater variety of a-a potentials,
all consistent with experimental phase shifts, but leading 
to greater variability in resonance energies and widths.
