

#
def generate_model_params(theta,eta,a1):

	'''
	Hola, soy un codigo escribido por "4b_document_fits.py". 
	~londonl@mit.edu/pilondon2@gmail.com 2020
	'''  

	# Import usefuls
	from numpy import cos, sqrt

	# Preliminaries
	u = cos(theta)
	u2 = u*u
	u3 = u2*u
	u4 = u3*u
	eta2 = eta*eta
	eta3 = eta2*eta
	delta = sqrt(1-4*eta)
	delta2 = delta*delta
	delta3 = delta2*delta

	# mu2
	mu2 = -3.41117957e+00*(a1) + -1.26114165e+01*(eta) + -1.78260476e+00*(u) + 1.14096945e+00 + 2.87119155e+00*(u*a1) + 1.74619307e+01*(u*eta) + 3.24528886e+01*(eta*eta) + 5.62963626e+01*(eta*a1) + -4.07011203e+01*(u*eta*eta) + -2.82964910e+01*(u*eta*a1) + -2.73998518e+02*(eta2*a1) + 3.10199443e+00*(u2*a1) + -2.35533179e+01*(u2*eta) + 4.45582933e+00*(u2*u) + 6.68689108e+01*(u*eta2*a1) + -3.90646042e+01*(u2*eta*a1) + 2.56411554e+02*(u2*eta*eta) + 4.15891046e+02*(eta3*a1) + -5.98335831e+01*(u3*eta) + -4.17750985e+00*(u3*a1) + 4.55511464e+01*(u3*eta*a1) + 2.53976464e+02*(u3*eta*eta) + -6.51185502e+02*(u2*eta2*eta) + 2.81747582e+01*(u4*eta) + -1.50564333e+00*(u4*a1) + 1.08936678e+02*(u2*eta2*a1) + -3.47186598e+02*(u3*eta2*eta) + -3.15151865e+02*(u4*eta*eta) + -1.14414257e+02*(u3*eta2*a1) + 8.15243101e+02*(u4*eta2*eta) + 1.46482764e+02*(u4*eta2*a1) + -5.08040843e+02*(u4*eta3*a1)

	# mu4
	mu4 = 2.55446068e-01*(eta) + -6.03579670e-02*(u) + 1.69491215e-01*(a1) + 8.56914198e-02*(delta) + -6.43252376e-02 + 4.21327055e-02*(u*delta) + 5.56652450e-02*(u*a1) + 7.75311113e-02*(u*u) + 2.52144757e-01*(u*eta) + -1.84521215e-01*(delta*a1) + -6.64903506e-01*(eta*a1) + -3.00520403e-01*(eta*delta) + -1.27802948e-01*(u*eta*delta) + -2.38747715e-01*(u*eta*a1) + -1.40526287e-01*(u2*a1) + -1.04990018e-01*(u2*delta) + -3.14899760e-01*(u2*eta) + 6.86535072e-01*(eta*delta*a1) + -1.18037962e-01*(u2*u) + 5.60294985e-01*(u2*eta*a1) + 3.86637885e-01*(u2*eta*delta) + 1.79806305e-01*(u2*delta*a1) + 4.50906259e-01*(u3*eta) + 1.52695734e-01*(u3*delta) + 2.36765772e-01*(u3*a1) + -4.32355250e-01*(u3*eta*delta) + -9.13839082e-01*(u3*eta*a1) + -3.14062854e-01*(u3*delta*a1) + -7.02878719e-01*(u2*eta*delta*a1) + 9.52370285e-01*(u3*eta*delta*a1)

	# nu4
	nu4 = -4.86966171e+00*(u) + 6.93863611e-01*(a1) + 4.37724131e-01*(delta) + 1.00119631e-01 + 4.38472579e+00*(u*delta) + 8.63522416e+00*(u*a1) + 3.86950384e-02*(u*u) + 1.94405999e+01*(u*eta) + -1.54686419e+00*(delta*a1) + -3.07883619e+00*(eta*a1) + -3.62179850e+00*(eta*delta) + -1.09865795e+01*(u*eta*delta) + -3.46312872e+01*(u*eta*a1) + -6.68201885e-01*(u2*delta) + 9.59431004e+00*(eta*delta*a1) + -7.82712377e+00*(u*delta*a1) + 2.22572127e+00*(u2*u) + 3.98679852e+00*(u2*eta*delta) + 1.08805705e+00*(u2*delta*a1) + 2.14012501e+01*(u*eta*delta*a1) + -8.83935683e+00*(u3*eta) + -1.19341771e+00*(u3*delta) + -7.84869328e+00*(u3*a1) + 3.14505468e+01*(u3*eta*a1) + 6.30404546e+00*(u3*delta*a1) + -8.47800646e+00*(u2*eta*delta*a1) + -1.25581055e+01*(u3*eta*delta*a1)

	# nu5
	nu5 = 4.47557999e-02 + 5.01608722e-02*(u) + -1.80499707e-01*(eta) + -9.64518529e-02*(a1) + -4.40717499e-02*(delta) + -7.12090858e-02*(u*delta) + -2.15294466e-01*(u*a1) + 3.96238219e-01*(eta*a1) + 1.23844728e-01*(eta*delta) + -2.02703922e-01*(u*eta) + 1.32977175e-01*(delta*a1) + 2.50059544e-01*(u*eta*delta) + -8.46678246e-03*(u2*delta) + -3.41052286e-01*(eta*delta*a1) + -1.76366685e-03*(u2*a1) + 2.81043848e-01*(u*delta*a1) + 8.63970308e-01*(u*eta*a1) + -8.58019225e-01*(u*eta*delta*a1) + -5.03147876e-03*(u2*delta*a1) + 3.51766382e-02*(u2*eta*delta)

	# nu6
	nu6 = -5.10099774e-02*(a1) + -3.41406518e-01*(eta) + 4.63607367e-02*(u) + 3.41214209e-02 + -8.69091039e-02*(u*a1) + -1.04941197e-02*(u*u) + -1.39915858e+00*(u*eta) + 8.31105133e-01*(eta*eta) + 7.24763360e-01*(eta*a1) + 1.00885861e+01*(u*eta*eta) + 2.52071866e+00*(u*eta*a1) + -2.74324393e+00*(eta2*a1) + 6.58501754e-02*(u2*eta) + 7.44862400e-02*(u2*a1) + -7.07972042e-02*(u2*u) + -1.76329799e+01*(u*eta2*a1) + -2.09338044e+01*(u*eta2*eta) + -1.64496232e+00*(u2*eta*a1) + 2.63280607e+00*(eta3*a1) + -6.60176328e-02*(u3*u) + 1.94387173e+00*(u3*eta) + 2.03636795e-01*(u3*a1) + -4.90063664e+00*(u3*eta*a1) + -1.37798881e+01*(u3*eta*eta) + 3.56989888e+01*(u*eta3*a1) + 1.22697396e+00*(u4*eta) + 5.74746293e-02*(u4*a1) + 9.92895162e+00*(u2*eta2*a1) + 2.85006523e+01*(u3*eta2*eta) + -6.91287151e+00*(u4*eta*eta) + -5.98405616e-01*(u4*eta*a1) + -1.84783675e+01*(u2*eta3*a1) + 3.29648998e+01*(u3*eta2*a1) + -6.63469132e+01*(u3*eta3*a1) + 1.18265916e+01*(u4*eta2*eta) + 1.56733960e+00*(u4*eta2*a1)

	# zeta2
	zeta2 = -1.50295387e+00*(a1) + -6.71407118e+00*(eta) + 1.73481823e+00*(u) + 8.99453688e-01 + 9.33172371e-01*(u*u) + -3.79524956e+01*(u*eta) + 1.16545488e+01*(eta*eta) + 2.27851348e+02*(u*eta*eta) + 4.90948016e+00*(u*eta*a1) + 7.64898916e+01*(eta2*a1) + -2.33729469e+01*(u2*eta) + -4.22451064e+00*(u2*a1) + -2.26857246e+00*(u2*u) + -5.38953014e+01*(u*eta2*a1) + -4.20343240e+02*(u*eta2*eta) + 8.19250662e+01*(u2*eta*a1) + 1.88737058e+02*(u2*eta*eta) + -2.15713487e+02*(eta3*a1) + -7.95691800e-01*(u3*u) + 7.26645929e+01*(u3*eta) + -5.59039359e+00*(u3*a1) + -5.83523803e+02*(u2*eta2*a1) + 5.33046400e+01*(u3*eta*a1) + -5.50410053e+02*(u3*eta*eta) + -4.16773170e+02*(u2*eta2*eta) + 1.39646048e+02*(u*eta3*a1) + 4.88682486e+00*(u4*eta) + 1.55689938e+00*(u4*a1) + 1.19417790e+03*(u3*eta2*eta) + -1.28317317e+01*(u4*eta*eta) + -6.11717393e+00*(u4*eta*a1) + 1.28888503e+03*(u2*eta3*a1) + -1.27035352e+02*(u3*eta2*a1)

	#
	return mu2,mu4,nu4,nu5,nu6,zeta2
