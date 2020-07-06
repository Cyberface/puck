

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

	# mu0
	mu0 = 2.03787682e-01*(eta) + 2.68064186e+00*(u) + 6.90853785e-01*(a1) + 4.04392602e-02 + -2.96186870e+00*(u*delta) + -3.75508723e+00*(u*a1) + -1.80353662e+00*(u*u) + -1.06886186e+01*(u*eta) + -6.63735528e-01*(delta*a1) + -3.00485104e+00*(eta*a1) + -5.63085183e-01*(eta*delta) + 8.96332041e+00*(u*eta*delta) + 1.49635292e+01*(u*eta*a1) + -2.73696727e-01*(u2*a1) + 2.00737577e+00*(u2*delta) + 7.17724931e+00*(u2*eta) + 2.48467528e+00*(eta*delta*a1) + 4.04394461e+00*(u*delta*a1) + -3.72617894e+00*(u2*u) + 1.08706733e+00*(u2*eta*a1) + -5.22585484e+00*(u2*eta*delta) + 1.89547088e+00*(u3*u) + -1.17653381e+01*(u*eta*delta*a1) + 1.47446652e+01*(u3*eta) + 4.07076883e+00*(u3*delta) + 5.68482936e+00*(u3*a1) + -1.16431174e+01*(u3*eta*delta) + -2.25513943e+01*(u3*eta*a1) + -6.28954648e+00*(u3*delta*a1) + -6.99226697e+00*(u4*eta) + -2.98812465e-01*(u2*eta*delta*a1) + -1.24125140e-01*(u4*a1) + -2.01226022e+00*(u4*delta) + 1.85039121e+01*(u3*eta*delta*a1) + 4.15678961e+00*(u4*eta*delta) + 4.32435969e-01*(u4*delta*a1)

	# mu1
	mu1 = 1.47061471e-03*(u) + 7.27715292e-04*(a1) + 5.55728573e-03*(delta) + -1.07687956e-03 + -2.02191930e-03*(u*a1) + -1.56246533e-03*(u*u) + -1.69464109e-02*(delta*delta) + -4.34855329e-03*(delta*a1) + 1.91286382e-02*(delta2*a1) + 3.66669275e-03*(u2*a1) + 1.16429059e-02*(u2*delta) + -3.95171151e-03*(u2*u) + -1.16708470e-02*(u*delta*a1) + 2.28911853e-02*(delta2*delta) + -3.70474359e-02*(delta3*a1) + 5.05327551e-02*(u2*delta*delta) + -3.93604217e-02*(u2*delta*a1) + -2.40471020e-03*(u3*u) + 1.57444749e-02*(u3*delta) + 6.17003861e-03*(u3*a1) + 3.44109893e-02*(u*delta2*a1) + -2.78806477e-02*(u*delta3*a1) + -5.62655353e-02*(u3*delta*a1) + -9.45418255e-02*(u2*delta2*delta) + 8.78144838e-02*(u2*delta3*a1) + -1.23458917e-02*(u4*delta*a1) + -2.06230817e-02*(u3*delta2*delta) + 7.46948777e-02*(u3*delta2*a1) + 2.18398885e-02*(u4*delta2*delta)

	# mu2
	mu2 = 7.52617641e-01 + -9.85670103e+00*(eta) + -1.01480977e+00*(u) + 4.27950370e-01*(a1) + 2.03450769e+01*(eta*eta) + -1.97176699e+00*(u*a1) + 3.83164478e+00*(u*u) + 4.51759956e+00*(u*eta) + -5.64495227e+01*(u2*eta) + 3.33064276e+00*(u2*u) + 2.98935216e+01*(u*eta*a1) + -1.01812102e+01*(u2*a1) + 1.21010226e+02*(u2*eta*a1) + -3.18312703e+01*(u3*eta) + 1.63186301e+02*(u2*eta*eta) + -8.99151077e+01*(u*eta2*a1) + 7.73227602e+01*(u3*eta*eta) + -3.24658654e+02*(u2*eta2*a1)

	# mu3
	mu3 = -8.60959397e-03*(a1) + -5.53668497e-02*(eta) + -7.68585745e-03*(u) + 8.00589682e-03 + -1.20336889e-02*(u*a1) + 1.96401577e-02*(u*u) + 1.60627973e-01*(u*eta) + -1.03165122e+00*(u*eta*eta) + 2.16651674e-01*(eta2*eta) + 4.28078169e-01*(eta2*a1) + -6.61922983e-02*(u2*a1) + -4.57233546e-01*(u2*eta) + 4.00579723e-02*(u2*u) + 8.48938586e-01*(u*eta2*a1) + 2.07384788e+00*(u*eta2*eta) + 1.39506623e+00*(u2*eta*a1) + 2.31598183e+00*(u2*eta*eta) + -1.02999288e+00*(eta3*a1) + -7.59609171e-01*(u3*eta) + -3.22124184e-02*(u3*a1) + 9.78237920e-01*(u3*eta*a1) + 4.70814374e+00*(u3*eta*eta) + -3.03998184e+00*(u2*eta2*eta) + -2.66197460e+00*(u*eta3*a1) + 1.55449768e-01*(u4*eta) + -7.91889503e+00*(u2*eta2*a1) + -9.23909713e+00*(u3*eta2*eta) + -7.24576506e-01*(u4*eta*eta) + -2.20932434e-01*(u4*eta*a1) + -7.44672439e+00*(u3*eta2*a1) + 1.32586169e+01*(u2*eta3*a1) + 1.62561683e+01*(u3*eta3*a1) + 1.02323968e+00*(u4*eta2*a1)

	# mu4
	mu4 = 3.55546031e-02*(a1) + 5.16550887e-03*(delta) + -2.98142967e-03 + -4.05905735e-02*(u*delta) + -1.53144941e-01*(u*u) + -1.23190631e-01*(eta*a1) + -4.32606864e-02*(eta*delta) + 1.67173841e-01*(u*eta*delta) + 2.47962912e-03*(u*eta*a1) + 2.98251687e-01*(u2*a1) + 1.74096994e-01*(u2*delta) + 6.04765304e-01*(u2*eta) + 6.08695422e-02*(eta*delta*a1) + 6.75123634e-02*(u*delta*a1) + -2.50605559e-01*(u2*u) + -1.19594397e+00*(u2*eta*a1) + -6.16274455e-01*(u2*eta*delta) + -3.93831468e-01*(u2*delta*a1) + -2.14722411e-01*(u*eta*delta*a1) + 1.01348020e+00*(u3*eta) + 3.29266858e-01*(u3*delta) + 1.80913935e-01*(u3*a1) + -1.10436405e+00*(u3*eta*delta) + -7.36341736e-01*(u3*eta*a1) + -2.84160523e-01*(u3*delta*a1) + 1.39448403e+00*(u2*eta*delta*a1) + 1.06496025e+00*(u3*eta*delta*a1)

	# nu4
	nu4 = -1.64934289e+00*(u) + 9.47790734e-02 + 7.55829884e-01*(a1) + 2.72021835e-01*(delta) + 1.50235188e+00*(u*delta) + 3.89482860e-01*(u*a1) + -3.22961842e+00*(eta*a1) + -2.88159240e+00*(eta*delta) + 6.61006390e+00*(u*eta) + -1.16387218e+00*(delta*a1) + 2.76872675e-02*(u*u) + -3.82904186e+00*(u*eta*delta) + -3.32192252e-01*(u2*delta) + 7.79470670e+00*(eta*delta*a1) + -1.61914959e+00*(u*eta*a1) + 1.09847919e+00*(u*eta*delta*a1) + 5.34915446e-01*(u2*delta*a1) + 2.02660452e+00*(u2*eta*delta) + -5.10026036e+00*(u2*eta*delta*a1)

	# nu5
	nu5 = 4.40492301e-02 + -1.77614449e-01*(eta) + -7.13763526e-02*(a1) + -4.22210032e-02*(delta) + -1.54015140e-02*(u*delta) + -1.41964999e-01*(u*a1) + 2.96035592e-01*(eta*a1) + 1.15318857e-01*(eta*delta) + -2.84119976e-03*(u*eta) + 1.02609585e-01*(delta*a1) + 8.37520491e-02*(u*eta*delta) + -1.01362609e-02*(u2*delta) + -2.44147059e-01*(eta*delta*a1) + -6.34768177e-02*(u2*a1) + 1.99600495e-01*(u*delta*a1) + 5.71798724e-01*(u*eta*a1) + -6.16465348e-01*(u*eta*delta*a1) + 6.73418910e-02*(u2*delta*a1) + 4.61836749e-02*(u2*eta*delta) + 2.46237820e-01*(u2*eta*a1) + -2.28548092e-01*(u2*eta*delta*a1)

	# nu6
	nu6 = -2.58952013e-02 + 1.08461738e-01*(eta) + 6.39265152e-02*(a1) + 4.52848813e-02*(delta) + -1.07984875e-02*(u*delta) + -7.27725570e-02*(u*a1) + -2.56188395e-01*(eta*a1) + -1.92275364e-01*(eta*delta) + -7.24147966e-02*(delta*a1) + 6.64856815e-02*(u*eta*delta) + 2.92513863e-01*(eta*delta*a1) + -1.22009461e-01*(u2*a1) + 9.91815770e-02*(u*delta*a1) + 2.90020271e-01*(u*eta*a1) + -3.07699566e-01*(u*eta*delta*a1) + 1.30323999e-01*(u2*delta*a1) + 4.88980538e-01*(u2*eta*a1) + -4.28828295e-01*(u2*eta*delta*a1)

	#
	return mu0,mu1,mu2,mu3,mu4,nu4,nu5,nu6
