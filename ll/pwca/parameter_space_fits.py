

#
def generate_model_params(theta,eta,chi_eff,chi_p):

	'''
	Hola, soy un codigo escribido por codigo. 
	~londonl@mit.edu/pilondon2@gmail.com 2020
	'''  

	# Import usefuls
	from numpy import cos

	# Preliminaries
	u = cos(theta)
	u2 = u*u
	u3 = u2*u
	u3 = u2*u
	eta2 = eta*eta
	eta3 = eta2*eta
	chi_eff2 = chi_eff*chi_eff
	chi_eff3 = chi_eff2*chi_eff
	chi_p2 = chi_p*chi_p
	chi_p3 = chi_p2*chi_p

	# mu1
	mu1 = 9.42694308e+01*(u) + -1.47491362e+00*(eta) + 1.00900250e+03*(chi_eff) + 1.40448142e-01 + 2.20256641e+00*(u*chi_eff) + -1.35975946e+03*(u*chi_p) + 1.09454970e+00*(u2) + -4.60157194e+02*(u*eta) + 2.11830270e+02*(chi_eff*chi_p) + 3.75319931e+00*(eta2) + 1.25027901e+00*(eta*chi_p) + 3.94408067e+02*(u*eta2) + 1.24076320e+01*(u*eta*chi_eff) + 3.80979319e+03*(u*eta*chi_p) + -5.02561800e+00*(eta2*chi_p) + -1.27513881e+04*(eta2*chi_eff) + -1.27560967e+01*(u2*chi_p) + -9.46882672e+02*(u2*chi_eff) + -1.41860747e+01*(u2*eta) + -7.70645596e+02*(eta*chi_eff*chi_p) + 8.84314717e+00*(u*chi_eff*chi_p) + -8.28199021e+01*(u*eta2*chi_eff) + 4.67801874e+03*(u*eta2*chi_p) + 9.98064490e+01*(u2*eta*chi_p) + 8.87080364e+02*(u2*eta*chi_eff) + 3.91819270e+01*(u2*eta2) + -2.21911045e+01*(eta2*chi_eff*chi_p) + 4.91038756e+02*(u2*chi_eff*chi_p) + -1.06031209e+02*(u*eta*chi_eff*chi_p) + 8.90979074e+03*(u2*eta2*chi_eff) + 2.81956674e+02*(u*eta2*chi_eff*chi_p) + -1.80067545e+03*(u2*eta*chi_eff*chi_p) + -1.96558617e+02*(u2*eta2*chi_p)

	# mu2
	mu2 = -3.17047091e+00*(eta) + -3.54703616e-01*(u) + 3.54416040e-01 + -6.30211042e-01*(chi_eff) + 1.02428722e+01*(u*chi_eff) + -2.24810010e+00*(u2) + -9.73560681e+00*(chi_eff2) + 1.29151513e+00*(u*chi_p) + 6.62433906e-01*(chi_p2) + -2.89917748e+01*(u*eta*chi_eff) + -4.66536160e+00*(eta*chi_p2) + 1.55836612e+00*(u*eta2) + 2.52562973e+01*(eta*chi_eff2) + 1.31810241e+01*(eta2*chi_p) + -9.88565681e-01*(u*chi_p2) + -9.15771665e-01*(chi_eff2*chi_p) + 6.40210250e+00*(u2*eta)

	# mu3
	mu3 = 1.40825629e-02*(eta) + 1.75505597e+01*(chi_eff) + -9.71374510e-04 + 2.33785797e-01*(u*chi_eff) + -7.54084250e+01*(u*chi_p) + -1.03143175e-01*(u2) + 4.97046913e+01*(u*eta) + 3.00905350e-01*(chi_eff2) + -8.81549118e+02*(chi_eff*chi_p) + -3.89131095e-02*(eta2) + -2.72014531e-03*(chi_p2) + -1.57617012e+02*(u*eta2) + -8.83661602e-01*(u*eta*chi_eff) + -8.69584261e+01*(u*eta*chi_p) + 1.10131462e+03*(u*chi_p2) + -1.52456521e+02*(chi_eff*chi_p2) + 4.49638253e-01*(u2*eta) + 2.94201584e+03*(eta*chi_eff*chi_p) + 9.00059974e+01*(u*chi_eff2) + 1.30344721e+03*(u*eta2*chi_p) + -4.23470308e+00*(eta*chi_eff2*chi_p) + -1.31662735e+02*(u2*eta*chi_eff) + -2.20093692e-01*(u2*eta2) + -4.13391382e+03*(u*eta*chi_p2) + 4.28237014e-02*(u2*chi_eff2) + 4.45523544e+02*(u2*chi_eff*chi_p) + -3.43375787e+02*(eta*chi_eff*chi_p2) + -9.47103537e-01*(chi_eff2*chi_p2) + -9.83686389e+01*(u*chi_eff2*chi_p2) + 7.86157615e+00*(eta2*chi_eff2*chi_p) + 5.43050530e+00*(eta*chi_eff2*chi_p2) + -4.92964911e+02*(u*eta2*chi_eff2) + 4.12313853e+03*(eta2*chi_eff*chi_p2) + -1.32086509e+03*(u2*eta*chi_eff*chi_p) + -1.18832374e+03*(u*eta2*chi_eff2*chi_p) + -3.47624917e+00*(u2*eta2*chi_eff2) + 1.95052133e+02*(u2*eta*chi_eff*chi_p2)

	# mu4
	mu4 = 8.21950311e-03 + 5.95876100e-01*(eta*chi_eff) + -1.84251160e-01*(u*eta) + 1.90871640e-01*(eta*chi_p2) + 1.92457722e-01*(u*eta2) + -1.49869850e-01*(u2*eta) + -8.50697425e-01*(eta*chi_eff*chi_p) + -7.81395498e-01*(eta2*chi_p) + 6.32866877e-02*(u*chi_p2) + -4.28462320e-01*(eta*chi_eff2*chi_p) + 2.25958528e-01*(u2*eta2) + 4.65212761e-02*(u2*chi_p2) + 1.06750735e-01*(u2*eta*chi_p) + -3.83731940e-01*(u2*eta*chi_eff)

	# nu4
	nu4 = 2.49862074e+00*(chi_p) + 2.78745561e+01*(eta) + 2.28566030e+00*(u) + 1.19839921e+01*(chi_eff) + -9.16166242e-01 + 4.41588514e-01*(u*chi_eff) + -1.86216993e+01*(u*chi_p) + -3.68079931e-01*(u2) + -1.86517101e+02*(eta2) + -6.99809788e+01*(eta*chi_p) + 6.53458900e+00*(eta*chi_p2) + -3.14409041e+01*(u*eta2) + 4.34323851e+00*(u*eta*chi_eff) + 3.98640948e+02*(eta3) + 4.48867607e+02*(eta2*chi_p) + -1.28234825e+02*(eta2*chi_eff) + -1.06462010e+00*(u2*chi_p) + -1.55668979e+01*(u2*chi_eff) + 8.03257439e+01*(eta*chi_eff*chi_p) + 1.82588181e+02*(u*eta2*chi_p) + 1.17126594e+01*(u2*eta*chi_p) + 1.87971327e+01*(u2*eta*chi_eff) + -4.90176450e+00*(u2*eta2) + 2.43137107e+01*(u*eta*chi_p2) + 3.04444128e+00*(u*chi_eff*chi_p2) + -4.99693722e+02*(eta2*chi_eff*chi_p) + -9.66381886e+02*(eta3*chi_p) + 1.34660259e+01*(u2*chi_eff*chi_p) + -2.22091266e+01*(u*eta*chi_eff*chi_p) + 5.32314858e+02*(eta3*chi_eff)

	# nu5
	nu5 = 9.14138071e-03 + -2.94567372e-03*(u) + -1.18074112e-02*(eta) + 7.67635205e-02*(chi_eff) + -9.11260225e-02*(u*chi_p) + 1.75457734e-02*(chi_p2) + 1.66158807e-01*(u2*eta) + -2.77168742e-01*(u*eta2) + -7.52532213e-01*(u*eta*chi_eff) + -2.89958876e-01*(eta2*chi_p) + 1.63804825e+00*(eta2*chi_eff) + 1.05742809e-01*(u*chi_p2) + 1.02744882e-01*(chi_eff2*chi_p) + 2.28269379e+00*(u*eta2*chi_eff) + -6.05187034e-01*(u2*eta*chi_eff) + -6.86935446e-01*(u2*eta2) + -7.79522093e-01*(eta*chi_eff*chi_p2)

	# nu6
	nu6 = -3.20851415e-01*(eta) + -5.44757611e-02*(chi_p) + 2.53922300e-02*(chi_eff) + 3.41882256e-02 + 2.45280954e-02*(u*chi_eff) + 7.54272075e-01*(eta2) + 6.10029612e-01*(eta*chi_p) + -5.12743822e-01*(eta*chi_eff) + -2.03119024e-01*(u*eta2) + -1.57204892e+00*(eta2*chi_p) + 2.54196220e+00*(eta2*chi_eff) + -1.08357681e-02*(u2*chi_eff) + 4.07229888e-02*(u2*eta) + 7.94970508e-01*(eta*chi_eff*chi_p) + -4.13931883e-01*(u*chi_eff*chi_p) + 9.75857559e-02*(u*chi_eff2) + 5.60522393e-01*(chi_eff2*chi_p) + -1.51867857e+00*(u*eta2*chi_eff) + 3.29995972e-01*(u*eta2*chi_p) + -2.86917291e+00*(eta*chi_eff2*chi_p) + 1.01821816e-01*(u2*eta*chi_eff) + 2.92162296e+00*(eta2*chi_eff2) + -4.13168313e+00*(eta2*chi_eff*chi_p) + -4.35665355e-02*(u2*chi_eff2) + -2.02383388e-01*(u2*chi_eff*chi_p) + 1.80421855e+00*(u*eta*chi_eff*chi_p)

	#
	return mu1,mu2,mu3,mu4,nu4,nu5,nu6
