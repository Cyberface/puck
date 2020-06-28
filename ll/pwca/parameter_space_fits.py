

#
def generate_model_params(theta,eta,chi_eff,chi_p):

	'''
	Hola, soy un codigo escribido por "4b_document_fits.py". 
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
	mu1 = -3.01407019e-01*(eta) + -1.37264203e-01*(u) + 8.21783650e-01*(chi_p) + 5.81746290e-01*(chi_eff) + -2.11470765e-01 + -2.77178535e+00*(u*chi_eff) + 2.58972934e-01*(u*chi_p) + 9.17067335e-01*(u2) + -3.26671477e-01*(u*eta) + 2.24885245e+00*(chi_eff2) + -1.38319820e+00*(chi_eff*chi_p) + 5.09687484e-01*(eta2) + 8.06726804e-01*(eta*chi_p) + 8.51549729e-01*(eta*chi_eff) + -7.08101629e-01*(chi_p2) + -4.22240808e-01*(eta*chi_p2) + 2.52361035e+00*(u*eta2) + 6.77922225e+00*(u*eta*chi_eff) + -7.22958673e-01*(u*eta*chi_p) + -8.26807061e-01*(eta2*chi_p) + -4.48410141e+00*(eta2*chi_eff) + 6.93749654e-01*(chi_eff*chi_p2) + -3.39765534e-01*(u2*chi_p) + -7.71881168e-02*(u2*chi_eff) + -4.07972854e+00*(eta*chi_eff2) + 9.23192182e-01*(eta*chi_eff*chi_p) + -1.03972061e-01*(u*chi_eff2) + -2.19814456e+00*(u2*eta)

	# mu2
	mu2 = -3.17047091e+00*(eta) + -3.54703616e-01*(u) + 3.54416040e-01 + -6.30211042e-01*(chi_eff) + 1.02428722e+01*(u*chi_eff) + -2.24810010e+00*(u2) + -9.73560681e+00*(chi_eff2) + 1.29151513e+00*(u*chi_p) + 6.62433906e-01*(chi_p2) + -2.89917748e+01*(u*eta*chi_eff) + -4.66536160e+00*(eta*chi_p2) + 1.55836612e+00*(u*eta2) + 2.52562973e+01*(eta*chi_eff2) + 1.31810241e+01*(eta2*chi_p) + -9.88565681e-01*(u*chi_p2) + -9.15771665e-01*(chi_eff2*chi_p) + 6.40210250e+00*(u2*eta)

	# mu3
	mu3 = 5.04031872e-03 + -4.35378868e-02*(eta) + -1.82668033e-02*(chi_eff) + 9.69908669e-02*(eta2) + 2.81733808e-02*(u*chi_p) + -6.63045912e-02*(eta*chi_eff) + -8.81683968e-03*(chi_eff2) + 2.52502946e-02*(u*chi_eff) + -5.51634884e-03*(u2) + -2.47515355e-02*(u*chi_eff*chi_p) + -4.27364782e-02*(u*eta2) + 6.67260731e-02*(u*eta*chi_eff) + -6.26950371e-02*(u2*chi_eff*chi_p) + -1.47774183e-01*(u*eta*chi_p2) + 9.29792805e-01*(eta2*chi_eff*chi_p) + -1.02458159e-01*(chi_eff2*chi_p2) + 7.49284664e-01*(u*eta*chi_eff2) + 8.97254239e-01*(u2*eta2*chi_eff) + 9.83173344e-01*(eta*chi_eff2*chi_p2) + -3.63160890e+00*(u*eta2*chi_eff2) + -1.92526063e+00*(eta2*chi_eff2*chi_p)

	# mu4
	mu4 = 8.21950311e-03 + 5.95876100e-01*(eta*chi_eff) + -1.84251160e-01*(u*eta) + 1.90871640e-01*(eta*chi_p2) + 1.92457722e-01*(u*eta2) + -1.49869850e-01*(u2*eta) + -8.50697425e-01*(eta*chi_eff*chi_p) + -7.81395498e-01*(eta2*chi_p) + 6.32866877e-02*(u*chi_p2) + -4.28462320e-01*(eta*chi_eff2*chi_p) + 2.25958528e-01*(u2*eta2) + 4.65212761e-02*(u2*chi_p2) + 1.06750735e-01*(u2*eta*chi_p) + -3.83731940e-01*(u2*eta*chi_eff)

	# nu4
	nu4 = -1.27683166e-01*(u) + 9.79369227e-02 + 4.28446033e-01*(u2*chi_eff) + -9.31671992e+00*(eta*chi_eff2) + 1.82501010e+01*(eta*chi_eff*chi_p) + 2.81828679e+00*(u*chi_eff*chi_p) + 2.07620911e+00*(u*chi_eff2) + -3.11640604e+00*(chi_eff2*chi_p) + 2.02345179e+00*(chi_eff*chi_p2) + -3.83354427e+00*(u*eta2*chi_p) + 4.11472807e+01*(eta*chi_eff2*chi_p) + 2.52459470e+00*(u2*eta2) + -2.41539084e+01*(eta*chi_eff*chi_p2) + -6.08238521e+00*(u2*chi_eff*chi_p) + -2.35709378e+00*(u2*chi_p2) + -3.64493423e+01*(eta2*chi_eff2)

	# nu5
	nu5 = 9.14138071e-03 + -2.94567372e-03*(u) + -1.18074112e-02*(eta) + 7.67635205e-02*(chi_eff) + -9.11260225e-02*(u*chi_p) + 1.75457734e-02*(chi_p2) + 1.66158807e-01*(u2*eta) + -2.77168742e-01*(u*eta2) + -7.52532213e-01*(u*eta*chi_eff) + -2.89958876e-01*(eta2*chi_p) + 1.63804825e+00*(eta2*chi_eff) + 1.05742809e-01*(u*chi_p2) + 1.02744882e-01*(chi_eff2*chi_p) + 2.28269379e+00*(u*eta2*chi_eff) + -6.05187034e-01*(u2*eta*chi_eff) + -6.86935446e-01*(u2*eta2) + -7.79522093e-01*(eta*chi_eff*chi_p2)

	# nu6
	nu6 = -3.20851415e-01*(eta) + -5.44757611e-02*(chi_p) + 2.53922300e-02*(chi_eff) + 3.41882256e-02 + 2.45280954e-02*(u*chi_eff) + 7.54272075e-01*(eta2) + 6.10029612e-01*(eta*chi_p) + -5.12743822e-01*(eta*chi_eff) + -2.03119024e-01*(u*eta2) + -1.57204892e+00*(eta2*chi_p) + 2.54196220e+00*(eta2*chi_eff) + -1.08357681e-02*(u2*chi_eff) + 4.07229888e-02*(u2*eta) + 7.94970508e-01*(eta*chi_eff*chi_p) + -4.13931883e-01*(u*chi_eff*chi_p) + 9.75857559e-02*(u*chi_eff2) + 5.60522393e-01*(chi_eff2*chi_p) + -1.51867857e+00*(u*eta2*chi_eff) + 3.29995972e-01*(u*eta2*chi_p) + -2.86917291e+00*(eta*chi_eff2*chi_p) + 1.01821816e-01*(u2*eta*chi_eff) + 2.92162296e+00*(eta2*chi_eff2) + -4.13168313e+00*(eta2*chi_eff*chi_p) + -4.35665355e-02*(u2*chi_eff2) + -2.02383388e-01*(u2*chi_eff*chi_p) + 1.80421855e+00*(u*eta*chi_eff*chi_p)

	#
	return mu1,mu2,mu3,mu4,nu4,nu5,nu6
