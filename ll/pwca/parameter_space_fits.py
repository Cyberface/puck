

#
def generate_model_params(theta,eta,a1):

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
	u4 = u3*u
	eta2 = eta*eta
	eta3 = eta2*eta

	# mu2
	mu2 = -8.03370686e-02  +  2.14824134e-01 * (  1.48722667e+01*(u) + 5.92151057e+01*(u*u) + -1.87204079e+01*(u3*eta) + -1.45601884e+01*(u3*a1) + -5.59838594e+02*(u2*eta) + 1.26172225e+03*(u2*eta*eta) + -4.12118759e+01*(u2*a1) + -8.97837100e+01*(u*eta) + -8.93423324e+00*(u*a1) + 1.88158204e+02*(eta) + -8.04125088e+02*(eta*eta) + 8.89624720e+01*(eta*a1) + -1.13596927e+01 ) / ( 1.0 +  -6.24506982e+01*(u3*eta) + 1.40029553e+01*(u3*a1) + -2.16161622e+02*(u2*eta*a1) + 2.86082684e+01*(u2*a1) + 3.04350043e+02*(u*eta*eta) + -4.34926236e+02*(u*eta2*eta) + 1.83049993e+02*(u*eta2*a1) + 3.10517740e+04*(eta3*a1) + -1.46114665e+04*(eta2*a1) + 2.38922547e+03*(eta*a1) + -1.18893783e+02*(a1) )

	# mu3
	mu3 = -1.24677786e+00*(eta) + -3.19062642e-02*(u) + -5.65811528e-02*(a1) + 6.16688589e-02 + -1.00390700e-02*(u*u) + 6.37446379e-01*(u*eta) + 7.96505517e+00*(eta*eta) + 1.27759521e+00*(eta*a1) + -3.65344242e+00*(u*eta*eta) + -1.15120106e-01*(u*eta*a1) + -1.60721847e+01*(eta2*eta) + -8.77523728e+00*(eta2*a1) + -1.51262297e-02*(u2*a1) + 1.79565182e-01*(u2*eta) + -3.30908259e-03*(u2*u) + 4.43721343e-01*(u*eta2*a1) + 6.57487130e+00*(u*eta2*eta) + 6.37824211e-02*(u2*eta*a1) + -6.09777504e-01*(u2*eta*eta) + 5.16718936e-03*(u3*u) + 6.81634785e-03*(u3*a1) + 1.85286162e+01*(eta3*a1)

	# mu4
	mu4 = 1.39701243e-04  +  1.27087477e-02 * (  -5.63022060e+01*(u*u) + 1.15835934e+01*(u3*u) + -1.31157596e+02*(u2*u2*eta) + 3.61267488e+00*(u2*u2*a1) + -1.57962622e+02*(u3*eta*eta) + -3.46827078e+01*(u3*eta*a1) + 1.14580701e+03*(u2*eta) + -7.12103057e+03*(u2*eta*eta) + 1.39263020e+04*(u2*eta2*eta) + -2.89610224e+02*(u*eta3*a1) + -5.86355518e+01*(eta) + 4.65699978e+02*(eta2*eta) + -5.86517155e+03*(eta3*a1) + 2.02343398e+03*(eta2*a1) + -1.27944819e+02*(eta*a1) + 5.66631638e+00 ) / ( 1.0 +  4.52737693e+00*(u) + -6.54191065e+00*(u2*u) + 2.03076229e+02*(u3*eta*eta) + -5.56536037e+01*(u2*eta*a1) + 1.44521852e+01*(u2*a1) + -1.03079613e+01*(u*a1) + 9.88041902e+01*(eta2*eta) + -2.16689604e+02*(eta2*a1) + 1.07441301e+02*(eta*a1) + -7.89707147e+00*(a1) )

	# nu4
	nu4 = 1.47581741e+02*(u) + 5.56472619e+02*(eta) + 7.79859437e+01*(a1) + -1.20367338e+01 + -1.48807873e+02*(u*a1) + -3.24304321e+03*(u*eta) + -2.04531479e+02*(u*u) + -4.68432076e+03*(eta*eta) + -2.03906689e+03*(eta*a1) + 2.04942854e+04*(u*eta*eta) + 3.65042450e+03*(u*eta*a1) + 1.09383146e+04*(eta2*eta) + 1.52325585e+04*(eta2*a1) + 1.62252253e+02*(u2*a1) + 3.08014467e+03*(u2*eta) + -5.15707000e+01*(u2*u) + -2.31526846e+04*(u*eta2*a1) + -3.97820981e+04*(u*eta2*eta) + -1.64728417e+03*(u2*eta*a1) + -1.49512417e+04*(u2*eta*eta) + 2.85841012e+02*(u3*u) + 1.38100620e+03*(u3*eta) + -3.35719406e+04*(eta3*a1) + 4.00596847e+03*(u2*eta2*a1) + 2.38415298e+04*(u2*eta2*eta) + -8.52787812e+02*(u3*eta*a1) + -9.49060082e+03*(u3*eta*eta) + 4.38174202e+04*(u*eta3*a1) + -4.96204957e+03*(u2*u2*eta) + -4.45118390e+02*(u2*u2*a1) + 1.92206160e+04*(u3*eta2*eta) + 7.12159801e+03*(u3*eta2*a1) + 2.78820525e+04*(u2*u2*eta*eta) + 7.92009596e+03*(u2*u2*eta*a1) + -1.55020445e+04*(u3*eta3*a1) + -5.00792698e+04*(u2*u2*eta2*eta) + -4.64700910e+04*(u2*u2*eta2*a1) + 8.64814441e+04*(u2*u2*eta3*a1)

	# nu5
	nu5 = 2.88393666e-03  +  1.15379345e-02 * (  -9.10049554e+00*(u2*u2*a1) + 2.37675566e+01*(u3*eta) + -3.34125728e+02*(u3*eta*eta) + 1.90868859e+02*(u2*eta*eta) + -8.52787956e+02*(u2*eta2*eta) + -2.68410171e+01*(u*eta) + 1.32023939e+02*(u*eta*eta) + -1.87413293e+03*(u*eta3*a1) + 1.00962716e+03*(u*eta2*a1) + -1.94212974e+02*(u*eta*a1) + 1.32224096e+01*(u*a1) + -1.19677106e+02*(eta*eta) + 4.64863719e+02*(eta2*eta) + -6.23722649e+02*(eta3*a1) + 2.40899888e+01*(eta*a1) + 3.82948929e+00*(a1) + 3.55090472e-02 ) / ( 1.0 +  5.03629896e+00*(u) + -4.13900660e+00*(u2*u2*a1) + 2.80843911e+01*(u2*eta) + -4.74961130e+01*(u*eta) + 4.31117359e+02*(u*eta2*eta) + 9.93083630e+01*(u*eta2*a1) + -5.54397347e+00*(u*a1) + 1.88956541e+02*(eta*eta) + -8.11966421e+02*(eta2*eta) + 1.42214128e+00*(a1) )

	# nu6
	nu6 = 1.65193130e-01*(u) + 1.96529075e-01*(eta) + 8.89560794e-03 + -2.45531301e-01*(u*a1) + -3.39901212e+00*(u*eta) + -2.52412114e+00*(eta*eta) + -4.99256620e-01*(eta*a1) + 2.09634475e+01*(u*eta*eta) + 5.26404661e+00*(u*eta*a1) + 6.49767139e+00*(eta2*eta) + 5.40347638e+00*(eta2*a1) + 6.64006638e-02*(u2*a1) + -8.89632532e-01*(u2*eta) + -2.23570551e-01*(u2*u) + -3.27643303e+01*(u*eta2*a1) + -4.00640084e+01*(u*eta2*eta) + 3.84020955e-01*(u2*eta*a1) + 8.84012578e+00*(u2*eta*eta) + -1.37213270e+01*(eta3*a1) + -4.46915134e-02*(u3*u) + 4.46789123e+00*(u3*eta) + 4.06043368e-01*(u3*a1) + -9.63694391e+00*(u2*eta2*a1) + -2.09695829e+01*(u2*eta2*eta) + -8.35137583e+00*(u3*eta*a1) + -2.74675690e+01*(u3*eta*eta) + 6.24539415e+01*(u*eta3*a1) + 1.83585311e+00*(u2*u2*eta) + 5.26143874e+01*(u3*eta2*eta) + 5.20883880e+01*(u3*eta2*a1) + -1.46867234e+01*(u2*u2*eta*eta) + -1.67352623e+00*(u2*u2*eta*a1) + 2.82025788e+01*(u2*eta3*a1) + -1.00435406e+02*(u3*eta3*a1) + 3.21820832e+01*(u2*u2*eta2*eta) + 1.69463671e+01*(u2*u2*eta2*a1) + -4.09910758e+01*(u2*u2*eta3*a1)

	#
	return mu2,mu3,mu4,nu4,nu5,nu6
