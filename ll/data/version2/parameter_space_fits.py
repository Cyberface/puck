

#
def generate_model_params(theta,eta,chi_eff,chi_p):

	'''
	Hola, soy un codigo escribido por codigo. ~londonl@mit.edu/pilondon2@gmail.com 2020
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
	mu1 = lambda u,eta,chi_eff,chi_p: 1.10537213e-01*(u) + -4.38019446e-01*(u*chi_p) + -9.44045919e-01*(u2) + 3.28897407e-02*(u*eta) + 5.09971956e+00*(eta2) + -1.24850868e+00*(eta*chi_eff) + -2.91941031e+01*(eta2*chi_p) + 3.89207583e-01*(u*chi_p2) + 2.15169710e-01*(u2*chi_eff) + 1.89358546e+01*(u2*eta) + 3.62194064e+00*(eta*chi_eff*chi_p) + -3.81181189e-02*(u*chi_eff*chi_p) + 3.56406791e+00*(u2*chi_p) + -1.31032507e+01*(eta3) + -7.30541129e-01*(u*eta3) + -7.21836073e+01*(u2*eta*chi_p) + -1.37413033e+02*(u2*eta2) + 1.53943581e+00*(eta2*chi_eff*chi_p) + -3.01346957e+00*(u2*chi_p2) + 8.60554145e+01*(eta3*chi_p) + 3.65954641e+01*(eta2*chi_p2) + -1.16858023e+02*(eta3*chi_p2) + 2.80070974e+02*(u2*eta3) + 5.42336786e+02*(u2*eta2*chi_p) + 6.20099998e+01*(u2*eta*chi_p2) + -1.36386424e+00*(u2*chi_eff*chi_p2) + -9.37938538e+00*(eta3*chi_eff*chi_p) + -1.12518583e+03*(u2*eta3*chi_p) + -4.84338920e+02*(u2*eta2*chi_p2) + 1.02643193e+03*(u2*eta3*chi_p2)

	# mu2
	mu2 = lambda u,eta,chi_eff,chi_p: 6.55030950e-01 + -4.50593505e+00*(eta) + -4.49780862e-01*(chi_p) + -2.20057275e+00*(chi_eff) + -1.15747920e+01*(eta2) + 5.49253634e+00*(chi_eff*chi_p) + -9.33964152e+00*(eta*chi_eff*chi_p) + 1.53613292e+01*(eta2*chi_eff) + 5.69698747e+01*(eta2*chi_p) + -7.07139293e-01*(u*chi_p2) + -5.32246860e+00*(chi_eff*chi_p2) + -3.55868537e+01*(eta2*chi_p2) + 3.14907103e+00*(u*eta*chi_p2) + 8.41697675e+00*(u*eta*chi_eff*chi_p2) + 4.55690308e+00*(u2*chi_eff*chi_p2) + 1.93030524e+01*(eta2*chi_eff2*chi_p)

	# mu3
	mu3 = lambda u,eta,chi_eff,chi_p: -2.39102315e-02 + -2.78245101e-01*(u*chi_eff) + -1.34573701e-01*(u2) + 4.38103300e-01*(u*eta) + -2.13414327e-01*(chi_eff*chi_p) + -1.13809660e+00*(eta*chi_eff) + -3.82634368e-01*(chi_p2) + 2.21684438e+00*(eta*chi_p2) + 4.24699810e+00*(u*eta*chi_eff) + -1.86856495e+00*(u*eta*chi_p) + 1.63565416e+00*(eta2*chi_p) + 3.26576390e+00*(eta2*chi_eff) + 1.05368697e-01*(u*chi_p2) + 8.17016421e-01*(u2*chi_p) + 7.56352218e-01*(u2*chi_eff) + 6.05445054e-01*(u2*eta) + 2.26292344e+00*(eta*chi_eff*chi_p) + 2.28458362e+00*(u*chi_eff2) + -1.11315744e+01*(u*eta2*chi_eff) + 1.75833187e+00*(u*eta2*chi_p) + -3.54959851e+00*(eta*chi_eff2*chi_p) + -1.66777203e+00*(u2*eta2) + 4.94232549e-01*(u*eta*chi_p2) + 4.14419198e-01*(u*chi_eff*chi_p2) + -1.87296126e+00*(eta2*chi_eff2) + 3.47453905e+00*(eta2*chi_eff*chi_p) + -3.76844529e-01*(u2*chi_p2) + -1.66681368e+00*(u2*eta*chi_p) + -6.35060190e-01*(u2*chi_eff2) + -1.40470591e+00*(u2*chi_eff*chi_p) + -1.99792758e+00*(u*eta*chi_eff*chi_p) + 3.62970135e+00*(u*eta*chi_eff2) + -3.29916488e+00*(eta*chi_eff*chi_p2) + -4.51234630e+00*(eta2*chi_p2) + -5.71633049e+00*(u*chi_eff2*chi_p) + 2.04697921e+00*(chi_eff2*chi_p2)

	# mu4
	mu4 = lambda u,eta,chi_eff,chi_p: 1.50751469e-02*(chi_p) + -2.66809634e-01*(eta2) + 7.46132339e-03*(u*chi_p) + 6.64365335e-03*(u2) + -3.31775516e-02*(u*eta) + -3.06991258e+00*(eta2*chi_eff) + 1.18469761e-01*(chi_eff*chi_p2) + 5.14121283e-01*(eta*chi_eff2*chi_p) + -8.47238517e-01*(u*eta*chi_eff2) + -7.96010583e-01*(eta*chi_eff*chi_p2) + -1.06248147e-01*(u2*chi_eff*chi_p) + 4.75324225e-01*(u*eta2*chi_eff) + 4.56125389e+00*(eta2*chi_eff*chi_p)

	# nu4
	nu4 = lambda u,eta,chi_eff,chi_p: -7.03539714e-01 + 2.72898145e+01*(eta) + 1.69264112e+01*(u*chi_eff) + -1.22548479e+01*(u*chi_p) + 1.88395158e+01*(chi_eff*chi_p) + -4.94180151e+00*(u*eta) + -9.46993023e+01*(eta2) + -4.31543853e+01*(eta*chi_eff) + -1.98072984e+02*(u*eta*chi_eff) + 2.13105749e+01*(u*eta*chi_p) + 3.57208039e+01*(eta2*chi_p) + 1.59389568e+02*(eta2*chi_eff) + 1.97347608e+01*(u2*chi_p) + 5.87685546e+01*(u2*chi_eff) + -1.26149126e+02*(u2*eta) + -5.36744347e+01*(u*chi_eff*chi_p) + 2.50994986e+02*(u*eta2*chi_eff) + 1.94054765e+02*(u*eta2*chi_p) + -3.69386692e+02*(u2*eta*chi_eff) + 4.84921041e+02*(u2*eta2) + -4.88054798e+02*(eta2*chi_eff*chi_p) + -2.39670483e+02*(u2*chi_eff*chi_p) + 1.73061279e+02*(u*eta*chi_eff*chi_p) + -7.36385334e+01*(u3*chi_eff) + 2.25964930e+01*(u3*chi_p) + -3.72284762e+02*(u2*eta2*chi_p) + 9.65926420e+02*(u2*eta2*chi_eff) + 1.16351741e+02*(u3*eta*chi_eff) + 5.64290993e+01*(u3*eta*chi_p) + 6.98705007e+01*(u3*eta2) + 1.18334008e+03*(u*eta2*chi_eff*chi_p) + 3.14959937e+02*(u3*chi_eff*chi_p) + 2.27003500e+03*(u2*eta*chi_eff*chi_p) + 7.08540335e+02*(u3*eta2*chi_eff) + -7.64157523e+02*(u3*eta2*chi_p) + -5.68261294e+03*(u2*eta2*chi_eff*chi_p) + -1.48495260e+03*(u3*eta*chi_eff*chi_p)

	# nu5
	nu5 = lambda u,eta,chi_eff,chi_p: -4.20609605e-03 + 1.31538883e-01*(chi_p) + -2.11277253e+00*(eta*chi_p) + 3.93432001e-01*(eta2) + -1.52212181e+00*(eta3) + 1.03455962e+01*(eta2*chi_p) + 3.69023184e+00*(eta2*chi_eff) + 6.26333388e-02*(u2*chi_p) + -1.63082057e-02*(u3) + -3.78630860e+00*(u*eta2*chi_eff) + -1.41697997e+00*(u*eta2*chi_p) + -1.52338017e+00*(u2*eta*chi_eff) + -6.43356213e+00*(eta2*chi_eff*chi_p) + -1.46320753e+01*(eta3*chi_p) + -2.74012620e+01*(eta3*chi_eff) + 3.16208504e-01*(u2*chi_eff*chi_p) + 6.25583097e-02*(u3*chi_eff) + 1.26313207e-01*(u3*chi_p) + 3.01028523e-01*(u2*eta3) + -1.05627524e+00*(u2*eta2*chi_p) + 6.15883523e+00*(u2*eta2*chi_eff) + -5.83939491e-01*(u3*eta*chi_p) + 1.95243686e+01*(u*eta3*chi_eff) + 8.28310672e+00*(u*eta3*chi_p) + 3.77407974e+01*(eta3*chi_eff*chi_p) + -6.97429461e-01*(u2*eta*chi_eff*chi_p)

	# nu6
	nu6 = lambda u,eta,chi_eff,chi_p: 2.31288184e-01*(eta) + 6.35552298e-02*(u) + 3.53270593e-02*(chi_p) + -1.90603964e-02 + 3.42350141e-01*(u*chi_eff) + -7.96275825e-02*(u*chi_p) + 1.94460509e-01*(u2) + -5.24684381e-01*(u*eta) + 1.24378783e-01*(chi_eff*chi_p) + -5.11865431e-01*(eta2) + -2.49325024e+00*(eta*chi_eff) + 2.61494121e+00*(u*eta2) + -1.77582054e+00*(u*eta*chi_eff) + 1.46965490e+01*(eta2*chi_eff) + -3.17867328e-01*(u2*chi_p) + 2.54760192e-01*(u2*chi_eff) + -1.54303978e+00*(u2*eta) + 1.79052457e+00*(eta*chi_eff*chi_p) + -4.04370669e-01*(u*chi_eff*chi_p) + 2.34013334e+00*(u*eta2*chi_eff) + 1.17704241e+00*(u*eta2*chi_p) + -5.81294290e+00*(u*eta3) + 1.06529755e+00*(u2*eta*chi_p) + 4.98165544e-01*(u2*eta*chi_eff) + 3.11319122e+00*(u2*eta2) + -6.68812957e+00*(eta2*chi_eff*chi_p) + -2.26557219e+00*(eta3*chi_p) + -2.38460514e+01*(eta3*chi_eff) + -7.28699181e-01*(u2*chi_eff*chi_p) + 1.10187291e+00*(u*eta*chi_eff*chi_p) + 1.40278885e-01*(u3*eta) + -8.44384626e-02*(u3*chi_eff) + -8.09188149e-02*(u3*chi_p)

	#
	return mu1,mu2,mu3,mu4,nu4,nu5,nu6
