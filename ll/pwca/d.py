'''
Imports and clones from PhenomD as copied from Sebastian's phenom python repo.
'''

#
def gamma1(eta,chi):
    """
    gamma 1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
    """
    xi = -1 + chi
    xi2 = xi*xi
    xi3 = xi2*xi
    eta2 = eta*eta

    return 0.006927402739328343 + 0.03020474290328911*eta \
    + (0.006308024337706171 - 0.12074130661131138*eta + 0.26271598905781324*eta2)*xi \
    + (0.0034151773647198794 - 0.10779338611188374*eta + 0.27098966966891747*eta2)*xi2 \
    + (0.0007374185938559283 - 0.02749621038376281*eta + 0.0733150789135702*eta2)*xi3

#
def gamma2(eta,chi):
    """
    gamma 2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
    """
    xi = -1 + chi
    xi2 = xi*xi
    xi3 = xi2*xi
    eta2 = eta*eta

    return 1.010344404799477 + 0.0008993122007234548*eta \
    + (0.283949116804459 - 4.049752962958005*eta + 13.207828172665366*eta2)*xi \
    + (0.10396278486805426 - 7.025059158961947*eta + 24.784892370130475*eta2)*xi2 \
    + (0.03093202475605892 - 2.6924023896851663*eta + 9.609374464684983*eta2)*xi3

#
def gamma3(eta,chi):
    """
    gamma 3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
    """
    xi = -1 + chi
    xi2 = xi*xi
    xi3 = xi2*xi
    eta2 = eta*eta

    return 1.3081615607036106 - 0.005537729694807678*eta \
    + (-0.06782917938621007 - 0.6689834970767117*eta + 3.403147966134083*eta2)*xi \
    + (-0.05296577374411866 - 0.9923793203111362*eta + 4.820681208409587*eta2)*xi2 \
    + (-0.006134139870393713 - 0.38429253308696365*eta + 1.7561754421985984*eta2)*xi3
    

#
def alpha1(eta,chi):
    """
    alpha 1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
    """
    xi = -1 + chi
    xi2 = xi*xi
    xi3 = xi2*xi
    eta2 = eta*eta

    return 43.31514709695348 + 638.6332679188081*eta \
    + (-32.85768747216059 + 2415.8938269370315*eta - 5766.875169379177*eta2)*xi \
    + (-61.85459307173841 + 2953.967762459948*eta - 8986.29057591497*eta2)*xi2 \
    + (-21.571435779762044 + 981.2158224673428*eta - 3239.5664895930286*eta2)*xi3

#
def alpha2(eta,chi):
    """
    alpha 2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
    """
    xi = -1 + chi
    xi2 = xi*xi
    xi3 = xi2*xi
    eta2 = eta*eta

    return -0.07020209449091723 - 0.16269798450687084*eta \
    + (-0.1872514685185499 + 1.138313650449945*eta - 2.8334196304430046*eta2)*xi \
    + (-0.17137955686840617 + 1.7197549338119527*eta - 4.539717148261272*eta2)*xi2 \
    + (-0.049983437357548705 + 0.6062072055948309*eta - 1.682769616644546*eta2)*xi3

#
def alpha3(eta,chi):
    """
    alpha 3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
    """
    xi = -1 + chi
    xi2 = xi*xi
    xi3 = xi2*xi
    eta2 = eta*eta

    return 9.5988072383479 - 397.05438595557433*eta \
    + (16.202126189517813 - 1574.8286986717037*eta + 3600.3410843831093*eta2)*xi \
    + (27.092429659075467 - 1786.482357315139*eta + 5152.919378666511*eta2)*xi2 \
    + (11.175710130033895 - 577.7999423177481*eta + 1808.730762932043*eta2)*xi3

#
def alpha4(eta,chi):
    """
    alpha 4 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
    """
    xi = -1 + chi
    xi2 = xi*xi
    xi3 = xi2*xi
    eta2 = eta*eta

    return -0.02989487384493607 + 1.4022106448583738*eta \
    + (-0.07356049468633846 + 0.8337006542278661*eta + 0.2240008282397391*eta2)*xi \
    + (-0.055202870001177226 + 0.5667186343606578*eta + 0.7186931973380503*eta2)*xi2 \
    + (-0.015507437354325743 + 0.15750322779277187*eta + 0.21076815715176228*eta2)*xi3

#
def alpha5(eta,chi):
    """
    alpha 5 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
    """
    xi = -1 + chi
    xi2 = xi*xi
    xi3 = xi2*xi
    eta2 = eta*eta

    return 0.9974408278363099 - 0.007884449714907203*eta \
    + (-0.059046901195591035 + 1.3958712396764088*eta - 4.516631601676276*eta2)*xi \
    + (-0.05585343136869692 + 1.7516580039343603*eta - 5.990208965347804*eta2)*xi2 \
    + (-0.017945336522161195 + 0.5965097794825992*eta - 2.0608879367971804*eta2)*xi3
