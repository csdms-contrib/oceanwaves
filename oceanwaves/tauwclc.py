import numpy as np

def tauwclc(ub,Tp):
    rho = 1020
    ks = 0.025
    #D50=0.075e-3; D84=2*D50;
    #ks=3*D84; 
    ab=ub * Tp / (2 * np.pi)
#    fw = np.zeros(ub.shape())
#    i100 = find(ab/ks>=100); fw(i100)=0.04*(ab(i100)/ks).^(-1/4);
#    i10=find(ab/ks<100&ab/ks>10); fw(i10)=0.4*(ab(i10)/ks).^(-3/4);
#    i1=find(ab/ks<10); fw(i1)=(0.4*10.^(-3/4));
    fw = 0.04 * (ab / ks) ** -0.25
    tauw = rho * fw / 2 * ub ** 2
    return tauw