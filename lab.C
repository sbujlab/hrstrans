#include <TMath.h>

double HRS = 0.0872665;//5 degree central angle


double a = sin(HRS);
double b = cos(HRS);
double c = tan(HRS);

//sign = -1 for LHRS and +1 for RHRS

double thlab(double th, double ph, int sign){

double cos_lab = ( b + sign*ph*a)/sqrt( 1 + th*th + ph*ph );

return acos(cos_lab);

}


double philab(double th, double ph, int sign){ 

//Using cosine definition

double num = -sign*th;

double den = sqrt( th*th + a*a - 2*sign*ph*c + ph*ph/b*b );

double phi_lab = acos(num/den);

return phi_lab;

}


double jacobian(double th, double ph, int sign){

double theta = thlab(th,ph,sign); double phi = philab(th,ph,sign);

double d = 1.0/tan(phi);
double e = -1.0/sin(phi);

double D1 = sqrt( 1 + th*th + ph*ph );
double D2 = sqrt( th*th + a*a - 2*sign*ph*c + ph*ph/b*b );

double A = -th*cos(theta)/D1*D1;
double B = e*( -sign*(1.0 - pow(cos(phi),2.0)) )/D2;
double C = ( sign*D1*a - ph*cos(theta) )/D1*D1;
double D = d*( ph/b*b - sign*c )/D2*D2;

return fabs(A*D - B*C);

}

double Q2(double E, double th, double ph, int A, int sign){

double scat_theta = thlab(th,ph,sign);

double amu = 931.494*1e-3; //amu to GeV
double M = A*amu;

double ratio = E/M;

double recoil = 1.0 + ratio*(1.0-cos(scat_theta));
double Eprime = Energy/recoil;

double q2 = 2*E*Eprime*(1.0-cos(scat_theta));

return q2;//This Q^2 in energy


}


double crsMott( double Z, double E, double th, double ph, int A, int sign){ 

double halfth = thlab(th,ph,sign)/2.0;
double qsq = Q2(E,th,ph,A,sign);
double amu = 931.494*1e-3; //amu to GeV
double M = A*amu;

double cth2 = cos(halfth)*cos(halfth);
double Ef = E*M/(M+E*(1.0-cos(thlab(th,ph,sign))));

double num = 4.0*Z*Z*(0.197*0.197)*Ef*Ef*Ef*cth2;
double den = qsq*qsq*137.0*137.0*E;


return num*1e-2/den;// from fm2 to barns



}
