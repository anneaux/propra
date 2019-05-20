#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <vector>

//g++ -std=c++11 -o main main.cpp

using namespace std;

const double TWOPI = 2.0 * acos(-1.0); //2*Pi
const double PI = TWOPI / 2.0; //Pi
const double G = 1.0; //Gravitationskonstante
const double M = 1.0; //Zentralmasse
const double a1 = 90.0; //Große Bahnhalbachse
const double a2 = 100.0;

// Integration functions adapted from Press et al. (Numerical recipes in C++)

double LogIntegrate(double A, double B,
                    const std::function<double(double)> &Integrand,
                    double prec = 1e-4) {
  double x,tnm,del, exp_x;
  double psum, ps,pst,post=0.0,pos=0.0;
  double f1, f2, a = log(A), b = log(B);
  int it,j,i, min = 3, max = 20;
  for (i=1;i<=max; ++i) {
    if (i == 1) {
      f1 = Integrand(A)*A; f2 = Integrand(B)*B;
      if (not std::isfinite(f1)) f1 = 0.0;
      if (not std::isfinite(f2)) f2 = 0.0;
      pst = 0.5*(b - a)*(f1 + f2);
    }
    else {
      for (it = 1, j = 1; j < i - 1; ++j) it <<= 1;
      tnm=it;
      del=(b-a)/tnm; // This is the spacing of the points to be added.
      x=a+0.5*del;
      for (psum=0.0,j=1;j<=it;++j,x+=del) {
        exp_x = exp(x);
        f1 = Integrand(exp_x)*exp_x;
        if (not std::isfinite(f1)) f1 = 0.0;
        psum += f1;
      }
      pst = 0.5*(post + (b-a)*psum/tnm);
    }
    ps=(4.0*pst-post)/3.0; // Compare equation (4.2.4), above.
    if (i > min) { // Avoid spurious early convergence.
      if ((fabs(ps-pos) < prec*fabs(pos)) || (ps == 0.0 && pos == 0.0))
        return ps;
    }
    pos=ps;
    post=pst;
  }
  return ps;
}

double Integrate(double a, double b,
                 const std::function<double(double)> &Integrand,
                 double prec = 1e-4) {
  double x,tnm,del;
  double psum, ps,pst,post=0.0,pos=0.0;
  double f1, f2;
  int it,j,i, min = 3, max = 20;
  for (i = 1; i <= max; ++i) {
    if (i == 1) {
      f1 = Integrand(a); f2 = Integrand(b);
      pst = 0.5*(b - a)*(f1 + f2);
    }
    else {
      for (it = 1, j = 1; j < i - 1; ++j) it <<= 1;
      tnm=it;
      del=(b - a)/tnm; // This is the spacing of the points to be added.
      x = a+0.5*del;
      for (psum = 0.0, j = 1; j <= it; ++j, x += del) {
        f1 = Integrand(x);
        psum += f1;
      }
      pst = 0.5*(post + (b - a)*psum/tnm);
    }
    ps=(4.0*pst - post)/3.0; // Compare equation (4.2.4), above.
    if (i > min) { // Avoid spurious early convergence.
      if ((fabs(ps - pos) < prec*fabs(pos)) || (ps == 0.0 && pos == 0.0))
        return ps;
    }
    pos=ps;
    post=pst;
  }
  return ps;
}

double g_omega(const double &theta1, const double &theta2, const double &ecc1, const double &ecc2, const double &omega1, const double &omega2) {
  
	double var1=1.0+ecc1*cos(theta1);
	double var2=1.0+ecc2*cos(theta2);
	double var12=var1*var1;
	double var22=var2*var2;
	double ecc12=ecc1*ecc1;
	double ecc22=ecc2*ecc2;
	double a12=a1*a1;
	double a22=a2*a2;
	double deltal=theta1-theta2+omega1-omega2;
	
	double integrand = -(G * M * var1 * var2 * sin(deltal)) / (4.0 * a12* a22 * sqrt(1.0 - ecc12) * sqrt(1.0 - ecc22) * PI*PI * abs(var1) * abs(var2) * pow(var12 / (a12 * (-1.0 + ecc12)*(-1.0 + ecc12)) + var22 / (a22 * (-1.0 + ecc22)*(-1.0 + ecc22)) - (2.0 * var1 * var2 * cos(deltal)) / (a1 * a2 * (-1.0 + ecc12) * (-1.0 + ecc22)), 1.5));
	
	return integrand;
}

double f_omega(const double &theta_1, const double &ecc_1, const double &ecc_2, const double &omega_1, const double &omega_2) {
	double out= Integrate(0.0, TWOPI, std::bind(g_omega, theta_1, placeholders::_1, ecc_1, ecc_2, omega_1, omega_2)); //hier passiert die theta2-Integration!
	return out;  
}

double R_omega(const double &ECCENTRICITY1, const double &ECCENTRICITY2, const double &OMEGA1, const double &OMEGA2) {
	double result = Integrate(0.0, TWOPI, std::bind(f_omega, placeholders::_1, ECCENTRICITY1, ECCENTRICITY2, OMEGA1, OMEGA2));  //hier passiert die theta1-Integration!
	return result;
}


double g_e(const double &theta1, const double &theta2, const double &ecc1, const double &ecc2, const double &omega1, const double &omega2) {
  
	double var1=1.0+ecc1*cos(theta1);
	double var2=1.0+ecc2*cos(theta2);
	double var12=var1*var1;
	double var22=var2*var2;
	double ecc12=ecc1*ecc1;
	double ecc22=ecc2*ecc2;
	double a12=a1*a1;
	double a22=a2*a2;
	double deltal=theta1-theta2+omega1-omega2;
	
	double integrand = -(G*M*var1*
      (a2*(1.0 - ecc22)*var1*
         (2.0*ecc1 + (1.0 + ecc12)*cos(theta1))*
         (a2 - a2*ecc22 - 
           a2*ecc1*(-1.0 + ecc22)*cos(theta1) + 
           a1*(-1.0 + ecc12)*var2*
            cos(deltal)) + 
        (1.0 - ecc12)*cos(theta1)*
         (a22*pow(-1.0 + ecc22,2)*
            var12 + 
           a12*pow(-1.0 + ecc12,2)*
            var22 - 
           2*a1*a2*(-1.0 + ecc12)*
            (-1.0 + ecc22)*var1*
            var2*cos(deltal)) + 
        ecc1*var1*
         (a22*pow(-1.0 + ecc22,2)*
            var12 + 
           a12*pow(-1.0 + ecc12,2)*
            var22 - 
           2.0*a1*a2*(-1.0 + ecc12)*
            (-1.0 + ecc22)*var1*
            var2*cos(deltal))))/
   (4.0*pow(a1,3)*pow(a2,3)*
     pow(1.0 - ecc12,2.5)*
     pow(1.0 - ecc22,1.5)*pow(PI,2)*
     pow(var12,1.5)*
     sqrt(var22)*
     pow(var12/
        (a12*pow(-1.0 + ecc12,2)) + 
       var22/
        (a22*pow(-1.0 + ecc22,2)) - 
       (2.0*var1*var2*
          cos(deltal))/
        (a1*a2*(-1.0 + ecc12)*(-1.0 + ecc22)),
      1.5));
	
	return integrand;
}

double f_e(const double &theta_1, const double &ecc_1, const double &ecc_2, const double &omega_1, const double &omega_2) {
	double out= Integrate(0.0, TWOPI, std::bind(g_e, theta_1, placeholders::_1, ecc_1, ecc_2, omega_1, omega_2)); //hier passiert die theta2-Integration!
	return out;  
}

double R_e(const double &ECCENTRICITY1, const double &ECCENTRICITY2, const double &OMEGA1, const double &OMEGA2) {
	double result = Integrate(0.0, TWOPI, std::bind(f_e, placeholders::_1, ECCENTRICITY1, ECCENTRICITY2, OMEGA1, OMEGA2));  //hier passiert die theta1-Integration!
	return result;
}


//+++ LOG +++
//[NOTIZ von Pierre am 04.05.2019, 18:20] Code läuft, ist aber extrem langsam, selbst mit sehr viel höherer Fehlertoleranz!
//[NOTIZ von Pierre am 06.05.2019, 01:24] Habe die Eingangsparameter mal auf andere Größenordnungen gesetzt (Sonnenmasse=1, 1 AE=100...), scheint den Code allerding auch nicht zu beschleunigen..

int main() {
  	cout << g_omega(0.1, 0.2, 0.2, 0.21, 0.34, 0.36) << endl;
	cout << g_e(0.1, 0.2, 0.2, 0.21, 0.34, 0.36) << endl;
	cout << R_omega(0.2, 0.21, 0.34, 0.36) << endl;
	cout << R_e(0.2, 0.21, 0.34, 0.36) << endl;
}

//  size_t Anzahl = 2;
//	std::vector<double> omega{ 0.0, 1.0, 2.0 }, ecc{ 0.0, 0.1, 0.2 };
//	std::vector<double> omega_dot = omega, ecc_dot = ecc;