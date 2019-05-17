#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <vector>

//g++ -std=c++11 -o main main.cpp

const double TWOPI = 2.0*acos(-1.0); //2*PI
const double PI=TWOPI/2;
const double G=6.673*pow(10,-11); //Gravitationskonstante
const double M=1.99*pow(10,30); //Zentralmasse
const double AE=1.49*pow(10,11); //eine Astronomische Einheit
const double a1=90*AE; //Große Bahnhalbachse
const double a2=100*AE;

using namespace std;

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

double g(const double &theta1, const double &theta2) {
  double e1=1.1;
  double e2=1.2;
  double omega1=0;
  double omega2=0;
   double integrand=-(G*M*(1 + e1*cos(theta1))*(1 + e2*cos(theta2))*sin(theta1 - theta2 + omega1 - omega2))/(4.*pow(a1,2)*pow(a2,2)*sqrt(1 - pow(e1,2))*sqrt(1 - pow(e2,2))*pow(PI,2)*abs(1 + e1*cos(theta1))*abs(1 + e2*cos(theta2))*pow(pow(1 + e1*cos(theta1),2)/(pow(a1,2)*pow(-1 + pow(e1,2),2)) + pow(1 + e2*cos(theta2),2)/(pow(a2,2)*pow(-1 + pow(e2,2),2)) - (2*(1 + e1*cos(theta1))*(1 + e2*cos(theta2))*cos(theta1 - theta2 + omega1 - omega2))/(a1*a2*(-1 + pow(e1,2))*(-1 + pow(e2,2))),1.5));
   
   return integrand;
}

double f(const double &theta1) {
  return Integrate(0.0, TWOPI, std::bind(g, theta1, std::placeholders::_1));
}



int main() {
  auto Ergebnis = Integrate(0.0, TWOPI, f);
  size_t Anzahl = 2;
  std::vector<double> omega {0.0, 1.0, 2.0}, ecc {0.0, 0.1, 0.2};
  std::vector<double> omega_dot = omega, ecc_dot = ecc;
  cout << Ergebnis << endl;
}