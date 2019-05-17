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
