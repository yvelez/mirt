

/////////////////////////////////////

#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

double a1_cpp(
    arma::colvec par,
    arma::mat Theta,
    const int D,
    const int C,
    arma::colvec Z,
    const int dref) {

  const int J = Z.n_rows ;
  const int M = 2*C+1 ;
  double tau = 0 ;
  double tau_prime = 0 ;
  double dentau = 0 ;
  double dentau_prime = 0 ;
  double sumtau = 0 ;
  double dist = 0 ;
  double sumdist = 0 ;
  double num1 = 0 ;
  double num2 = 0 ;
  double num1_prime = 0 ;
  double num2_prime = 0 ;
  double num1_z_prime = 0 ;
  double num2_z_prime = 0 ;
  double a = 0 ;
  double a_prime = 0 ;
  double b = 0 ;
  double b_prime = 0 ;
  double c_prime = 0 ;
  double cnum_prime = 0 ;
  double a1_ln = 0 ;
  double com_a1_ln = 0 ;
  double x1 = 0 ;
  double x2 = 0 ;
  arma::colvec num =  arma::colvec(C+1) ;
  arma::colvec num_prime = arma::colvec(C+1) ;

  for (int j=0; j<=(J-1); j++) {

    tau = 0 ;
    tau_prime = 0 ;
    dentau = 0 ;
    dentau_prime = 0 ;
    sumtau = 0 ;
    dist = 0 ;
    sumdist = 0 ;
    num1 = 0 ;
    num2 = 0 ;
    num1_prime = 0 ;
    num2_prime = 0 ;
    num1_z_prime = 0 ;
    num2_z_prime = 0 ;
    a = 0 ;
    a_prime = 0 ;
    b = 0 ;
    b_prime = 0 ;
    c_prime = 0 ;
    cnum_prime = 0 ;
    a1_ln = 0 ;
    x1 = 0 ;
    x2 = 0 ;

    for (int d=0; d<=(D-1);d++) {
      dist = pow(as_scalar(par.row(d)),2)*pow(as_scalar(Theta(j,d)-par.row(D+d)),2) ;
      sumdist = dist+sumdist ;
    }
    dist = sqrt(sumdist) ;

    for (int w=0; w<=C;w++) {
      x1 = w*dist ;
      x2 = (M-w)*dist ;

      if (w>0) {
        for (int d=0; d<=(D-1);d++) {
          tau = as_scalar(par.row(d)*par.row(w+2*D-1)) ;
          sumtau = tau + sumtau ;
        }
      }

      if (w==as_scalar(Z.row(j))) {
        num1_z_prime = as_scalar(Z.row(j)*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x1)) ;
        num2_z_prime = as_scalar((M-Z.row(j))*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x2)) ;
        a = exp(x1)+exp(x2) ;
        a_prime = (num1_z_prime+num2_z_prime)/dist ;
        c_prime = 0 ;
        if (w>0) {
          for (int t=0; t<=(as_scalar(Z.row(j))-1);t++) {
            cnum_prime = as_scalar(par.row(t+2*D)) ;
            c_prime = c_prime + cnum_prime ;
          }
        }
      }

      dentau = exp(sumtau) ;
      dentau_prime = 0 ;
      if (w>0) {
        for (int t=0; t<=(w-1);t++) {
          tau_prime = as_scalar(par.row(t+2*D)) ;
          dentau_prime = tau_prime + dentau_prime ;
        }
      }
      num1 = exp(x1) ;
      num2 = exp(x2) ;

      num1_prime = as_scalar(w*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x1)) ;
      num2_prime = as_scalar((M-w)*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x2)) ;

      num.row(w) = dentau*(num1+num2) ;
      num_prime.row(w) = (dentau*((num1_prime+num2_prime)/dist) + (num1+num2)*dentau*dentau_prime) ;
    }

    b = sum(num) ;
    b_prime = sum(num_prime) ;
    a1_ln = c_prime + (a_prime/a) - (b_prime/b) ;
    com_a1_ln = com_a1_ln + a1_ln ;

  }

  return(com_a1_ln);

}

/////////////////////////////////////////////////////////////////////

#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

double d1_cpp(
    arma::colvec par,
    arma::mat Theta,
    const int D,
    const int C,
    arma::colvec Z,
    const int dref) {

  const int J = Z.n_rows ;
  const int M = 2*C+1 ;
  double tau = 0 ;
  // double tau_prime = 0 ;
  double dentau = 0 ;
  // double dentau_prime = 0 ;
  double sumtau = 0 ;
  double dist = 0 ;
  double sumdist = 0 ;
  double num1 = 0 ;
  double num2 = 0 ;
  double num1_prime = 0 ;
  double num2_prime = 0 ;
  double num1_z_prime = 0 ;
  double num2_z_prime = 0 ;
  double a = 0 ;
  double a_prime = 0 ;
  double b = 0 ;
  double b_prime = 0 ;
  // double c_prime = 0 ;
  // double cnum_prime = 0 ;
  double d1_ln = 0 ;
  double com_d1_ln = 0 ;
  double x1 = 0 ;
  double x2 = 0 ;
  arma::colvec num =  arma::colvec(C+1) ;
  arma::colvec num_prime = arma::colvec(C+1) ;

  for (int j=0; j<=(J-1); j++) {

    tau = 0 ;
    // tau_prime = 0 ;
    dentau = 0 ;
    //dentau_prime = 0 ;
    sumtau = 0 ;
    dist = 0 ;
    sumdist = 0 ;
    num1 = 0 ;
    num2 = 0 ;
    num1_prime = 0 ;
    num2_prime = 0 ;
    num1_z_prime = 0 ;
    num2_z_prime = 0 ;
    a = 0 ;
    a_prime = 0 ;
    b = 0 ;
    b_prime = 0 ;
    // c_prime = 0 ;
    // cnum_prime = 0 ;
    d1_ln = 0 ;
    x1 = 0 ;
    x2 = 0 ;

    for (int d=0; d<=(D-1);d++) {
      dist = pow(as_scalar(par.row(d)),2)*pow(as_scalar(Theta(j,d)-par.row(D+d)),2) ;
      sumdist = dist+sumdist ;
    }
    dist = sqrt(sumdist) ;

    for (int w=0; w<=C;w++) {
      x1 = w*dist ;
      x2 = (M-w)*dist ;

      if (w>0) {
        for (int d=0; d<=(D-1);d++) {
          tau = as_scalar(par.row(d)*par.row(w+2*D-1)) ;
          sumtau = tau + sumtau ;
        }
      }

      if (w==as_scalar(Z.row(j))) {
        num1_z_prime = as_scalar(Z.row(j)*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1))*exp(x1)) ;
        num2_z_prime = as_scalar((M-Z.row(j))*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1))*exp(x2)) ;
        a = exp(x1)+exp(x2) ;
        a_prime = (num1_z_prime+num2_z_prime)/dist ;
        //         c_prime = 0 ;
        //         if (w>0) {
        //           for (int t=0; t<=(as_scalar(Z.row(j))-1);t++) {
        //             cnum_prime = as_scalar(par.row(t+2*D)) ;
        //             c_prime = c_prime + cnum_prime ;
        //           }
        //         }
      }

      dentau = exp(sumtau) ;
      // dentau_prime = 0 ;
      //       if (w>0) {
      //         for (int t=0; t<=(w-1);t++) {
      //           tau_prime = as_scalar(par.row(t+2*D)) ;
      //           dentau_prime = tau_prime + dentau_prime ;
      //         }
      //       }
      num1 = exp(x1) ;
      num2 = exp(x2) ;

      num1_prime = as_scalar(w*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1))*exp(x1)) ;
      num2_prime = as_scalar((M-w)*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1))*exp(x2)) ;

      num.row(w) = dentau*(num1+num2) ;
      num_prime.row(w) = dentau*(num1_prime+num2_prime)/dist ;
    }

    b = sum(num) ;
    b_prime = sum(num_prime) ;
    d1_ln = (a_prime/a) - (b_prime/b) ;
    com_d1_ln = com_d1_ln + d1_ln ;

  }

  return (com_d1_ln) ;

}

////////////////////////////////////////////////////////////////////////////

#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

double t1_cpp(
    arma::colvec par,
    arma::mat Theta,
    const int D,
    const int C,
    arma::colvec Z,
    const int tauref) {

  const int J = Z.n_rows ;
  const int M = 2*C+1 ;
  double tau = 0 ;
  double tau_prime = 0 ;
  double dentau = 0 ;
  double dentau_prime = 0 ;
  double sumtau = 0 ;
  double dist = 0 ;
  double sumdist = 0 ;
  double num1 = 0 ;
  double num2 = 0 ;
  double b = 0 ;
  double b_prime = 0 ;
  double c_prime = 0 ;
  double cnum_prime = 0 ;
  double t1_ln = 0 ;
  double com_t1_ln = 0 ;
  double x1 = 0 ;
  double x2 = 0 ;
  double U_z = 0 ;
  double U_w = 0 ;
  arma::colvec num =  arma::colvec(C+1) ;
  arma::colvec num_prime = arma::colvec(C+1) ;

  for (int j=0; j<=(J-1); j++) {

    tau = 0 ;
    tau_prime = 0 ;
    dentau = 0 ;
    dentau_prime = 0 ;
    sumtau = 0 ;
    dist = 0 ;
    sumdist = 0 ;
    num1 = 0 ;
    num2 = 0 ;
    b = 0 ;
    b_prime = 0 ;
    c_prime = 0 ;
    cnum_prime = 0 ;
    t1_ln = 0 ;
    x1 = 0 ;
    x2 = 0 ;

    for (int d=0; d<=(D-1);d++) {
      dist = pow(as_scalar(par.row(d)),2)*pow(as_scalar(Theta(j,d)-par.row(D+d)),2) ;
      sumdist = dist+sumdist ;
    }
    dist = sqrt(sumdist) ;

    for (int w=0; w<=C;w++) {
      x1 = w*dist ;
      x2 = (M-w)*dist ;

      if (w>0) {
        for (int d=0; d<=(D-1);d++) {
          tau = as_scalar(par.row(d)*par.row(w+2*D-1)) ;
          sumtau = tau + sumtau ;
        }
      }

      if (w==as_scalar(Z.row(j))) {
        c_prime = 0 ;
        if (w>0) {
          for (int d=0; d<=(D-1);d++) {
            cnum_prime = as_scalar(par.row(d)) ;
            c_prime = c_prime + cnum_prime ;
          } }

        U_z = 0 ;

        if (tauref <= as_scalar(Z.row(j))) {U_z = 1 ; }

        c_prime = U_z*c_prime ;
      }

      dentau = exp(sumtau) ;
      dentau_prime = 0 ;
      if (w>0) {
        for (int d=0; d<=(D-1);d++) {
          tau_prime = as_scalar(par.row(d)) ;
          dentau_prime = tau_prime + dentau_prime ;
        }
      }
      num1 = exp(x1) ;
      num2 = exp(x2) ;

      U_w = 0 ;

      if (tauref <= w) {U_w = 1 ; }

      num.row(w) = dentau*(num1+num2) ;
      num_prime.row(w) = U_w*dentau_prime*dentau*(num1+num2) ;
    }

    b = sum(num) ;
    b_prime = sum(num_prime) ;
    t1_ln = c_prime - (b_prime/b) ;
    com_t1_ln = com_t1_ln + t1_ln ;

  }

  return (com_t1_ln) ;

  //   return List::create(Named("c_prime")=c_prime,Named("a_prime")=a_prime,
  //                       Named("a")=a,Named("b_prime")=b_prime,
  //                       Named("b")=b,Named("dist")=dist,
  //                       Named("x1")=x1,Named("x2")=x2) ;

}

#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

// [[Rcpp::export]]

NumericVector ggum_grad (
    arma::colvec par,
    arma::mat Theta,
    int D,
    int C,
    arma::colvec Z) {

  int dref ;
  int tauref ;
  int ind ;
  int ind2 ;
  //  arma::vec grad =  arma::vec(2*D+C) ;
  NumericVector grad(2*D+C) ;

  for (int d=0; d<=(D-1);d++) {
    // grad.row(d) = a1_cpp(par=par,Theta=Theta,D=D,C=C,Z=Z,dref=(d+1)) ;
    grad(d) = a1_cpp(par=par,Theta=Theta,D=D,C=C,Z=Z,dref=(d+1)) ;
    ind = D+d ;
    grad(ind) = d1_cpp(par=par,Theta=Theta,D=D,C=C,Z=Z,dref=(d+1)) ;
  }

  for (int t=0; t<=(C-1);t++) {
    ind2 = 2*D+t ;
    grad(ind2) = t1_cpp(par=par,Theta=Theta,D=D,C=C,Z=Z,tauref=(t+1)) ;
  }

  // grad2=grad ;

  return(grad) ;

}

///////////////////////////////////////////////////////////////////////

// #include <RcppArmadillo.h>
//
// // [[Rcpp::depends("RcppArmadillo")]]
//
// using namespace Rcpp;
//
// double a2_cpp(
//     arma::colvec par,
//     arma::mat Theta,
//     const int D,
//     const int C,
//     arma::colvec Z,
//     const int dref,
//     const int dref2) {
//
//   const int J = Z.n_rows ;
//   const int M = 2*C+1 ;
//   double tau = 0 ;
//   double tau_prime = 0 ;
//   double dentau = 0 ;
//   double dentau_prime = 0 ;
//   double sumtau = 0 ;
//   double dist = 0 ;
//   double dist_prime = 0 ;
//   double sumdist = 0 ;
//   double num1 = 0 ;
//   double num2 = 0 ;
//   double num1_prime = 0 ;
//   double num2_prime = 0 ;
//   double num1_z_prime = 0 ;
//   double num2_z_prime = 0 ;
//   double num1_prime2 = 0 ;
//   double a = 0 ;
//   double a_prime = 0 ;
//   double b = 0 ;
//   double b_prime = 0 ;
//   double c_prime = 0 ;
//   double cnum_prime = 0 ;
//   double e = 0 ;
//   double e_prime = 0 ;
//   double f = 0 ;
//   double f_prime = 0 ;
//   double g = 0 ;
//   double g_prime = 0 ;
//   double h = 0 ;
//   double h_prime = 0 ;
//   double a2_ln = 0 ;
//   double com_a2_ln = 0 ;
//   double x1 = 0 ;
//   double x2 = 0 ;
//   arma::colvec num =  arma::colvec(C+1) ;
//   arma::colvec num_prime = arma::colvec(C+1) ;
//   arma::colvec num_prime2 = arma::colvec(C+1) ;
//   arma::colvec num_prime_test = arma::colvec(C+1) ;
//
//   for (int j=0; j<=(J-1); j++) {
//
//     tau = 0 ;
//     tau_prime = 0 ;
//     dentau = 0 ;
//     dentau_prime = 0 ;
//     sumtau = 0 ;
//     dist = 0 ;
//     dist_prime = 0 ;
//     sumdist = 0 ;
//     num1 = 0 ;
//     num2 = 0 ;
//     num1_prime = 0 ;
//     num2_prime = 0 ;
//     num1_z_prime = 0 ;
//     num2_z_prime = 0 ;
//     num1_prime2 = 0 ;
//     a = 0 ;
//     a_prime = 0 ;
//     b = 0 ;
//     b_prime = 0 ;
//     c_prime = 0 ;
//     cnum_prime = 0 ;
//     e = 0 ;
//     e_prime = 0 ;
//     f = 0 ;
//     f_prime = 0 ;
//     g = 0 ;
//     g_prime = 0 ;
//     h = 0 ;
//     h_prime = 0 ;
//     a2_ln = 0 ;
//     x1 = 0 ;
//     x2 = 0 ;
//
//     for (int d=0; d<=(D-1);d++) {
//       dist = pow(as_scalar(par.row(d)),2)*pow(as_scalar(Theta(j,d)-par.row(D+d)),2) ;
//       sumdist = dist+sumdist ;
//     }
//     dist = sqrt(sumdist) ;
//     dist_prime = as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist ;
//
//     for (int w=0; w<=C;w++) {
//       x1 = w*dist ;
//       x2 = (M-w)*dist ;
//
//       if (w>0) {
//         for (int d=0; d<=(D-1);d++) {
//           tau = as_scalar(par.row(d)*par.row(w+2*D-1)) ;
//           sumtau = tau + sumtau ;
//         }
//       }
//
//       if (w==as_scalar(Z.row(j))) {
//         num1_z_prime = as_scalar(Z.row(j)*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x1)) ;
//         num2_z_prime = as_scalar((M-Z.row(j))*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x2)) ;
//
//         a = exp(x1)+exp(x2) ;
//         a_prime = (num1_z_prime+num2_z_prime)/dist ;
//
//         e = a_prime*dist ;
//         e_prime = as_scalar(Z.row(j)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1) +
//           (M-Z.row(j))*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2))+
//           ((pow(as_scalar(Z.row(j)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x1)) +
//           (pow(as_scalar((M-Z.row(j))*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x2)))/dist ;
//
//         f = a*dist ;
//         f_prime = a*(as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist) +
//           dist*a_prime ;
//         dist_prime = as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist ;
//
//         c_prime = 0 ;
//         if (w>0) {
//           for (int t=0; t<=(as_scalar(Z.row(j))-1);t++) {
//             cnum_prime = as_scalar(par.row(t+2*D)) ;
//             c_prime = c_prime + cnum_prime ;
//           }
//         }
//       }
//
//       dentau = exp(sumtau) ;
//       dentau_prime = 0 ;
//       if (w>0) {
//         for (int t=0; t<=(w-1);t++) {
//           tau_prime = as_scalar(par.row(t+2*D)) ;
//           dentau_prime = tau_prime + dentau_prime ;
//         }
//       }
//       num1 = exp(x1) ;
//       num2 = exp(x2) ;
//
//       num1_prime = as_scalar(w*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x1)) ;
//       num2_prime = as_scalar((M-w)*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x2)) ;
//
//       num1_prime2 = (w*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1) +
//         (M-w)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2))+
//         ((pow(as_scalar(w*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x1)) +
//         (pow(as_scalar((M-w)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x2)))/dist ;
//
//       num.row(w) = dentau*(num1+num2) ;
//       num_prime.row(w) = (dentau*((num1_prime+num2_prime)/dist) + (num1+num2)*dentau*dentau_prime) ;
//       num_prime2.row(w) = dentau*num1_prime2+dentau*dentau_prime*(num1_prime+num2_prime)+
//         (num1+num2)*(dentau*dentau_prime*dist_prime+dentau*dentau_prime*dentau_prime*dist)+
//         ((num1_prime+num2_prime)/dist)*(dentau*dentau_prime*dist) ;
//     }
//
//     b = sum(num) ;
//     b_prime = sum(num_prime) ;
//
//     g = b_prime*dist ;
//     g_prime = sum(num_prime2) ;
//
//     h = b*dist ;
//     h_prime = b*dist_prime + b_prime*dist ;
//
//     a2_ln = (f*e_prime-e*f_prime)/pow(f,2) - (h*g_prime-g*h_prime)/pow(h,2) ;
//     com_a2_ln = com_a2_ln + a2_ln ;
//
//   }
//
//   return (com_a2_ln) ;
//
//   //   return List::create(Named("e")=e,Named("e_prime")=e_prime,
//   //                       Named("f")=f,Named("f_prime")=f_prime,
//   //                       Named("g")=g,Named("g_prime")=g_prime,
//   //                       Named("h")=h,Named("h_prime")=h_prime) ;
//
// }
//
//
//
//
// ///////////////////////////////////////////////
//
// #include <RcppArmadillo.h>
//
// // [[Rcpp::depends("RcppArmadillo")]]
//
// using namespace Rcpp;
//
// double d2_cpp(
//     arma::colvec par,
//     arma::mat Theta,
//     const int D,
//     const int C,
//     arma::colvec Z,
//     const int dref,
//     const int dref2) {
//
//   const int J = Z.n_rows ;
//   const int M = 2*C+1 ;
//   double tau = 0 ;
//   // double tau_prime = 0 ;
//   double dentau = 0 ;
//   // double dentau_prime = 0 ;
//   double sumtau = 0 ;
//   double dist = 0 ;
//   double dist_prime = 0 ;
//   double sumdist = 0 ;
//   double num1 = 0 ;
//   double num2 = 0 ;
//   double num1_prime = 0 ;
//   double num2_prime = 0 ;
//   double num1_z_prime = 0 ;
//   double num2_z_prime = 0 ;
//   double num1_prime2 = 0 ;
//   double a = 0 ;
//   double a_prime = 0 ;
//   double b = 0 ;
//   double b_prime = 0 ;
//   // double c_prime = 0 ;
//   // double cnum_prime = 0 ;
//   double e = 0 ;
//   double e_prime = 0 ;
//   double f = 0 ;
//   double f_prime = 0 ;
//   double g = 0 ;
//   double g_prime = 0 ;
//   double h = 0 ;
//   double h_prime = 0 ;
//   double d2_ln = 0 ;
//   double com_d2_ln = 0 ;
//   double x1 = 0 ;
//   double x2 = 0 ;
//   // double test1 = 0 ;
//   //  double test2 = 0 ;
//   arma::colvec num =  arma::colvec(C+1) ;
//   arma::colvec num_prime = arma::colvec(C+1) ;
//   arma::colvec num_prime2 = arma::colvec(C+1) ;
//   arma::colvec num_prime_test = arma::colvec(C+1) ;
//
//   for (int j=0; j<=(J-1); j++) {
//
//     tau = 0 ;
//     // tau_prime = 0 ;
//     dentau = 0 ;
//     // dentau_prime = 0 ;
//     sumtau = 0 ;
//     dist = 0 ;
//     dist_prime = 0 ;
//     sumdist = 0 ;
//     num1 = 0 ;
//     num2 = 0 ;
//     num1_prime = 0 ;
//     num2_prime = 0 ;
//     num1_z_prime = 0 ;
//     num2_z_prime = 0 ;
//     num1_prime2 = 0 ;
//     a = 0 ;
//     a_prime = 0 ;
//     b = 0 ;
//     b_prime = 0 ;
//     // c_prime = 0 ;
//     // cnum_prime = 0 ;
//     e = 0 ;
//     e_prime = 0 ;
//     f = 0 ;
//     f_prime = 0 ;
//     g = 0 ;
//     g_prime = 0 ;
//     h = 0 ;
//     h_prime = 0 ;
//     d2_ln = 0 ;
//     x1 = 0 ;
//     x2 = 0 ;
//     //  test1 = 0 ;
//     //  test2 = 0 ;
//
//     for (int d=0; d<=(D-1);d++) {
//       dist = pow(as_scalar(par.row(d)),2)*pow(as_scalar(Theta(j,d)-par.row(D+d)),2) ;
//       sumdist = dist+sumdist ;
//     }
//     dist = sqrt(sumdist) ;
//     //dist_prime = as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist ;
//
//     for (int w=0; w<=C;w++) {
//       x1 = w*dist ;
//       x2 = (M-w)*dist ;
//
//       if (w>0) {
//         for (int d=0; d<=(D-1);d++) {
//           tau = as_scalar(par.row(d)*par.row(w+2*D-1)) ;
//           sumtau = tau + sumtau ;
//         }
//       }
//
//       if (w==as_scalar(Z.row(j))) {
//         num1_z_prime = as_scalar((Z.row(j)*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1)))*exp(x1)) ;
//         num2_z_prime = as_scalar(((M-Z.row(j))*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1)))*exp(x2)) ;
//
//         a = exp(x1)+exp(x2) ;
//         a_prime = (num1_z_prime+num2_z_prime)/dist ;
//
//         e = a_prime*dist ;
//         //      pow1 = as_scalar(pow((Z.row(j)*pow(par.row(dref2-1),2)*(par.row(dref2-1+D)-Theta(j,dref2-1))),2)*exp(x1)) ;
//         //      pow2 = as_scalar((pow((M-Z.row(j))*as_scalar(pow(par.row(dref2-1),2))*(par.row(dref2-1+D)-Theta(j,dref2-1))),2)*exp(x1)))/dist ;
//
//         e_prime = as_scalar(Z.row(j)*pow(par.row(dref2-1),2)*exp(x1) +
//           (M-Z.row(j))*pow(par.row(dref2-1),2)*exp(x2)) +
//           (as_scalar(pow((Z.row(j)*pow(par.row(dref2-1),2)*(par.row(dref2-1+D)-Theta(j,dref2-1))),2)*exp(x1)) +
//           as_scalar(pow(((M-Z.row(j))*as_scalar(pow(par.row(dref2-1),2))*(par.row(dref2-1+D)-Theta(j,dref2-1))),2)*exp(x2)))/dist ;
//
//         //      test1 = as_scalar((Z.row(j)*pow(as_scalar(par.row(dref2-1)),2)*pow(as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)),2))*exp(x1)) ;
//         //      test2 = as_scalar(((M-Z.row(j))*pow(as_scalar(par.row(dref2-1)),2)*pow(as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)),2))*exp(x2)) ;
//
//
//
//         f = a*dist ;
//         f_prime = a*((pow(as_scalar(par.row(dref2-1)),2)*as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)))/dist) +
//           dist*a_prime ;
//         dist_prime = (pow(as_scalar(par.row(dref2-1)),2)*as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)))/dist ;
//
//         //       c_prime = 0 ;
//         //       if (w>0) {
//         //         for (int t=0; t<=(as_scalar(Z.row(j))-1);t++) {
//         //           cnum_prime = as_scalar(par.row(t+2*D)) ;
//         //           c_prime = c_prime + cnum_prime ;
//         //         }
//         //       }
//       }
//
//       dentau = exp(sumtau) ;
//       //     dentau_prime = 0 ;
//       //       if (w>0) {
//       //         for (int t=0; t<=(w-1);t++) {
//       //           tau_prime = as_scalar(par.row(t+2*D)) ;
//       //           dentau_prime = tau_prime + dentau_prime ;
//       //         }
//       //       }
//       num1 = exp(x1) ;
//       num2 = exp(x2) ;
//
//       num1_prime = as_scalar((w*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1)))*exp(x1)) ;
//       num2_prime = as_scalar(((M-w)*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1)))*exp(x2)) ;
//
//       num1_prime2 = as_scalar(w*pow(par.row(dref2-1),2)*exp(x1) +
//         (M-w)*pow(par.row(dref2-1),2)*exp(x2)) +
//         (as_scalar(pow((w*pow(par.row(dref2-1),2)*(par.row(dref2-1+D)-Theta(j,dref2-1))),2)*exp(x1)) +
//         as_scalar(pow(((M-w)*as_scalar(pow(par.row(dref2-1),2))*(par.row(dref2-1+D)-Theta(j,dref2-1))),2)*exp(x2)))/dist ;
//
//       num.row(w) = dentau*(num1+num2) ;
//       num_prime.row(w) = dentau*((num1_prime+num2_prime)/dist) ;
//       num_prime2.row(w) = dentau*num1_prime2 ;
//     }
//
//     b = sum(num) ;
//     b_prime = sum(num_prime) ;
//
//     g = b_prime*dist ;
//     g_prime = sum(num_prime2) ;
//
//     h = b*dist ;
//     h_prime = b*dist_prime + b_prime*dist ;
//
//     d2_ln = (f*e_prime-e*f_prime)/pow(f,2) - (h*g_prime-g*h_prime)/pow(h,2) ;
//     com_d2_ln = com_d2_ln + d2_ln ;
//
//   }
//
//   return (com_d2_ln) ;
//
//   //   return List::create(Named("e")=e,Named("e_prime")=e_prime,
//   //                       Named("f")=f,Named("f_prime")=f_prime,
//   //                       Named("g")=g,Named("g_prime")=g_prime,
//   //                       Named("h")=h,Named("h_prime")=h_prime) ;
//
//
// }
//
//
//
//
// ////////////////////////////////////////
//
// #include <RcppArmadillo.h>
//
// // [[Rcpp::depends("RcppArmadillo")]]
//
// using namespace Rcpp;
//
// double t2_cpp(
//     arma::colvec par,
//     arma::mat Theta,
//     const int D,
//     const int C,
//     arma::colvec Z,
//     const int tauref) {
//
//   const int J = Z.n_rows ;
//   const int M = 2*C+1 ;
//   double tau = 0 ;
//   double tau_prime = 0 ;
//   double dentau = 0 ;
//   double dentau_prime = 0 ;
//   double sumtau = 0 ;
//   double dist = 0 ;
//   double sumdist = 0 ;
//   double num1 = 0 ;
//   double num2 = 0 ;
//   double b = 0 ;
//   double b_prime = 0 ;
//   double c_prime = 0 ;
//   double cnum_prime = 0 ;
//   double g = 0 ;
//   double g_prime = 0 ;
//   double h = 0 ;
//   double h_prime = 0 ;
//   double t1_ln = 0 ;
//   double com_t1_ln = 0 ;
//   double t2_ln = 0 ;
//   double com_t2_ln = 0 ;
//   double x1 = 0 ;
//   double x2 = 0 ;
//   double U_z = 0 ;
//   double U_w = 0 ;
//   arma::colvec num =  arma::colvec(C+1) ;
//   arma::colvec num_prime = arma::colvec(C+1) ;
//   arma::colvec num_prime2 = arma::colvec(C+1) ;
//
//   for (int j=0; j<=(J-1); j++) {
//
//     tau = 0 ;
//     tau_prime = 0 ;
//     dentau = 0 ;
//     dentau_prime = 0 ;
//     sumtau = 0 ;
//     dist = 0 ;
//     sumdist = 0 ;
//     num1 = 0 ;
//     num2 = 0 ;
//     b = 0 ;
//     b_prime = 0 ;
//     g = 0 ;
//     g_prime = 0 ;
//     h = 0 ;
//     h_prime = 0 ;
//     c_prime = 0 ;
//     cnum_prime = 0 ;
//     t1_ln = 0 ;
//     t2_ln = 0 ;
//     x1 = 0 ;
//     x2 = 0 ;
//
//     for (int d=0; d<=(D-1);d++) {
//       dist = pow(as_scalar(par.row(d)),2)*pow(as_scalar(Theta(j,d)-par.row(D+d)),2) ;
//       sumdist = dist+sumdist ;
//     }
//     dist = sqrt(sumdist) ;
//
//     for (int w=0; w<=C;w++) {
//       x1 = w*dist ;
//       x2 = (M-w)*dist ;
//
//       if (w>0) {
//         for (int d=0; d<=(D-1);d++) {
//           tau = as_scalar(par.row(d)*par.row(w+2*D-1)) ;
//           sumtau = tau + sumtau ;
//         }
//       }
//
//       if (w==as_scalar(Z.row(j))) {
//         c_prime = 0 ;
//         if (w>0) {
//           for (int d=0; d<=(D-1);d++) {
//             cnum_prime = as_scalar(par.row(d)) ;
//             c_prime = c_prime + cnum_prime ;
//           } }
//
//         U_z = 0 ;
//
//         if (tauref <= as_scalar(Z.row(j))) {U_z = 1 ; }
//
//         c_prime = U_z*c_prime ;
//       }
//
//       dentau = exp(sumtau) ;
//       dentau_prime = 0 ;
//       if (w>0) {
//         for (int d=0; d<=(D-1);d++) {
//           tau_prime = as_scalar(par.row(d)) ;
//           dentau_prime = tau_prime + dentau_prime ;
//         }
//       }
//       num1 = exp(x1) ;
//       num2 = exp(x2) ;
//
//       U_w = 0 ;
//
//       if (tauref <= w) {U_w = 1 ; }
//
//       num.row(w) = dentau*(num1+num2) ;
//       num_prime.row(w) = U_w*dentau_prime*dentau*(num1+num2) ;
//       num_prime2.row(w) = U_w*pow(dentau_prime,2)*dentau*(num1+num2) ;
//     }
//
//     b = sum(num) ;
//     b_prime = sum(num_prime) ;
//     t1_ln = c_prime - (b_prime/b) ;
//     com_t1_ln = com_t1_ln + t1_ln ;
//
//     g = b_prime ;
//     g_prime = sum(num_prime2) ;
//     h = b ;
//     h_prime = b_prime ;
//
//     t2_ln = -(h*g_prime - g*h_prime)/pow(h,2) ;
//     com_t2_ln = com_t2_ln + t2_ln ;
//
//   }
//
//   return (com_t2_ln) ;
//
//   //   return List::create(Named("c_prime")=c_prime,Named("a_prime")=a_prime,
//   //                       Named("a")=a,Named("b_prime")=b_prime,
//   //                       Named("b")=b,Named("dist")=dist,
//   //                       Named("x1")=x1,Named("x2")=x2) ;
//
// }
//
//
//
//
// //////////////////////////////////////////
//
// #include <RcppArmadillo.h>
//
// // [[Rcpp::depends("RcppArmadillo")]]
//
// using namespace Rcpp;
//
// double a2mix_cpp(
//     arma::colvec par,
//     arma::mat Theta,
//     const int D,
//     const int C,
//     arma::colvec Z,
//     const int dref,
//     const int dref2) {
//
//   const int J = Z.n_rows ;
//   const int M = 2*C+1 ;
//   double tau = 0 ;
//   double tau_prime = 0 ;
//   double dentau = 0 ;
//   double dentau_prime = 0 ;
//   double sumtau = 0 ;
//   double dist = 0 ;
//   double dist_prime = 0 ;
//   double dist_prime_mix = 0 ;
//   double sumdist = 0 ;
//   double num1 = 0 ;
//   double num2 = 0 ;
//   double num1_prime = 0 ;
//   double num2_prime = 0 ;
//   double num1_prime_mix = 0 ;
//   double num2_prime_mix = 0 ;
//   double num1_z_prime = 0 ;
//   double num2_z_prime = 0 ;
//   double num1_prime2 = 0 ;
//   // double num1_prime2_mix = 0 ;
//   double a = 0 ;
//   double a_prime = 0 ;
//   double a_prime_mix = 0 ;
//   double b = 0 ;
//   double b_prime = 0 ;
//   double b_prime_mix = 0 ;
//   double c_prime = 0 ;
//   double cnum_prime = 0 ;
//   double e = 0 ;
//   // double e_prime = 0 ;
//   double e_prime_mix = 0 ;
//   double f = 0 ;
//   // double f_prime = 0 ;
//   double f_prime_mix = 0 ;
//   double g = 0 ;
//   // double g_prime = 0 ;
//   double g_prime_mix = 0 ;
//   double h = 0 ;
//   // double h_prime = 0 ;
//   double h_prime_mix = 0 ;
//   // double a2_ln = 0 ;
//   double a2mix_ln = 0 ;
//   // double com_a2_ln = 0 ;
//   double com_a2mix_ln = 0 ;
//   double x1 = 0 ;
//   double x2 = 0 ;
//   arma::colvec num =  arma::colvec(C+1) ;
//   arma::colvec num_prime = arma::colvec(C+1) ;
//   arma::colvec num_prime_mix = arma::colvec(C+1) ;
//   arma::colvec num_prime2 = arma::colvec(C+1) ;
//   arma::colvec num_prime2_mix = arma::colvec(C+1) ;
//   arma::colvec num_prime_test = arma::colvec(C+1) ;
//
//   for (int j=0; j<=(J-1); j++) {
//
//     tau = 0 ;
//     tau_prime = 0 ;
//     dentau = 0 ;
//     dentau_prime = 0 ;
//     sumtau = 0 ;
//     dist = 0 ;
//     dist_prime = 0 ;
//     dist_prime_mix = 0 ;
//     sumdist = 0 ;
//     num1 = 0 ;
//     num2 = 0 ;
//     num1_prime = 0 ;
//     num2_prime = 0 ;
//     num1_prime_mix = 0 ;
//     num2_prime_mix = 0 ;
//     num1_z_prime = 0 ;
//     num2_z_prime = 0 ;
//     num1_prime2 = 0 ;
//     // num1_prime2_mix = 0 ;
//     a = 0 ;
//     a_prime = 0 ;
//     a_prime_mix = 0 ;
//     b = 0 ;
//     b_prime = 0 ;
//     b_prime_mix = 0 ;
//     c_prime = 0 ;
//     cnum_prime = 0 ;
//     e = 0 ;
//     // e_prime = 0 ;
//     e_prime_mix = 0 ;
//     f = 0 ;
//     // f_prime = 0 ;
//     f_prime_mix = 0 ;
//     g = 0 ;
//     // g_prime = 0 ;
//     g_prime_mix = 0 ;
//     h = 0 ;
//     // h_prime = 0 ;
//     h_prime_mix = 0 ;
//     // a2_ln = 0 ;
//     a2mix_ln = 0 ;
//     x1 = 0 ;
//     x2 = 0 ;
//
//     for (int d=0; d<=(D-1);d++) {
//       dist = pow(as_scalar(par.row(d)),2)*pow(as_scalar(Theta(j,d)-par.row(D+d)),2) ;
//       sumdist = dist+sumdist ;
//     }
//     dist = sqrt(sumdist) ;
//     dist_prime = as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist ;
//     dist_prime_mix = dist_prime ;
//
//
//     for (int w=0; w<=C;w++) {
//       x1 = w*dist ;
//       x2 = (M-w)*dist ;
//
//       if (w>0) {
//         for (int d=0; d<=(D-1);d++) {
//           tau = as_scalar(par.row(d)*par.row(w+2*D-1)) ;
//           sumtau = tau + sumtau ;
//         }
//       }
//
//       if (w==as_scalar(Z.row(j))) {
//         num1_z_prime = as_scalar(Z.row(j)*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x1)) ;
//         num2_z_prime = as_scalar((M-Z.row(j))*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x2)) ;
//
//         a = exp(x1)+exp(x2) ;
//         a_prime = (num1_z_prime+num2_z_prime)/dist ;
//         a_prime_mix = (as_scalar(Z.row(j)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*exp(x1) +
//           as_scalar((M-Z.row(j))*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*exp(x2))/dist ;
//
//         e = a_prime*dist ;
//         //       e_prime = as_scalar(Z.row(j)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1) +
//         //         (M-Z.row(j))*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2))+
//         //         ((pow(as_scalar(Z.row(j)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x1)) +
//         //         (pow(as_scalar((M-Z.row(j))*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x2)))/dist ;
//
//         e_prime_mix = (as_scalar(num1_z_prime*(Z.row(j)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))) +
//           as_scalar(num2_z_prime*((M-Z.row(j))*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))))/dist ;
//
//         f = a*dist ;
//         //       f_prime = a*(as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist) +
//         //                   dist*a_prime ;
//         f_prime_mix = a*dist_prime_mix + dist*a_prime_mix ;
//
//         c_prime = 0 ;
//         if (w>0) {
//           for (int t=0; t<=(as_scalar(Z.row(j))-1);t++) {
//             cnum_prime = as_scalar(par.row(t+2*D)) ;
//             c_prime = c_prime + cnum_prime ;
//           }
//         }
//       }
//
//       dentau = exp(sumtau) ;
//       dentau_prime = 0 ;
//       if (w>0) {
//         for (int t=0; t<=(w-1);t++) {
//           tau_prime = as_scalar(par.row(t+2*D)) ;
//           dentau_prime = tau_prime + dentau_prime ;
//         }
//       }
//       num1 = exp(x1) ;
//       num2 = exp(x2) ;
//
//       num1_prime = as_scalar(w*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x1)) ;
//       num2_prime = as_scalar((M-w)*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x2)) ;
//
//       num1_prime_mix = as_scalar(w*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1)) ;
//       num2_prime_mix = as_scalar((M-w)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2)) ;
//
//       num1_prime2 = (w*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1) +
//         (M-w)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2))+
//         ((pow(as_scalar(w*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x1)) +
//         (pow(as_scalar((M-w)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x2)))/dist ;
//
//       num.row(w) = dentau*(num1+num2) ;
//       num_prime.row(w) = (dentau*((num1_prime+num2_prime)/dist) + (num1+num2)*dentau*dentau_prime) ;
//
//       num_prime_mix.row(w) = (dentau*((num1_prime_mix+num2_prime_mix)/dist) + (num1+num2)*dentau*dentau_prime) ;
//
//       num_prime2.row(w) = dentau*num1_prime2+dentau*dentau_prime*(num1_prime+num2_prime)+
//         (num1+num2)*(dentau*dentau_prime*dist_prime+dentau*dentau_prime*dentau_prime*dist)+
//         ((num1_prime+num2_prime)/dist)*(dentau*dentau_prime*dist) ;
//
//       num_prime2_mix.row(w) = dentau*(num1_prime*(w*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)) +
//         num2_prime*((M-w)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)))/dist +
//         dentau_prime*dentau*(num1_prime+num2_prime)+
//         (num1+num2)*(dentau*dentau_prime*dist_prime_mix+dentau*dentau_prime*dentau_prime*dist) +
//         (((w*(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)))*num1 +
//         ((M-w)*(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)))*num2)/dist)*(dentau*dentau_prime*dist) ;
//
//     }
//
//     b = sum(num) ;
//     b_prime = sum(num_prime) ;
//     b_prime_mix = sum(num_prime_mix) ;
//
//     g = b_prime*dist ;
//     // g_prime = sum(num_prime2) ;
//     g_prime_mix = sum(num_prime2_mix) ;
//
//     h = b*dist ;
//     // h_prime = b*dist_prime + b_prime*dist ;
//     h_prime_mix = b*dist_prime_mix + b_prime_mix*dist ;
//
//     // a2_ln = (f*e_prime-e*f_prime)/pow(f,2) - (h*g_prime-g*h_prime)/pow(h,2) ;
//     a2mix_ln = (f*e_prime_mix-e*f_prime_mix)/pow(f,2) - (h*g_prime_mix-g*h_prime_mix)/pow(h,2) ;
//     // com_a2_ln = com_a2_ln + a2_ln ;
//     com_a2mix_ln = com_a2mix_ln + a2mix_ln ;
//
//   }
//
//   return (com_a2mix_ln) ;
//
//   //   return List::create(Named("e")=e,Named("e_prime_mix")=e_prime_mix,
//   //                       Named("f")=f,Named("f_prime_mix")=f_prime_mix,
//   //                       Named("g")=g,Named("g_prime_mix")=g_prime_mix,
//   //                       Named("h")=h,Named("h_prime_mix")=h_prime_mix) ;
//
//   //   return List::create(Named("a")=a,Named("dist_prime_mix")=dist_prime_mix,
//   //                       Named("dist")=dist,Named("a_prime_mix")=a_prime_mix) ;
//
// }
//
//
//
//
// //////////////////////////////////////////
//
// #include <RcppArmadillo.h>
//
// // [[Rcpp::depends("RcppArmadillo")]]
//
// using namespace Rcpp;
//
// double d2mix_cpp(
//     arma::colvec par,
//     arma::mat Theta,
//     const int D,
//     const int C,
//     arma::colvec Z,
//     const int dref,
//     const int dref2) {
//
//   const int J = Z.n_rows ;
//   const int M = 2*C+1 ;
//   double tau = 0 ;
//   // double tau_prime = 0 ;
//   double dentau = 0 ;
//   // double dentau_prime = 0 ;
//   double sumtau = 0 ;
//   double dist = 0 ;
//   double dist_prime = 0 ;
//   double sumdist = 0 ;
//   double num1 = 0 ;
//   double num2 = 0 ;
//   double num1_prime = 0 ;
//   double num2_prime = 0 ;
//   double num1_z_prime = 0 ;
//   double num2_z_prime = 0 ;
//   //   double num1_prime2 = 0 ;
//   double num1_prime_mix = 0 ;
//   double num2_prime_mix = 0 ;
//   // double num1_z_prime_mix = 0 ;
//   // double num2_z_prime_mix = 0 ;
//   double num1_prime2_mix = 0 ;
//   double a = 0 ;
//   double a_prime = 0 ;
//   double b = 0 ;
//   double b_prime = 0 ;
//   double b_prime_mix = 0 ;
//   // double c_prime = 0 ;
//   // double cnum_prime = 0 ;
//   double e = 0 ;
//   // double e_prime = 0 ;
//   double e_prime_mix = 0 ;
//   double f = 0 ;
//   // double f_prime = 0 ;
//   double f_prime_mix = 0 ;
//   double g = 0 ;
//   // double g_prime = 0 ;
//   double g_prime_mix = 0 ;
//   double h = 0 ;
//   // double h_prime = 0 ;
//   double h_prime_mix = 0 ;
//   // double d2_ln = 0 ;
//   double d2mix_ln = 0 ;
//   // double com_d2_ln = 0 ;
//   double com_d2mix_ln = 0 ;
//   double x1 = 0 ;
//   double x2 = 0 ;
//   arma::colvec num =  arma::colvec(C+1) ;
//   arma::colvec num_prime = arma::colvec(C+1) ;
//   arma::colvec num_prime_mix = arma::colvec(C+1) ;
//   arma::colvec num_prime2 = arma::colvec(C+1) ;
//   arma::colvec num_prime2_mix = arma::colvec(C+1) ;
//   // arma::colvec num_prime_test = arma::colvec(C+1) ;
//
//   for (int j=0; j<=(J-1); j++) {
//
//     tau = 0 ;
//     // tau_prime = 0 ;
//     dentau = 0 ;
//     // dentau_prime = 0 ;
//     sumtau = 0 ;
//     dist = 0 ;
//     dist_prime = 0 ;
//     sumdist = 0 ;
//     num1 = 0 ;
//     num2 = 0 ;
//     num1_prime = 0 ;
//     num2_prime = 0 ;
//     num1_z_prime = 0 ;
//     num2_z_prime = 0 ;
//     //     num1_prime2 = 0 ;
//     num1_prime_mix = 0 ;
//     num2_prime_mix = 0 ;
//     // num1_z_prime_mix = 0 ;
//     // num2_z_prime_mix = 0 ;
//     num1_prime2_mix = 0 ;
//     a = 0 ;
//     a_prime = 0 ;
//     b = 0 ;
//     b_prime = 0 ;
//     b_prime_mix = 0 ;
//     // c_prime = 0 ;
//     // cnum_prime = 0 ;
//     e = 0 ;
//     // e_prime = 0 ;
//     e_prime_mix = 0 ;
//     f = 0 ;
//     // f_prime = 0 ;
//     f_prime_mix = 0 ;
//     g = 0 ;
//     // g_prime = 0 ;
//     g_prime_mix = 0 ;
//     h = 0 ;
//     // h_prime = 0 ;
//     h_prime_mix = 0 ;
//     // d2_ln = 0 ;
//     d2mix_ln = 0 ;
//     x1 = 0 ;
//     x2 = 0 ;
//
//     for (int d=0; d<=(D-1);d++) {
//       dist = pow(as_scalar(par.row(d)),2)*pow(as_scalar(Theta(j,d)-par.row(D+d)),2) ;
//       sumdist = dist+sumdist ;
//     }
//     dist = sqrt(sumdist) ;
//     //dist_prime = as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist ;
//
//     for (int w=0; w<=C;w++) {
//       x1 = w*dist ;
//       x2 = (M-w)*dist ;
//
//       if (w>0) {
//         for (int d=0; d<=(D-1);d++) {
//           tau = as_scalar(par.row(d)*par.row(w+2*D-1)) ;
//           sumtau = tau + sumtau ;
//         }
//       }
//
//       if (w==as_scalar(Z.row(j))) {
//         num1_z_prime = as_scalar((Z.row(j)*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1)))*exp(x1)) ;
//         num2_z_prime = as_scalar(((M-Z.row(j))*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1)))*exp(x2)) ;
//
//         a = exp(x1)+exp(x2) ;
//         a_prime = (num1_z_prime+num2_z_prime)/dist ;
//
//         e = a_prime*dist ;
//         //       e_prime = as_scalar(Z.row(j)*pow(par.row(dref2-1),2)*exp(x1) +
//         //         (M-Z.row(j))*pow(par.row(dref2-1),2)*exp(x2))+
//         //         as_scalar((((Z.row(j)*pow(as_scalar(par.row(dref2-1)),2)*pow(as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)),2))*exp(x1)) +
//         //         ((M-Z.row(j))*pow(as_scalar(par.row(dref2-1)),2)*pow(as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)),2))*exp(x2)))/dist ;
//
//         e_prime_mix = (num1_z_prime*as_scalar(Z.row(j)*pow(par.row(dref2-1),2)*(par.row(dref2-1+D)-Theta(j,dref2-1))) +
//           num2_z_prime*as_scalar((M-Z.row(j))*pow(par.row(dref2-1),2)*(par.row(dref2-1+D)-Theta(j,dref2-1))))/dist ;
//
//         f = a*dist ;
//         //       f_prime = a*((pow(as_scalar(par.row(dref2-1)),2)*pow(as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)),2))/dist) +
//         //                   dist*a_prime ;
//
//         f_prime_mix = a*((pow(as_scalar(par.row(dref2-1)),2)*as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)))/dist) +
//           as_scalar((Z.row(j)*pow(par.row(dref2-1),2)*(par.row(dref2-1+D)-Theta(j,dref2-1)))*exp(x1)) +
//           as_scalar(((M-Z.row(j))*pow(par.row(dref2-1),2)*(par.row(dref2-1+D)-Theta(j,dref2-1)))*exp(x2)) ;
//
//         dist_prime = (pow(as_scalar(par.row(dref2-1)),2)*as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)))/dist ;
//
//         //       c_prime = 0 ;
//         //       if (w>0) {
//         //         for (int t=0; t<=(as_scalar(Z.row(j))-1);t++) {
//         //           cnum_prime = as_scalar(par.row(t+2*D)) ;
//         //           c_prime = c_prime + cnum_prime ;
//         //         }
//         //       }
//       }
//
//       dentau = exp(sumtau) ;
//       //     dentau_prime = 0 ;
//       //       if (w>0) {
//       //         for (int t=0; t<=(w-1);t++) {
//       //           tau_prime = as_scalar(par.row(t+2*D)) ;
//       //           dentau_prime = tau_prime + dentau_prime ;
//       //         }
//       //       }
//       num1 = exp(x1) ;
//       num2 = exp(x2) ;
//
//       num1_prime = as_scalar((w*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1)))*exp(x1)) ;
//       num2_prime = as_scalar(((M-w)*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1)))*exp(x2)) ;
//
//       num1_prime_mix = as_scalar((w*pow(par.row(dref2-1),2)*(par.row(dref2-1+D)-Theta(j,dref2-1)))*exp(x1)) ;
//       num2_prime_mix = as_scalar(((M-w)*pow(par.row(dref2-1),2)*(par.row(dref2-1+D)-Theta(j,dref2-1)))*exp(x2)) ;
//
//       //       num1_prime2 = w*pow(as_scalar(par.row(dref2-1)),2)*exp(x1) +
//       //         (M-w)*pow(as_scalar(par.row(dref2-1)),2)*exp(x2) +
//       //         ((w*pow(as_scalar(par.row(dref2-1)),2)*pow(as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)),2)*exp(x1)) +
//       //         ((M-w)*pow(as_scalar(par.row(dref2-1)),2)*pow(as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)),2)*exp(x2)))/dist ;
//
//       num1_prime2_mix = (num1_prime*(w*pow(as_scalar(par.row(dref2-1)),2)*as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1))) +
//         num2_prime*((M-w)*pow(as_scalar(par.row(dref2-1)),2)*as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1))))/dist ;
//
//
//       num.row(w) = dentau*(num1+num2) ;
//       num_prime.row(w) = dentau*((num1_prime+num2_prime)/dist) ;
//
//       num_prime_mix.row(w) = dentau*((num1_prime_mix+num2_prime_mix)/dist) ;
//
//       // num_prime2.row(w) = dentau*num1_prime2 ;
//
//       num_prime2_mix.row(w) = dentau*num1_prime2_mix ;
//     }
//
//     b = sum(num) ;
//     b_prime = sum(num_prime) ;
//     b_prime_mix = sum(num_prime_mix) ;
//
//     g = b_prime*dist ;
//     // g_prime = sum(num_prime2) ;
//     g_prime_mix = sum(num_prime2_mix) ;
//
//     h = b*dist ;
//     // h_prime = b*dist_prime + b_prime*dist ;
//     h_prime_mix = b*dist_prime + b_prime_mix*dist ;
//
//     d2mix_ln = (f*e_prime_mix-e*f_prime_mix)/pow(f,2) - (h*g_prime_mix-g*h_prime_mix)/pow(h,2) ;
//     com_d2mix_ln = com_d2mix_ln + d2mix_ln ;
//
//   }
//
//   return (com_d2mix_ln) ;
//
//   //   return List::create(Named("b")=b,Named("dist_prime")=dist_prime,
//   //                       Named("b_prime_mix")=b_prime_mix,Named("dist")=dist,
//   //                       Named("e")=e,Named("e_prime_mix")=e_prime_mix,
//   //                       Named("f")=f,Named("f_prime_mix")=f_prime_mix,
//   //                       Named("g")=g,Named("g_prime_mix")=g_prime_mix,
//   //                       Named("h")=h,Named("h_prime_mix")=h_prime_mix) ;
//
// }
//
//
//
//
// ///////////////////////////////////////////
//
// #include <RcppArmadillo.h>
//
// // [[Rcpp::depends("RcppArmadillo")]]
//
// using namespace Rcpp;
//
// double t2mix_cpp(
//     arma::colvec par,
//     arma::mat Theta,
//     const int D,
//     const int C,
//     arma::colvec Z,
//     const int tauref,
//     const int tauref2) {
//
//   const int J = Z.n_rows ;
//   const int M = 2*C+1 ;
//   double tau = 0 ;
//   double tau_prime = 0 ;
//   double dentau = 0 ;
//   double dentau_prime = 0 ;
//   double sumtau = 0 ;
//   double dist = 0 ;
//   double sumdist = 0 ;
//   double num1 = 0 ;
//   double num2 = 0 ;
//   double b = 0 ;
//   double b_prime = 0 ;
//   double c_prime = 0 ;
//   double cnum_prime = 0 ;
//   double g = 0 ;
//   // double g_prime = 0 ;
//   double g_prime_mix = 0 ;
//   double h = 0 ;
//   // double h_prime = 0 ;
//   double h_prime_mix = 0 ;
//   //double t1_ln = 0 ;
//   // double com_t1_ln = 0 ;
//   // double t2_ln = 0 ;
//   // double com_t2_ln = 0 ;
//   double t2mix_ln = 0 ;
//   double com_t2mix_ln = 0 ;
//   double x1 = 0 ;
//   double x2 = 0 ;
//   double U_z = 0 ;
//   double U_w = 0 ;
//   double U_w2 = 0 ;
//   arma::colvec num =  arma::colvec(C+1) ;
//   arma::colvec num_prime = arma::colvec(C+1) ;
//   arma::colvec num_prime_mix = arma::colvec(C+1) ;
//   arma::colvec num_prime2 = arma::colvec(C+1) ;
//   arma::colvec num_prime2_mix = arma::colvec(C+1) ;
//
//   for (int j=0; j<=(J-1); j++) {
//
//     tau = 0 ;
//     tau_prime = 0 ;
//     dentau = 0 ;
//     dentau_prime = 0 ;
//     sumtau = 0 ;
//     dist = 0 ;
//     sumdist = 0 ;
//     num1 = 0 ;
//     num2 = 0 ;
//     b = 0 ;
//     b_prime = 0 ;
//     g = 0 ;
//     // g_prime = 0 ;
//     g_prime_mix = 0 ;
//     h = 0 ;
//     // h_prime = 0 ;
//     h_prime_mix = 0 ;
//     c_prime = 0 ;
//     cnum_prime = 0 ;
//     // t1_ln = 0 ;
//     // t2_ln = 0 ;
//     t2mix_ln = 0 ;
//     x1 = 0 ;
//     x2 = 0 ;
//
//     for (int d=0; d<=(D-1);d++) {
//       dist = pow(as_scalar(par.row(d)),2)*pow(as_scalar(Theta(j,d)-par.row(D+d)),2) ;
//       sumdist = dist+sumdist ;
//     }
//     dist = sqrt(sumdist) ;
//
//     for (int w=0; w<=C;w++) {
//       x1 = w*dist ;
//       x2 = (M-w)*dist ;
//
//       if (w>0) {
//         for (int d=0; d<=(D-1);d++) {
//           tau = as_scalar(par.row(d)*par.row(w+2*D-1)) ;
//           sumtau = tau + sumtau ;
//         }
//       }
//
//       if (w==as_scalar(Z.row(j))) {
//         c_prime = 0 ;
//         if (w>0) {
//           for (int d=0; d<=(D-1);d++) {
//             cnum_prime = as_scalar(par.row(d)) ;
//             c_prime = c_prime + cnum_prime ;
//           } }
//
//         U_z = 0 ;
//
//         if (tauref <= as_scalar(Z.row(j))) {U_z = 1 ; }
//
//         c_prime = U_z*c_prime ;
//       }
//
//       dentau = exp(sumtau) ;
//       dentau_prime = 0 ;
//       if (w>0) {
//         for (int d=0; d<=(D-1);d++) {
//           tau_prime = as_scalar(par.row(d)) ;
//           dentau_prime = tau_prime + dentau_prime ;
//         }
//       }
//       num1 = exp(x1) ;
//       num2 = exp(x2) ;
//
//       U_w = 0 ;
//
//       if (tauref <= w) {U_w = 1 ; }
//
//       U_w2 = 0 ;
//
//       if (tauref2 <= w) {U_w2 = 1 ; }
//
//       num.row(w) = dentau*(num1+num2) ;
//       num_prime.row(w) = U_w*dentau_prime*dentau*(num1+num2) ;
//       num_prime_mix.row(w) = U_w2*dentau_prime*dentau*(num1+num2) ;
//       // num_prime2.row(w) = U_w*pow(dentau_prime,2)*dentau*(num1+num2) ;
//       num_prime2_mix.row(w) = U_w*dentau_prime*U_w2*dentau_prime*dentau*(num1+num2) ;
//     }
//
//     b = sum(num) ;
//     b_prime = sum(num_prime) ;
//     // t1_ln = c_prime - (b_prime/b) ;
//     // com_t1_ln = com_t1_ln + t1_ln ;
//
//     g = b_prime ;
//     // g_prime = sum(num_prime2) ;
//     g_prime_mix = sum(num_prime2_mix) ;
//     h = b ;
//     // h_prime = b_prime ;
//     h_prime_mix = sum(num_prime_mix) ;
//
//     //t2_ln = -(h*g_prime - g*h_prime)/pow(h,2) ;
//     //com_t2_ln = com_t2_ln + t2_ln ;
//     t2mix_ln = -(h*g_prime_mix - g*h_prime_mix)/pow(h,2) ;
//     com_t2mix_ln = com_t2mix_ln + t2mix_ln ;
//
//   }
//
//   return (com_t2mix_ln) ;
//
//   //   return List::create(Named("c_prime")=c_prime,Named("a_prime")=a_prime,
//   //                       Named("a")=a,Named("b_prime")=b_prime,
//   //                       Named("b")=b,Named("dist")=dist,
//   //                       Named("x1")=x1,Named("x2")=x2) ;
//
//   //   return List::create(Named("g")=g,Named("g_prime_mix")=g_prime_mix,
//   //                       Named("h")=h,Named("h_prime_mix")=h_prime_mix) ;
//
//
// }
//
//
//
//
// ////////////////////////////////////////////
//
// #include <RcppArmadillo.h>
//
// // [[Rcpp::depends("RcppArmadillo")]]
//
// using namespace Rcpp;
//
// double da2cross_cpp(
//     arma::colvec par,
//     arma::mat Theta,
//     const int D,
//     const int C,
//     arma::colvec Z,
//     const int dref,
//     const int dref2) {
//
//   const int J = Z.n_rows ;
//   const int M = 2*C+1 ;
//   double tau = 0 ;
//   double tau_prime = 0 ;
//   double dentau = 0 ;
//   double dentau_prime = 0 ;
//   double sumtau = 0 ;
//   double dist = 0 ;
//   double dist_prime = 0 ;
//   double dist_prime_mix = 0 ;
//   double dist_prime_cross = 0 ;
//   double sumdist = 0 ;
//   double num1 = 0 ;
//   double num2 = 0 ;
//   double num1_prime = 0 ;
//   double num2_prime = 0 ;
//   double num1_prime_mix = 0 ;
//   double num2_prime_mix = 0 ;
//   double num1_prime_cross = 0 ;
//   double num2_prime_cross = 0 ;
//   double num1_z_prime = 0 ;
//   double num2_z_prime = 0 ;
//   double num1_prime2 = 0 ;
//   // double num1_prime2_mix = 0 ;
//   double num1_prime2_cross = 0 ;
//   double a = 0 ;
//   double a_prime = 0 ;
//   // double a_prime_mix = 0 ;
//   double a_prime_cross = 0 ;
//   double b = 0 ;
//   double b_prime = 0 ;
//   // double b_prime_mix = 0 ;
//   double b_prime_cross = 0 ;
//   double c_prime = 0 ;
//   double cnum_prime = 0 ;
//   double e = 0 ;
//   // double e_prime = 0 ;
//   // double e_prime_mix = 0 ;
//   double e_prime_cross = 0 ;
//   double f = 0 ;
//   // double f_prime = 0 ;
//   // double f_prime_mix = 0 ;
//   double f_prime_cross = 0 ;
//   double g = 0 ;
//   // double g_prime = 0 ;
//   // double g_prime_mix = 0 ;
//   double g_prime_cross = 0 ;
//   double h = 0 ;
//   // double h_prime = 0 ;
//   // double h_prime_mix = 0 ;
//   double h_prime_cross = 0 ;
//   // double a2_ln = 0 ;
//   // double a2mix_ln = 0 ;
//   double da2cross_ln = 0 ;
//   // double com_a2_ln = 0 ;
//   // double com_a2mix_ln = 0 ;
//   double com_da2cross_ln = 0 ;
//   double x1 = 0 ;
//   double x2 = 0 ;
//   arma::colvec num =  arma::colvec(C+1) ;
//   arma::colvec num_prime = arma::colvec(C+1) ;
//   arma::colvec num_prime_mix = arma::colvec(C+1) ;
//   arma::colvec num_prime_cross = arma::colvec(C+1) ;
//   arma::colvec num_prime2 = arma::colvec(C+1) ;
//   arma::colvec num_prime2_mix = arma::colvec(C+1) ;
//   arma::colvec num_prime2_cross = arma::colvec(C+1) ;
//   arma::colvec num_prime_test = arma::colvec(C+1) ;
//
//   for (int j=0; j<=(J-1); j++) {
//
//     tau = 0 ;
//     tau_prime = 0 ;
//     dentau = 0 ;
//     dentau_prime = 0 ;
//     sumtau = 0 ;
//     dist = 0 ;
//     dist_prime = 0 ;
//     dist_prime_mix = 0 ;
//     dist_prime_cross = 0 ;
//     sumdist = 0 ;
//     num1 = 0 ;
//     num2 = 0 ;
//     num1_prime = 0 ;
//     num2_prime = 0 ;
//     num1_prime_mix = 0 ;
//     num2_prime_mix = 0 ;
//     num1_prime_cross = 0 ;
//     num2_prime_cross = 0 ;
//     num1_z_prime = 0 ;
//     num2_z_prime = 0 ;
//     num1_prime2 = 0 ;
//     // num1_prime2_mix = 0 ;
//     num1_prime2_cross = 0 ;
//     a = 0 ;
//     a_prime = 0 ;
//     // a_prime_mix = 0 ;
//     a_prime_cross = 0 ;
//     b = 0 ;
//     b_prime = 0 ;
//     // b_prime_mix = 0 ;
//     b_prime_cross = 0 ;
//     c_prime = 0 ;
//     cnum_prime = 0 ;
//     e = 0 ;
//     // e_prime = 0 ;
//     // e_prime_mix = 0 ;
//     e_prime_cross = 0 ;
//     f = 0 ;
//     // f_prime = 0 ;
//     // f_prime_mix = 0 ;
//     f_prime_cross = 0 ;
//     g = 0 ;
//     // g_prime = 0 ;
//     // g_prime_mix = 0 ;
//     g_prime_cross = 0 ;
//     h = 0 ;
//     // h_prime = 0 ;
//     // h_prime_mix = 0 ;
//     h_prime_cross = 0 ;
//     // a2_ln = 0 ;
//     // a2mix_ln = 0 ;
//     da2cross_ln = 0 ;
//     x1 = 0 ;
//     x2 = 0 ;
//
//     for (int d=0; d<=(D-1);d++) {
//       dist = pow(as_scalar(par.row(d)),2)*pow(as_scalar(Theta(j,d)-par.row(D+d)),2) ;
//       sumdist = dist+sumdist ;
//     }
//     dist = sqrt(sumdist) ;
//     dist_prime = as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist ;
//     dist_prime_mix = dist_prime ;
//     dist_prime_cross = as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist ;
//
//
//     for (int w=0; w<=C;w++) {
//       x1 = w*dist ;
//       x2 = (M-w)*dist ;
//
//       if (w>0) {
//         for (int d=0; d<=(D-1);d++) {
//           tau = as_scalar(par.row(d)*par.row(w+2*D-1)) ;
//           sumtau = tau + sumtau ;
//         }
//       }
//
//       if (w==as_scalar(Z.row(j))) {
//         num1_z_prime = as_scalar(Z.row(j)*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1))*exp(x1)) ;
//         num2_z_prime = as_scalar((M-Z.row(j))*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1))*exp(x2)) ;
//
//         a = exp(x1)+exp(x2) ;
//         a_prime = (num1_z_prime+num2_z_prime)/dist ;
//         //       a_prime_mix = (as_scalar(Z.row(j)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*exp(x1) +
//         //         as_scalar((M-Z.row(j))*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*exp(x2))/dist ;
//
//         a_prime_cross = (as_scalar(Z.row(j)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*exp(x1) +
//           as_scalar((M-Z.row(j))*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*exp(x2))/dist ;
//
//
//         e = a_prime*dist ;
//         //       e_prime = as_scalar(Z.row(j)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1) +
//         //         (M-Z.row(j))*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2))+
//         //         ((pow(as_scalar(Z.row(j)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x1)) +
//         //         (pow(as_scalar((M-Z.row(j))*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x2)))/dist ;
//
//         //       e_prime_mix = (as_scalar(num1_z_prime*(Z.row(j)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))) +
//         //         as_scalar(num2_z_prime*((M-Z.row(j))*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))))/dist ;
//
//         e_prime_cross = (as_scalar(Z.row(j)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*num1_z_prime +
//           as_scalar((M-Z.row(j))*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*num2_z_prime)/dist +
//           as_scalar((2*Z.row(j)*par.row(dref2-1)*(par.row(dref2-1+D)-Theta(j,dref2-1)))*exp(x1)) +
//           as_scalar((2*(M-Z.row(j))*par.row(dref2-1)*(par.row(dref2-1+D)-Theta(j,dref2-1)))*exp(x2)) ;
//
//         f = a*dist ;
//         //       f_prime = a*(as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist) +
//         //                   dist*a_prime ;
//         //       f_prime_mix = a*dist_prime_mix + dist*a_prime_mix ;
//
//         f_prime_cross = a*dist_prime_cross + dist*a_prime_cross ;
//
//         c_prime = 0 ;
//         if (w>0) {
//           for (int t=0; t<=(as_scalar(Z.row(j))-1);t++) {
//             cnum_prime = as_scalar(par.row(t+2*D)) ;
//             c_prime = c_prime + cnum_prime ;
//           }
//         }
//       }
//
//       dentau = exp(sumtau) ;
//       dentau_prime = 0 ;
//       if (w>0) {
//         for (int t=0; t<=(w-1);t++) {
//           tau_prime = as_scalar(par.row(t+2*D)) ;
//           dentau_prime = tau_prime + dentau_prime ;
//         }
//       }
//       num1 = exp(x1) ;
//       num2 = exp(x2) ;
//
//       num1_prime = as_scalar(w*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1))*exp(x1)) ;
//       num2_prime = as_scalar((M-w)*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1))*exp(x2)) ;
//
//       num1_prime_mix = as_scalar(w*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1)) ;
//       num2_prime_mix = as_scalar((M-w)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2)) ;
//
//       num1_prime_cross = as_scalar(w*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1)) ;
//       num2_prime_cross = as_scalar((M-w)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2)) ;
//
//
//       num1_prime2 = (w*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1) +
//         (M-w)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2))+
//         ((pow(as_scalar(w*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x1)) +
//         (pow(as_scalar((M-w)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x2)))/dist ;
//
//       num1_prime2_cross = (as_scalar(w*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*num1_prime +
//         as_scalar((M-w)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*num2_prime)/dist +
//         as_scalar((2*w*par.row(dref2-1)*(par.row(dref2-1+D)-Theta(j,dref2-1)))*exp(x1)) +
//         as_scalar((2*(M-w)*par.row(dref2-1)*(par.row(dref2-1+D)-Theta(j,dref2-1)))*exp(x2)) ;
//
//       num.row(w) = dentau*(num1+num2) ;
//       num_prime.row(w) = dentau*((num1_prime+num2_prime)/dist) ;
//
//       num_prime_mix.row(w) = (dentau*((num1_prime_mix+num2_prime_mix)/dist) + (num1+num2)*dentau*dentau_prime) ;
//
//       num_prime_cross.row(w) = (dentau*((num1_prime_cross+num2_prime_cross)/dist) + (num1+num2)*dentau*dentau_prime) ;
//
//       num_prime2.row(w) = dentau*num1_prime2+dentau*dentau_prime*(num1_prime+num2_prime)+
//         (num1+num2)*(dentau*dentau_prime*dist_prime+dentau*dentau_prime*dentau_prime*dist)+
//         ((num1_prime+num2_prime)/dist)*(dentau*dentau_prime*dist) ;
//
//       num_prime2_mix.row(w) = dentau*(num1_prime*(w*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)) +
//         num2_prime*((M-w)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)))/dist +
//         dentau_prime*dentau*(num1_prime+num2_prime)+
//         (num1+num2)*(dentau*dentau_prime*dist_prime_mix+dentau*dentau_prime*dentau_prime*dist) +
//         (((w*(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)))*num1 +
//         ((M-w)*(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)))*num2)/dist)*(dentau*dentau_prime*dist) ;
//
//       num_prime2_cross.row(w) = dentau*num1_prime2_cross+dentau*dentau_prime*(num1_prime+num2_prime) ;
//
//     }
//
//     b = sum(num) ;
//     b_prime = sum(num_prime) ;
//     // b_prime_mix = sum(num_prime_mix) ;
//     b_prime_cross = sum(num_prime_cross) ;
//
//     g = b_prime*dist ;
//     // g_prime = sum(num_prime2) ;
//     // g_prime_mix = sum(num_prime2_mix) ;
//     g_prime_cross = sum(num_prime2_cross) ;
//
//     h = b*dist ;
//     // h_prime = b*dist_prime + b_prime*dist ;
//     // h_prime_mix = b*dist_prime_mix + b_prime_mix*dist ;
//     h_prime_cross = b*dist_prime_cross + b_prime_cross*dist ;
//
//     // a2_ln = (f*e_prime-e*f_prime)/pow(f,2) - (h*g_prime-g*h_prime)/pow(h,2) ;
//     // a2mix_ln = (f*e_prime_mix-e*f_prime_mix)/pow(f,2) - (h*g_prime_mix-g*h_prime_mix)/pow(h,2) ;
//     da2cross_ln = (f*e_prime_cross-e*f_prime_cross)/pow(f,2) - (h*g_prime_cross-g*h_prime_cross)/pow(h,2) ;
//     // com_a2_ln = com_a2_ln + a2_ln ;
//     // com_a2mix_ln = com_a2mix_ln + a2mix_ln ;
//     com_da2cross_ln = com_da2cross_ln + da2cross_ln ;
//
//   }
//
//   return (com_da2cross_ln) ;
//
//   //   return List::create(Named("e")=e,Named("e_prime_cross")=e_prime_cross,
//   //                       Named("f")=f,Named("f_prime_cross")=f_prime_cross,
//   //                       Named("g")=g,Named("g_prime_cross")=g_prime_cross,
//   //                       Named("h")=h,Named("h_prime_cross")=h_prime_cross) ;
//
//   //  return List::create(Named("a_prime")=a_prime) ;
//
// }
//
//
//
//
// ////////////////////////////////////////
//
// #include <RcppArmadillo.h>
//
// // [[Rcpp::depends("RcppArmadillo")]]
//
// using namespace Rcpp;
//
// double da2crossmix_cpp(
//     arma::colvec par,
//     arma::mat Theta,
//     const int D,
//     const int C,
//     arma::colvec Z,
//     const int dref,
//     const int dref2) {
//
//   const int J = Z.n_rows ;
//   const int M = 2*C+1 ;
//   double tau = 0 ;
//   double tau_prime = 0 ;
//   double dentau = 0 ;
//   double dentau_prime = 0 ;
//   double sumtau = 0 ;
//   double dist = 0 ;
//   double dist_prime = 0 ;
//   // double dist_prime_mix = 0 ;
//   double dist_prime_cross = 0 ;
//   // double dist_prime_cross2 = 0 ;
//   double sumdist = 0 ;
//   double num1 = 0 ;
//   double num2 = 0 ;
//   double num1_prime = 0 ;
//   double num2_prime = 0 ;
//   // double num1_prime_mix = 0 ;
//   // double num2_prime_mix = 0 ;
//   double num1_prime_cross = 0 ;
//   double num2_prime_cross = 0 ;
//   // double num1_prime_cross2 = 0 ;
//   // double num2_prime_cross2 = 0 ;
//   double num1_z_prime = 0 ;
//   double num2_z_prime = 0 ;
//   double num1_prime2 = 0 ;
//   // double num1_prime2_mix = 0 ;
//   double num1_prime2_cross = 0 ;
//   double num1_prime2_cross2 = 0 ;
//   double a = 0 ;
//   double a_prime = 0 ;
//   // double a_prime_mix = 0 ;
//   double a_prime_cross = 0 ;
//   // double a_prime_cross2 = 0 ;
//   double b = 0 ;
//   double b_prime = 0 ;
//   // double b_prime_mix = 0 ;
//   double b_prime_cross = 0 ;
//   // double b_prime_cross2 = 0 ;
//   double c_prime = 0 ;
//   double cnum_prime = 0 ;
//   double e = 0 ;
//   // double e_prime = 0 ;
//   // double e_prime_mix = 0 ;
//   // double e_prime_cross = 0 ;
//   double e_prime_cross2 = 0 ;
//   double f = 0 ;
//   // double f_prime = 0 ;
//   // double f_prime_mix = 0 ;
//   double f_prime_cross = 0 ;
//   // double f_prime_cross2 = 0 ;
//   double g = 0 ;
//   // double g_prime = 0 ;
//   // double g_prime_mix = 0 ;
//   // double g_prime_cross = 0 ;
//   double g_prime_cross2 = 0 ;
//   double h = 0 ;
//   // double h_prime = 0 ;
//   // double h_prime_mix = 0 ;
//   double h_prime_cross = 0 ;
//   // double h_prime_cross2 = 0 ;
//   // double a2_ln = 0 ;
//   // double a2mix_ln = 0 ;
//   double da2crossmix_ln = 0 ;
//   // double com_a2_ln = 0 ;
//   // double com_a2mix_ln = 0 ;
//   double com_da2crossmix_ln = 0 ;
//   double x1 = 0 ;
//   double x2 = 0 ;
//   arma::colvec num =  arma::colvec(C+1) ;
//   arma::colvec num_prime = arma::colvec(C+1) ;
//   // arma::colvec num_prime_mix = arma::colvec(C+1) ;
//   arma::colvec num_prime_cross = arma::colvec(C+1) ;
//   arma::colvec num_prime_cross2 = arma::colvec(C+1) ;
//   arma::colvec num_prime2 = arma::colvec(C+1) ;
//   // arma::colvec num_prime2_mix = arma::colvec(C+1) ;
//   arma::colvec num_prime2_cross = arma::colvec(C+1) ;
//   arma::colvec num_prime2_cross2 = arma::colvec(C+1) ;
//   // arma::colvec num_prime_test = arma::colvec(C+1) ;
//
//   for (int j=0; j<=(J-1); j++) {
//
//     tau = 0 ;
//     tau_prime = 0 ;
//     dentau = 0 ;
//     dentau_prime = 0 ;
//     sumtau = 0 ;
//     dist = 0 ;
//     dist_prime = 0 ;
//     // dist_prime_mix = 0 ;
//     // dist_prime_cross2 = 0 ;
//     sumdist = 0 ;
//     num1 = 0 ;
//     num2 = 0 ;
//     num1_prime = 0 ;
//     num2_prime = 0 ;
//     //     num1_prime_mix = 0 ;
//     //     num2_prime_mix = 0 ;
//     num1_prime_cross = 0 ;
//     num2_prime_cross = 0 ;
//     //     num1_prime_cross2 = 0 ;
//     //     num2_prime_cross2 = 0 ;
//     num1_z_prime = 0 ;
//     num2_z_prime = 0 ;
//     num1_prime2 = 0 ;
//     // num1_prime2_mix = 0 ;
//     num1_prime2_cross = 0 ;
//     num1_prime2_cross2 = 0 ;
//     a = 0 ;
//     a_prime = 0 ;
//     // a_prime_mix = 0 ;
//     a_prime_cross = 0 ;
//     // a_prime_cross2 = 0 ;
//     b = 0 ;
//     b_prime = 0 ;
//     // b_prime_mix = 0 ;
//     b_prime_cross = 0 ;
//     // b_prime_cross2 = 0 ;
//     c_prime = 0 ;
//     cnum_prime = 0 ;
//     e = 0 ;
//     // e_prime = 0 ;
//     // e_prime_mix = 0 ;
//     // e_prime_cross = 0 ;
//     e_prime_cross2 = 0 ;
//     f = 0 ;
//     // f_prime = 0 ;
//     // f_prime_mix = 0 ;
//     f_prime_cross = 0 ;
//     // f_prime_cross2 = 0 ;
//     g = 0 ;
//     // g_prime = 0 ;
//     // g_prime_mix = 0 ;
//     // g_prime_cross = 0 ;
//     g_prime_cross2 = 0 ;
//     h = 0 ;
//     // h_prime = 0 ;
//     // h_prime_mix = 0 ;
//     h_prime_cross = 0 ;
//     // h_prime_cross2 = 0 ;
//     // a2_ln = 0 ;
//     // a2mix_ln = 0 ;
//     da2crossmix_ln = 0 ;
//     x1 = 0 ;
//     x2 = 0 ;
//
//     for (int d=0; d<=(D-1);d++) {
//       dist = pow(as_scalar(par.row(d)),2)*pow(as_scalar(Theta(j,d)-par.row(D+d)),2) ;
//       sumdist = dist+sumdist ;
//     }
//     dist = sqrt(sumdist) ;
//     dist_prime = as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist ;
//     // dist_prime_mix = dist_prime ;
//     dist_prime_cross = as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist ;
//
//     for (int w=0; w<=C;w++) {
//       x1 = w*dist ;
//       x2 = (M-w)*dist ;
//
//       if (w>0) {
//         for (int d=0; d<=(D-1);d++) {
//           tau = as_scalar(par.row(d)*par.row(w+2*D-1)) ;
//           sumtau = tau + sumtau ;
//         }
//       }
//
//       if (w==as_scalar(Z.row(j))) {
//         num1_z_prime = as_scalar(Z.row(j)*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1))*exp(x1)) ;
//         num2_z_prime = as_scalar((M-Z.row(j))*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1))*exp(x2)) ;
//
//         a = exp(x1)+exp(x2) ;
//         a_prime = (num1_z_prime+num2_z_prime)/dist ;
//         //       a_prime_mix = (as_scalar(Z.row(j)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*exp(x1) +
//         //         as_scalar((M-Z.row(j))*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*exp(x2))/dist ;
//
//         a_prime_cross = (as_scalar(Z.row(j)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*exp(x1) +
//           as_scalar((M-Z.row(j))*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*exp(x2))/dist ;
//
//
//         e = a_prime*dist ;
//
//         //       e_prime = as_scalar(Z.row(j)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1) +
//         //         (M-Z.row(j))*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2))+
//         //         ((pow(as_scalar(Z.row(j)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x1)) +
//         //         (pow(as_scalar((M-Z.row(j))*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x2)))/dist ;
//
//         //       e_prime_mix = (as_scalar(num1_z_prime*(Z.row(j)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))) +
//         //         as_scalar(num2_z_prime*((M-Z.row(j))*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))))/dist ;
//
//         //       e_prime_cross = (as_scalar(Z.row(j)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*num1_z_prime +
//         //           as_scalar((M-Z.row(j))*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*num2_z_prime)/dist +
//         //           as_scalar((2*Z.row(j)*par.row(dref2-1)*(par.row(dref2-1+D)-Theta(j,dref2-1)))*exp(x1)) +
//         //           as_scalar((2*(M-Z.row(j))*par.row(dref2-1)*(par.row(dref2-1+D)-Theta(j,dref2-1)))*exp(x2)) ;
//
//         e_prime_cross2 = (as_scalar(Z.row(j)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*num1_z_prime +
//           as_scalar((M-Z.row(j))*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*num2_z_prime)/dist ;
//
//         f = a*dist ;
//         //       f_prime = a*(as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist) +
//         //                   dist*a_prime ;
//         //       f_prime_mix = a*dist_prime_mix + dist*a_prime_mix ;
//
//         f_prime_cross = a*dist_prime_cross + dist*a_prime_cross ;
//         //      f_prime_cross2 = a*dist_prime_cross + dist*a_prime_cross ;
//
//         c_prime = 0 ;
//         if (w>0) {
//           for (int t=0; t<=(as_scalar(Z.row(j))-1);t++) {
//             cnum_prime = as_scalar(par.row(t+2*D)) ;
//             c_prime = c_prime + cnum_prime ;
//           }
//         }
//       }
//
//       dentau = exp(sumtau) ;
//       dentau_prime = 0 ;
//       if (w>0) {
//         for (int t=0; t<=(w-1);t++) {
//           tau_prime = as_scalar(par.row(t+2*D)) ;
//           dentau_prime = tau_prime + dentau_prime ;
//         }
//       }
//       num1 = exp(x1) ;
//       num2 = exp(x2) ;
//
//       num1_prime = as_scalar(w*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1))*exp(x1)) ;
//       num2_prime = as_scalar((M-w)*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1))*exp(x2)) ;
//
//       // num1_prime_mix = as_scalar(w*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1)) ;
//       // num2_prime_mix = as_scalar((M-w)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2)) ;
//
//       num1_prime_cross = as_scalar(w*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1)) ;
//       num2_prime_cross = as_scalar((M-w)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2)) ;
//
//
//       num1_prime2 = (w*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1) +
//         (M-w)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2))+
//         ((pow(as_scalar(w*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x1)) +
//         (pow(as_scalar((M-w)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x2)))/dist ;
//
//       num1_prime2_cross = (as_scalar(w*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*num1_prime +
//         as_scalar((M-w)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*num2_prime)/dist +
//         as_scalar((2*w*par.row(dref2-1)*(par.row(dref2-1+D)-Theta(j,dref2-1)))*exp(x1)) +
//         as_scalar((2*(M-w)*par.row(dref2-1)*(par.row(dref2-1+D)-Theta(j,dref2-1)))*exp(x2)) ;
//
//       num1_prime2_cross2 = (as_scalar(w*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*num1_prime +
//         as_scalar((M-w)*par.row(dref2-1)*pow((Theta(j,dref2-1)-par.row(dref2-1+D)),2))*num2_prime)/dist ;
//
//       num.row(w) = dentau*(num1+num2) ;
//       num_prime.row(w) = dentau*((num1_prime+num2_prime)/dist) ;
//
//       // num_prime_mix.row(w) = (dentau*((num1_prime_mix+num2_prime_mix)/dist) + (num1+num2)*dentau*dentau_prime) ;
//
//       num_prime_cross.row(w) = (dentau*((num1_prime_cross+num2_prime_cross)/dist) + (num1+num2)*dentau*dentau_prime) ;
//
//       num_prime2.row(w) = dentau*num1_prime2+dentau*dentau_prime*(num1_prime+num2_prime)+
//         (num1+num2)*(dentau*dentau_prime*dist_prime+dentau*dentau_prime*dentau_prime*dist)+
//         ((num1_prime+num2_prime)/dist)*(dentau*dentau_prime*dist) ;
//
//       //       num_prime2_mix.row(w) = dentau*(num1_prime*(w*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)) +
//       //           num2_prime*((M-w)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)))/dist +
//       //           dentau_prime*dentau*(num1_prime+num2_prime)+
//       //         (num1+num2)*(dentau*dentau_prime*dist_prime_mix+dentau*dentau_prime*dentau_prime*dist) +
//       //         (((w*(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)))*num1 +
//       //         ((M-w)*(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)))*num2)/dist)*(dentau*dentau_prime*dist) ;
//
//       num_prime2_cross.row(w) = dentau*num1_prime2_cross+dentau*dentau_prime*(num1_prime+num2_prime) ;
//       num_prime2_cross2.row(w) = dentau*num1_prime2_cross2+dentau*dentau_prime*(num1_prime+num2_prime) ;
//
//     }
//
//     b = sum(num) ;
//     b_prime = sum(num_prime) ;
//     // b_prime_mix = sum(num_prime_mix) ;
//     b_prime_cross = sum(num_prime_cross) ;
//
//     g = b_prime*dist ;
//     // g_prime = sum(num_prime2) ;
//     // g_prime_mix = sum(num_prime2_mix) ;
//     // g_prime_cross = sum(num_prime2_cross) ;
//     g_prime_cross2 = sum(num_prime2_cross2) ;
//
//     h = b*dist ;
//     // h_prime = b*dist_prime + b_prime*dist ;
//     // h_prime_mix = b*dist_prime_mix + b_prime_mix*dist ;
//     h_prime_cross = b*dist_prime_cross + b_prime_cross*dist ;
//
//     // a2_ln = (f*e_prime-e*f_prime)/pow(f,2) - (h*g_prime-g*h_prime)/pow(h,2) ;
//     // a2mix_ln = (f*e_prime_mix-e*f_prime_mix)/pow(f,2) - (h*g_prime_mix-g*h_prime_mix)/pow(h,2) ;
//     da2crossmix_ln = (f*e_prime_cross2-e*f_prime_cross)/pow(f,2) - (h*g_prime_cross2-g*h_prime_cross)/pow(h,2) ;
//     // com_a2_ln = com_a2_ln + a2_ln ;
//     // com_a2mix_ln = com_a2mix_ln + a2mix_ln ;
//     com_da2crossmix_ln = com_da2crossmix_ln + da2crossmix_ln ;
//
//   }
//
//   return (com_da2crossmix_ln) ;
//
//   //   return List::create(Named("e")=e,Named("e_prime_cross")=e_prime_cross,
//   //                       Named("f")=f,Named("f_prime_cross")=f_prime_cross,
//   //                       Named("g")=g,Named("g_prime_cross")=g_prime_cross,
//   //                       Named("h")=h,Named("h_prime_cross")=h_prime_cross) ;
//
//   //  return List::create(Named("a_prime")=a_prime) ;
//
// }
//
//
//
//
// ///////////////////////////////////////////
//
// #include <RcppArmadillo.h>
//
// // [[Rcpp::depends("RcppArmadillo")]]
//
// using namespace Rcpp;
//
// double dt2cross_cpp(
//     arma::colvec par,
//     arma::mat Theta,
//     const int D,
//     const int C,
//     arma::colvec Z,
//     const int dref,
//     const int dref2,
//     const int tauref) {
//
//   const int J = Z.n_rows ;
//   const int M = 2*C+1 ;
//   double tau = 0 ;
//   double tau_prime = 0 ;
//   double dentau = 0 ;
//   double dentau_prime = 0 ;
//   double sumtau = 0 ;
//   double dist = 0 ;
//   // double dist_prime = 0 ;
//   double sumdist = 0 ;
//   double num1 = 0 ;
//   double num2 = 0 ;
//   double num1_prime = 0 ;
//   double num2_prime = 0 ;
//   //   double num1_z_prime = 0 ;
//   //   double num2_z_prime = 0 ;
//   double num1_prime2 = 0 ;
//   // double a = 0 ;
//   // double a_prime = 0 ;
//   double b = 0 ;
//   double b_prime = 0 ;
//   // double c_prime = 0 ;
//   // double cnum_prime = 0 ;
//   // double e = 0 ;
//   // double e_prime = 0 ;
//   // double f = 0 ;
//   // double f_prime = 0 ;
//   double g = 0 ;
//   // double g_prime = 0 ;
//   double g_prime_cross2 = 0 ;
//   double h = 0 ;
//   // double h_prime = 0 ;
//   double h_prime_cross2 = 0 ;
//   double dt2cross_ln = 0 ;
//   double com_dt2cross_ln = 0 ;
//   double x1 = 0 ;
//   double x2 = 0 ;
//   double U_w = 0 ;
//   arma::colvec num =  arma::colvec(C+1) ;
//   arma::colvec num_prime = arma::colvec(C+1) ;
//   arma::colvec num_prime_cross2 = arma::colvec(C+1) ;
//   arma::colvec num_prime2 = arma::colvec(C+1) ;
//   arma::colvec num_prime2_cross2 = arma::colvec(C+1) ;
//   // arma::colvec num_prime_test = arma::colvec(C+1) ;
//
//   for (int j=0; j<=(J-1); j++) {
//
//     tau = 0 ;
//     tau_prime = 0 ;
//     dentau = 0 ;
//     dentau_prime = 0 ;
//     sumtau = 0 ;
//     dist = 0 ;
//     //   dist_prime = 0 ;
//     sumdist = 0 ;
//     num1 = 0 ;
//     num2 = 0 ;
//     num1_prime = 0 ;
//     num2_prime = 0 ;
//     //     num1_z_prime = 0 ;
//     //     num2_z_prime = 0 ;
//     num1_prime2 = 0 ;
//     //     a = 0 ;
//     //     a_prime = 0 ;
//     b = 0 ;
//     b_prime = 0 ;
//     // c_prime = 0 ;
//     // cnum_prime = 0 ;
//     //     e = 0 ;
//     //     e_prime = 0 ;
//     //     f = 0 ;
//     //     f_prime = 0 ;
//     g = 0 ;
//     // g_prime = 0 ;
//     g_prime_cross2 = 0 ;
//     h = 0 ;
//     // h_prime = 0 ;
//     h_prime_cross2 = 0 ;
//     dt2cross_ln = 0 ;
//     x1 = 0 ;
//     x2 = 0 ;
//     U_w = 0 ;
//
//     for (int d=0; d<=(D-1);d++) {
//       dist = pow(as_scalar(par.row(d)),2)*pow(as_scalar(Theta(j,d)-par.row(D+d)),2) ;
//       sumdist = dist+sumdist ;
//     }
//     dist = sqrt(sumdist) ;
//     //dist_prime = as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist ;
//
//     for (int w=0; w<=C;w++) {
//       x1 = w*dist ;
//       x2 = (M-w)*dist ;
//
//       if (w>0) {
//         for (int d=0; d<=(D-1);d++) {
//           tau = as_scalar(par.row(d)*par.row(w+2*D-1)) ;
//           sumtau = tau + sumtau ;
//         }
//       }
//
//       if (w==as_scalar(Z.row(j))) {
//         //       num1_z_prime = as_scalar((Z.row(j)*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1)))*exp(x1)) ;
//         //       num2_z_prime = as_scalar(((M-Z.row(j))*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1)))*exp(x2)) ;
//
//         //    a = exp(x1)+exp(x2) ;
//         //    a_prime = (num1_z_prime+num2_z_prime)/dist ;
//
//         //  e = a_prime*dist ;
//         //       e_prime = as_scalar(Z.row(j)*pow(par.row(dref2-1),2)*exp(x1) +
//         //         (M-Z.row(j))*pow(par.row(dref2-1),2)*exp(x2))+
//         //         as_scalar((((Z.row(j)*pow(as_scalar(par.row(dref2-1)),2)*pow(as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)),2))*exp(x1)) +
//         //         ((M-Z.row(j))*pow(as_scalar(par.row(dref2-1)),2)*pow(as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)),2))*exp(x2)))/dist ;
//
//         //       f = a*dist ;
//         //       f_prime = a*((pow(as_scalar(par.row(dref2-1)),2)*pow(as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)),2))/dist) +
//         //                   dist*a_prime ;
//         //   dist_prime = ((pow(as_scalar(par.row(dref2-1)),2)*pow(as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)),2))/dist) ;
//
//         //       c_prime = 0 ;
//         //       if (w>0) {
//         //         for (int t=0; t<=(as_scalar(Z.row(j))-1);t++) {
//         //           cnum_prime = as_scalar(par.row(t+2*D)) ;
//         //           c_prime = c_prime + cnum_prime ;
//         //         }
//         //       }
//       }
//
//       dentau = exp(sumtau) ;
//       dentau_prime = 0 ;
//       if (w>0) {
//         for (int d=0; d<=(D-1);d++) {
//           tau_prime = as_scalar(par.row(d)) ;
//           dentau_prime = tau_prime + dentau_prime ;
//         }
//       }
//       num1 = exp(x1) ;
//       num2 = exp(x2) ;
//
//       U_w = 0 ;
//
//       if (tauref <= w) {U_w = 1 ; }
//
//       num1_prime = as_scalar((w*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1)))*exp(x1)) ;
//       num2_prime = as_scalar(((M-w)*pow(par.row(dref-1),2)*(par.row(dref-1+D)-Theta(j,dref-1)))*exp(x2)) ;
//
//       num1_prime2 = w*pow(as_scalar(par.row(dref2-1)),2)*exp(x1) +
//         (M-w)*pow(as_scalar(par.row(dref2-1)),2)*exp(x2) +
//         ((w*pow(as_scalar(par.row(dref2-1)),2)*pow(as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)),2)*exp(x1)) +
//         ((M-w)*pow(as_scalar(par.row(dref2-1)),2)*pow(as_scalar(par.row(dref2-1+D)-Theta(j,dref2-1)),2)*exp(x2)))/dist ;
//
//       num.row(w) = dentau*(num1+num2) ;
//       num_prime.row(w) = dentau*((num1_prime+num2_prime)/dist) ;
//       num_prime_cross2.row(w) = U_w*dentau_prime*dentau*(num1+num2)*dist ;
//       num_prime2.row(w) = dentau*num1_prime2 ;
//       num_prime2_cross2.row(w) = U_w*dentau_prime*dentau*(num1_prime+num2_prime) ;
//     }
//
//     b = sum(num) ;
//     b_prime = sum(num_prime) ;
//
//     g = b_prime*dist ;
//     // g_prime = sum(num_prime2) ;
//     g_prime_cross2 = sum(num_prime2_cross2) ;
//
//     h = b*dist ;
//     // h_prime = b*dist_prime + b_prime*dist ;
//     h_prime_cross2 = sum(num_prime_cross2) ;
//
//     dt2cross_ln = -(h*g_prime_cross2-g*h_prime_cross2)/pow(h,2) ;
//     com_dt2cross_ln = com_dt2cross_ln + dt2cross_ln ;
//
//   }
//
//   return (com_dt2cross_ln) ;
//
//   //   return List::create(Named("g")=g,Named("g_prime_cross2")=g_prime_cross2,
//   //                       Named("h")=h,Named("h_prime_cross2")=h_prime_cross2) ;
//
// }
//
//
//
//
// ///////////////////////////////////////////
//
// #include <RcppArmadillo.h>
//
// // [[Rcpp::depends("RcppArmadillo")]]
//
// using namespace Rcpp;
//
// double at2cross_cpp(
//     arma::colvec par,
//     arma::mat Theta,
//     const int D,
//     const int C,
//     arma::colvec Z,
//     const int dref,
//     const int dref2,
//     const int tauref) {
//
//   const int J = Z.n_rows ;
//   const int M = 2*C+1 ;
//   double tau = 0 ;
//   double tau_prime = 0 ;
//   double dentau = 0 ;
//   double dentau_prime = 0 ;
//   double dentau_prime_a = 0 ;
//   double dentau_prime_t = 0 ;
//   double sumtau = 0 ;
//   double dist = 0 ;
//   double dist_prime = 0 ;
//   double sumdist = 0 ;
//   double num1 = 0 ;
//   double num2 = 0 ;
//   double num1_prime = 0 ;
//   double num2_prime = 0 ;
//   //  double num1_z_prime = 0 ;
//   //  double num2_z_prime = 0 ;
//   double num1_prime2 = 0 ;
//   //   double a = 0 ;
//   //   double a_prime = 0 ;
//   double b = 0 ;
//   double b_prime = 0 ;
//   double c_prime = 0 ;
//   double cnum_prime = 0 ;
//   // double e = 0 ;
//   // double e_prime = 0 ;
//   // double f = 0 ;
//   // double f_prime = 0 ;
//   double g = 0 ;
//   // double g_prime = 0 ;
//   double g_prime_cross = 0 ;
//   double h = 0 ;
//   // double h_prime = 0 ;
//   double h_prime_cross = 0 ;
//   double at2cross_ln = 0 ;
//   double com_at2cross_ln = 0 ;
//   double x1 = 0 ;
//   double x2 = 0 ;
//   double l_prime = 0 ;
//   double lw_prime = 0 ;
//   double U_w = 0 ;
//   arma::colvec num =  arma::colvec(C+1) ;
//   arma::colvec num_prime = arma::colvec(C+1) ;
//   arma::colvec num_prime_cross = arma::colvec(C+1) ;
//   arma::colvec num_prime2 = arma::colvec(C+1) ;
//   arma::colvec num_prime2_cross = arma::colvec(C+1) ;
//   // arma::colvec num_prime_test = arma::colvec(C+1) ;
//
//   for (int j=0; j<=(J-1); j++) {
//
//     tau = 0 ;
//     tau_prime = 0 ;
//     dentau = 0 ;
//     dentau_prime = 0 ;
//     dentau_prime_a = 0 ;
//     dentau_prime_t = 0 ;
//     sumtau = 0 ;
//     dist = 0 ;
//     dist_prime = 0 ;
//     sumdist = 0 ;
//     num1 = 0 ;
//     num2 = 0 ;
//     num1_prime = 0 ;
//     num2_prime = 0 ;
//     //    num1_z_prime = 0 ;
//     //    num2_z_prime = 0 ;
//     num1_prime2 = 0 ;
//     //    a = 0 ;
//     //    a_prime = 0 ;
//     b = 0 ;
//     b_prime = 0 ;
//     c_prime = 0 ;
//     cnum_prime = 0 ;
//     //     e = 0 ;
//     //     e_prime = 0 ;
//     //     f = 0 ;
//     //     f_prime = 0 ;
//     g = 0 ;
//     //    g_prime = 0 ;
//     g_prime_cross = 0 ;
//     h = 0 ;
//     //    h_prime = 0 ;
//     h_prime_cross = 0 ;
//     at2cross_ln = 0 ;
//     x1 = 0 ;
//     x2 = 0 ;
//     l_prime = 0 ;
//     lw_prime = 0 ;
//     U_w = 0 ;
//
//     for (int d=0; d<=(D-1);d++) {
//       dist = pow(as_scalar(par.row(d)),2)*pow(as_scalar(Theta(j,d)-par.row(D+d)),2) ;
//       sumdist = dist+sumdist ;
//     }
//     dist = sqrt(sumdist) ;
//     dist_prime = as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist ;
//
//     for (int w=0; w<=C;w++) {
//       x1 = w*dist ;
//       x2 = (M-w)*dist ;
//
//       if (w>0) {
//         for (int d=0; d<=(D-1);d++) {
//           tau = as_scalar(par.row(d)*par.row(w+2*D-1)) ;
//           sumtau = tau + sumtau ;
//         }
//       }
//
//       if (w==as_scalar(Z.row(j))) {
//         //      num1_z_prime = as_scalar(Z.row(j)*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x1)) ;
//         //      num2_z_prime = as_scalar((M-Z.row(j))*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x2)) ;
//
//         //      a = exp(x1)+exp(x2) ;
//         //      a_prime = (num1_z_prime+num2_z_prime)/dist ;
//
//         //       e = a_prime*dist ;
//         //       e_prime = as_scalar(Z.row(j)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1) +
//         //         (M-Z.row(j))*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2))+
//         //         ((pow(as_scalar(Z.row(j)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x1)) +
//         //         (pow(as_scalar((M-Z.row(j))*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x2)))/dist ;
//
//         //       f = a*dist ;
//         //       f_prime = a*(as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist) +
//         //                   dist*a_prime ;
//         dist_prime = as_scalar(par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2))/dist ;
//
//         c_prime = 0 ;
//         if (w>0) {
//           for (int t=0; t<=(as_scalar(Z.row(j))-1);t++) {
//             cnum_prime = as_scalar(par.row(t+2*D)) ;
//             c_prime = c_prime + cnum_prime ;
//           }
//         }
//         l_prime = 0 ;
//         if (tauref <= as_scalar(Z.row(j))) {l_prime = 1 ; }
//       }
//
//       lw_prime = 0 ;
//       if (tauref <= w) {lw_prime = 1 ; }
//
//       dentau = exp(sumtau) ;
//
//       dentau_prime_a = 0 ;
//       if (w>0) {
//         for (int d=0; d<=(D-1);d++) {
//           tau_prime = as_scalar(par.row(d)) ;
//           dentau_prime_a = tau_prime + dentau_prime_a ;
//         }
//       }
//
//       U_w = 0 ;
//       if (tauref <= w) {U_w = 1 ; }
//
//       dentau_prime_t = 0 ;
//       if (w>0) {
//         for (int t=0; t<=(w-1);t++) {
//           tau_prime = as_scalar(par.row(t+2*D)) ;
//           dentau_prime_t = tau_prime + dentau_prime_t ;
//         }
//       }
//       num1 = exp(x1) ;
//       num2 = exp(x2) ;
//
//       num1_prime = as_scalar(w*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x1)) ;
//       num2_prime = as_scalar((M-w)*par.row(dref-1)*pow((Theta(j,dref-1)-par.row(dref-1+D)),2)*exp(x2)) ;
//
//       num1_prime2 = (w*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x1) +
//         (M-w)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)*exp(x2))+
//         ((pow(as_scalar(w*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x1)) +
//         (pow(as_scalar((M-w)*par.row(dref2-1)*pow(as_scalar(Theta(j,dref2-1)-par.row(dref2-1+D)),2)),2)*exp(x2)))/dist ;
//
//       num.row(w) = dentau*(num1+num2) ;
//       num_prime.row(w) = (dentau*((num1_prime+num2_prime)/dist) + (num1+num2)*dentau*dentau_prime_t) ;
//       num_prime_cross.row(w) = U_w*dentau_prime_a*dentau*(num1+num2)*dist ;
//
//       num_prime2.row(w) = dentau*num1_prime2+dentau*dentau_prime*(num1_prime+num2_prime)+
//         (num1+num2)*(dentau*dentau_prime*dist_prime+dentau*dentau_prime*dentau_prime*dist)+
//         ((num1_prime+num2_prime)/dist)*(dentau*dentau_prime*dist) ;
//       num_prime2_cross.row(w) = U_w*dentau_prime_a*dentau*(num1_prime+num2_prime) +
//         (num1+num2)*dist*(lw_prime*dentau+U_w*dentau_prime_t*dentau_prime_a*dentau) ;
//     }
//
//     b = sum(num) ;
//     b_prime = sum(num_prime) ;
//
//     g = b_prime*dist ;
//     // g_prime = sum(num_prime2) ;
//     g_prime_cross = sum(num_prime2_cross) ;
//
//     h = b*dist ;
//     // h_prime = b*dist_prime + b_prime*dist ;
//     h_prime_cross = sum(num_prime_cross) ;
//
//     at2cross_ln = l_prime - (h*g_prime_cross-g*h_prime_cross)/pow(h,2) ;
//     com_at2cross_ln = com_at2cross_ln + at2cross_ln ;
//
//   }
//
//   return (com_at2cross_ln) ;
//
//   //   return List::create(Named("l_prime")=l_prime,
//   //                       Named("g")=g,Named("g_prime_cross")=g_prime_cross,
//   //                       Named("h")=h,Named("h_prime_cross")=h_prime_cross) ;
//
// }
//
//
//
//
// /////////////////////////////////////////////
//
// #include <RcppArmadillo.h>
//
// // [[Rcpp::depends("RcppArmadillo")]]
//
// using namespace Rcpp;
//
// arma::mat ggum_hess (
//     arma::colvec par,
//     arma::mat Theta,
//     int D,
//     int C,
//     arma::colvec Z) {
//
//   int dref ;
//   int dref2 ;
//   int tauref ;
//   int tauref2 ;
//   int ind ;
//   int ind2 ;
//   int dim = (2*D)+C ;
//   arma::mat hess =  arma::mat(dim,dim) ;
//   //int d=1 ;
//   //hess(1,1) = a2_cpp(par,Theta,D,C,Z,dref=d,dref2=d) ;
//
//   for (int d=0; d<=(D-1);d++) {
//     ind = D+d ;
//     hess(d,d) = a2_cpp(par,Theta,D,C,Z,dref=(d+1),dref2=(d+1)) ;
//     hess(d,ind) = da2cross_cpp(par,Theta,D,C,Z,dref=(d+1),dref2=(d+1)) ;
//     hess(ind,ind) = d2_cpp(par,Theta,D,C,Z,dref=(d+1),dref2=(d+1)) ;
//   }
//
//   for (int d1=0; d1<=(D-1);d1++) {
//     for (int d2=0; d2<=(D-1);d2++) {
//       if (d1 < d2) {
//         ind = D+d1 ;
//         ind2 = D+d2 ;
//         hess(d1,d2) = a2mix_cpp(par,Theta,D,C,Z,dref=(d1+1),dref2=(d2+1)) ;
//         hess(ind,ind2) = d2mix_cpp(par,Theta,D,C,Z,dref=(d1+1),dref2=(d2+1)) ;
//       }
//     }
//   }
//
//   for (int a1=0; a1<=(D-1);a1++) {
//     for (int d1=0; d1<=(D-1);d1++) {
//       if (a1 != d1) {
//         ind = D+d1 ;
//         hess(a1,ind) = da2crossmix_cpp(par,Theta,D,C,Z,dref=(d1+1),dref2=(a1+1)) ;
//       }
//     }
//   }
//
//   for (int t=0; t<=(C-1);t++) {
//     ind = 2*D+t ;
//     hess(ind,ind) = t2_cpp(par,Theta,D,C,Z,tauref=(t+1)) ;
//   }
//
//   for (int t1=0; t1<=(C-1);t1++) {
//     for (int t2=0; t2<=(C-1);t2++) {
//       if (t1 < t2) {
//         ind = 2*D+t1 ;
//         ind2 = 2*D+t2 ;
//         hess(ind,ind2) = t2mix_cpp(par,Theta,D,C,Z,tauref=(t1+1),tauref2=(t2+1)) ;
//       }
//     }
//   }
//
//   for (int d=0; d<=(D-1);d++) {
//     for (int t=0; t<=(C-1);t++) {
//       ind = 2*D+t ;
//       ind2 = D+d ;
//       hess(d,ind) = at2cross_cpp(par,Theta,D,C,Z,dref=(d+1),dref2=(d+1),tauref=(t+1)) ;
//       hess(ind2,ind) = dt2cross_cpp(par,Theta,D,C,Z,dref=(d+1),dref2=(d+1),tauref=(t+1)) ;
//     }
//   }
//
//   hess = symmatu(hess) ;
//
//   return(hess) ;
//
// }
//
// /////////////////////////////////////////////////////////////////////////////
//
// #include <RcppArmadillo.h>
//
// // [[Rcpp::depends("RcppArmadillo")]]
//
// using namespace Rcpp;
//
// // [[Rcpp::export]]
//
// List ggum_derivs (
//     arma::colvec par,
//     arma::mat Theta,
//     int D,
//     int C,
//     arma::colvec Z) {
//
//   int dim = (2*D)+C ;
//   NumericVector grad(dim) ;
//   arma::mat hess =  arma::mat(dim,dim) ;
//
//   grad = ggum_grad(par=par,Theta=Theta,D=D,C=C,Z=Z) ;
//   hess = ggum_hess(par=par,Theta=Theta,D=D,C=C,Z=Z) ;
//
//   return List::create(Named("grad")=grad,Named("hess")=hess) ;
//
// }
