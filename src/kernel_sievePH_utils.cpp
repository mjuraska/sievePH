#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
double Epankercplusplus(double tk, double tvalue, double hband, double delt) {
  double result = 0.0;

  if (std::abs(tk - tvalue) < hband) {
    result = 3.0 / (4.0 * hband) * delt * (1 - std::pow((tk - tvalue) / hband, 2));
  }

  return result;
}




// [[Rcpp::export]]
Rcpp::List estpcomcplusplus(double tau,int KK, arma::ivec N, int NP,
                   arma::mat X, arma::cube ZT, arma::mat DELTA,
                   arma::mat WGHT, arma::vec BETA0) {
  const int mxn = max(N);

  arma::mat BZt(mxn, 1, arma::fill::zeros);
  arma::mat S0(KK, mxn, arma::fill::zeros);
  arma::mat S1(NP, mxn, arma::fill::zeros);
  arma::mat S2(NP, NP, arma::fill::zeros);
  arma::vec U(NP, arma::fill::zeros);
  arma::mat F2(NP, NP, arma::fill::zeros);
  arma::vec BETA(NP);
  arma::mat F(NP, NP, arma::fill::zeros);
  arma::mat var(NP, NP, arma::fill::zeros);

  int KL = 6;
  for(int j = 0; j < NP; j++){
    BETA(j) = BETA0(j);
  }
  for (int Kiter = 1; Kiter <= KL; Kiter++) {
    F.zeros();
    F2.zeros();
    U.zeros();

    for (int ks = 0; ks < KK; ks++) {
      for (int i = 0; i < N(ks); i++) {
        BZt(i, 0) = 0.0;
        for (int j = 0; j < NP; j++) {
          BZt(i, 0) += BETA(j) * ZT(ks, j, i);
        }
      }

      for (int i = 0; i < N(ks); i++) {
        if (X(ks, i) <= tau) {
          S0(ks, i) = 0.0;
          S1.col(i).zeros();
          S2.zeros();

          for (int l = 0; l < N(ks); l++) {
            if (X(ks, l) >= X(ks, i)) {
              S0(ks, i) += exp(BZt(l, 0)) * WGHT(ks, l);

              for (int j = 0; j < NP; j++) {
                S1(j, i) += exp(BZt(l, 0)) * ZT(ks, j, l) * WGHT(ks, l);

                for (int k = 0; k < NP; k++) {
                  S2(j, k) += exp(BZt(l, 0)) * ZT(ks, j, l) * ZT(ks, k, l) * WGHT(ks, l);
                }
              }
            }
          }

          if (S0(ks, i) > 0.00000001) {
            for (int j = 0; j < NP; j++) {
              U(j) += DELTA(ks, i) * (ZT(ks, j, i) - S1(j, i) / S0(ks, i)) * WGHT(ks, i);

              for (int k = 0; k < NP; k++) {
                F(j, k) += DELTA(ks, i) * (S2(j, k) / S0(ks, i) - S1(j, i) * S1(k, i) / pow(S0(ks, i), 2)) * WGHT(ks, i);
                F2(j, k) += pow(DELTA(ks, i), 2) * (ZT(ks, j, i) - S1(j, i) / S0(ks, i)) * (ZT(ks, k, i) - S1(k, i) / S0(ks, i)) * pow(WGHT(ks, i), 2);
              }
            }
          }
        }
      }
    }

    BETA += arma::solve(F, U);
  }

  var = arma::inv(F) * F2 * arma::inv(F);

  return Rcpp::List::create(
    Rcpp::Named("BETA") = BETA,
    Rcpp::Named("F") = F,
    Rcpp::Named("var") = var
  );
}



// [[Rcpp::export]]
Rcpp::List estpipwcplusplus(double tau, double tstep, int ntgrid, double TBAND, int KK, arma::ivec N, int NP,
                   arma::mat X, arma::cube ZT, arma::mat CENSOR, arma::mat DELTA,
                   arma::mat WGHT, arma::vec BETA0, int maxit, int estBaseLamInd) {
  const int mxn = max(N);

  arma::mat BZt(mxn, 1, arma::fill::zeros);
  arma::mat S0(KK, mxn, arma::fill::zeros);
  arma::mat S1(NP, mxn, arma::fill::zeros);
  arma::mat S2(NP, NP, arma::fill::zeros);
  arma::vec U(NP, arma::fill::zeros);
  arma::mat F(NP, NP, arma::fill::zeros);
  arma::mat invF(NP, NP, arma::fill::zeros);
  arma::mat F2(NP, NP, arma::fill::zeros);
  arma::mat var(NP, NP, arma::fill::zeros);
  arma::mat LAMBDA0(KK, mxn, arma::fill::zeros);
  arma::vec BETA (NP);
  arma::mat LAMBDAk0(KK, ntgrid, arma::fill::zeros);
  arma::vec change(NP, arma::fill::ones);
  int Kiter = 1;
  BETA = BETA0;
  
  while (arma::sum(arma::abs(change)) > 0.00001 && Kiter <= maxit) {
  //for (int Kiter = 1; Kiter <= maxit; Kiter++) {
    F.zeros();
    F2.zeros();
    U.zeros();
    S0.zeros();
    for (int ks = 0; ks < KK; ks++) {
      for (int i = 0; i < N(ks); i++) {
        BZt(i, 0) = 0.0;
        for (int j = 0; j < NP; j++) {
          BZt(i, 0) += BETA(j) * ZT(ks, j, i);
        }
      }
      
      for (int i = 0; i < N(ks); i++) {
        if (X(ks, i) <= tau) {
          S0(ks, i) = 0.0;
          S1.col(i).zeros();
          S2.zeros();
          
          for (int l = 0; l < N(ks); l++) {
            if (X(ks, l) >= X(ks, i)) {
              S0(ks, i) += exp(BZt(l, 0)) * WGHT(ks, l);
              
              for (int j = 0; j < NP; j++) {
                S1(j, i) += exp(BZt(l, 0)) * ZT(ks, j, l) * WGHT(ks, l);
                
                for (int k = 0; k < NP; k++) {
                  S2(j, k) += exp(BZt(l, 0)) * ZT(ks, j, l) * ZT(ks, k, l) * WGHT(ks, l);
                }
              }
            }
          }
          
          if (S0(ks, i) > 0.00000001) {
            for (int j = 0; j < NP; j++) {
              U(j) += DELTA(ks, i) * (ZT(ks, j, i) - S1(j, i) / S0(ks, i)) * WGHT(ks, i);
              
              for (int k = 0; k < NP; k++) {
                F(j, k) += DELTA(ks, i) * (S2(j, k) / S0(ks, i) - S1(j, i) * S1(k, i) / pow(S0(ks, i), 2)) * WGHT(ks, i);
                F2(j, k) += pow(DELTA(ks, i), 2) * (ZT(ks, j, i) - S1(j, i) / S0(ks, i)) * (ZT(ks, k, i) - S1(k, i) / S0(ks, i)) * pow(WGHT(ks, i), 2);
              }
            }
          }
        }
      }}
      Kiter = Kiter +1;
      change = arma::solve(F, U);
      BETA += change;
      //Rcpp::Rcout << "F" << F << std::endl;
      //Rcpp::Rcout << "U" << U << std::endl;
      
    };
  

  var = arma::inv(F) * F2 * arma::inv(F);
  invF  = arma::inv(F);
  
  for (int ks = 0; ks < KK; ks++) {
    for (int i = 0; i < N(ks); i++) {
      double TEMPB = 0.0;
      for (int ii = 0; ii < N(ks); ii++) {
        if (X(ks, ii) <= tau && S0(ks, ii) > 0.00000001) {
          TEMPB += Epankercplusplus(X(ks, ii), X(ks, i), TBAND, CENSOR(ks, ii)) * DELTA(ks, ii) * WGHT(ks, ii) / S0(ks, ii);
        }
      }
      LAMBDA0(ks, i) = TEMPB;
    }
  }

  if(estBaseLamInd == 1){
    for (int ks = 0; ks < KK; ks++) {
      for (int Itgrid = 0; Itgrid < ntgrid; Itgrid++) {
        double tvalue = tstep * Itgrid;
        double TEMPB = 0.0;
        for (int ii = 0; ii < N(ks); ii++) {
          if (X(ks, ii) <= tau && S0(ks,ii) > 0.00000001) {
              TEMPB += Epankercplusplus(X(ks, ii), tvalue, TBAND, CENSOR(ks, ii)) * (DELTA(ks, ii) * WGHT(ks, ii)/ S0(ks,ii));
          }
        }
        LAMBDAk0(ks, Itgrid) = TEMPB;
      }
    } 
  }
  
  
  return Rcpp::List::create(
    Rcpp::Named("BETA") = BETA,
    Rcpp::Named("F") = invF,
    Rcpp::Named("var") = var,
    Rcpp::Named("LAMBDA0") = LAMBDA0,
    Rcpp::Named("LAMBDAk0") = LAMBDAk0
  );
}



// [[Rcpp::export]]
Rcpp::List estpaugcplusplus(double tau, double tstep, int ntgrid, double TBAND, int KK, arma::ivec N, int NP,
             arma::mat X, arma::cube ZT, arma::mat CENSOR, arma::mat DELTA,
             arma::mat WGHT, arma::mat DRHOipw, arma::vec BETA0, int maxit, int estBaseLamInd) {
  const int mxn = max(N);
  arma::mat BZt(mxn, 1, arma::fill::zeros);
  arma::mat S0(KK, mxn, arma::fill::zeros);
  arma::mat S1(NP, mxn, arma::fill::zeros);
  arma::mat S2(NP, NP, arma::fill::zeros);
  arma::vec U(NP, arma::fill::zeros);
  arma::mat F2(NP, NP, arma::fill::zeros);
  arma::mat FF2(NP, NP, arma::fill::zeros);
  arma::mat LAMBDA0(KK, mxn, arma::fill::zeros);
  arma::vec BETA(NP);
  arma::mat F(NP, NP, arma::fill::zeros);
  arma::mat invF(NP, NP, arma::fill::zeros);
  arma::mat var(NP, NP, arma::fill::zeros);
  arma::mat LAMBDAk0(KK, ntgrid, arma::fill::zeros);
  arma::vec change(NP, arma::fill::ones);
  int Kiter = 1; 
  
  for(int j = 0; j < NP; j++){
    BETA(j) = BETA0(j);
  }
  while (arma::sum(arma::abs(change)) > 0.00001 && Kiter <= maxit) {
  //for (int Kiter = 1; Kiter <= maxit; Kiter++) {
    F.zeros();
    F2.zeros();
    S0.zeros();
    U.zeros();
    for (int ks = 0; ks < KK; ks++) {
      for (int i = 0; i < N(ks); i++) {
        BZt.row(i).zeros();
        for (int j = 0; j < NP; j++) {
          BZt(i, 0) += BETA(j) * ZT(ks, j, i);
        }
      }

      for (int i = 0; i < N(ks); i++) {
        if (X(ks, i) <= tau) {
          S0(ks, i) = 0.0;
          S1.col(i).zeros();
          S2.zeros();

          for (int L = 0; L < N(ks); L++) {
            if (X(ks, L) >= X(ks, i)) {
              S0(ks, i) += exp(BZt(L, 0));
              for (int j = 0; j < NP; j++) {
                S1(j, i) += exp(BZt(L, 0)) * ZT(ks, j, L);
                for (int k = 0; k < NP; k++) {
                  S2(j, k) += exp(BZt(L, 0)) * ZT(ks, j, L) * ZT(ks, k, L);
                }
              }
            }
          }

          if (S0(ks, i) > 0.00000001) {
            for (int j = 0; j < NP; j++) {
              U(j) += (ZT(ks, j, i) - S1(j, i) / S0(ks, i)) * (WGHT(ks, i) * DELTA(ks, i) + (1.0 - WGHT(ks, i)) * CENSOR(ks, i) * DRHOipw(ks, i));

              for (int k = 0; k < NP; k++) {
                F(j, k) += (S2(j, k) / S0(ks, i) - S1(j, i) * S1(k, i) / pow(S0(ks, i), 2)) * (WGHT(ks, i) * DELTA(ks, i) + (1.0 - WGHT(ks, i)) * CENSOR(ks, i) * DRHOipw(ks, i));
                F2(j, k) += (ZT(ks, j, i) - S1(j, i) / S0(ks, i)) * (ZT(ks, k, i) - S1(k, i) / S0(ks, i)) * pow((WGHT(ks, i) * DELTA(ks, i) + (1.0 - WGHT(ks, i)) * CENSOR(ks, i) * DRHOipw(ks, i)), 2);
              }
            }
          }
        }
      }
    }
    Kiter = Kiter +1;
    change = arma::solve(F, U);
    BETA += change;
  }

  invF = arma::inv(F);
  var = invF * F2 * invF;
  
  
  for (int ks = 0; ks < KK; ks++) {
    for (int i = 0; i < N(ks); i++) {
      double TEMPB = 0.0;
      for (int ii = 0; ii < N(ks); ii++) {
        if (X(ks, ii) <= tau && S0(ks, ii) > 0.00000001) {
          TEMPB += Epankercplusplus(X(ks, ii), X(ks, i), TBAND, CENSOR(ks, ii)) * (WGHT(ks, ii) * DELTA(ks, ii) + (1.0 - WGHT(ks, ii)) * CENSOR(ks, ii) * DRHOipw(ks, ii)) / S0(ks, ii);
        }
      }
      LAMBDA0(ks, i) = TEMPB;
    }
  }
  
  if(estBaseLamInd == 1){
    for (int ks = 0; ks < KK; ks++) {
      for (int Itgrid = 0; Itgrid < ntgrid; Itgrid++) {
        double tvalue = tstep * Itgrid;
        double TEMPB = 0.0;
        for (int ii = 0; ii < N(ks); ii++) {
          if (X(ks, ii) <= tau) {
            if (S0(ks, ii) > 0.00000001) {
              TEMPB += Epankercplusplus(X(ks, ii), tvalue, TBAND, CENSOR(ks, ii)) * (WGHT(ks, ii) * DELTA(ks, ii) / S0(ks, ii) + (1.0 - WGHT(ks, ii)) * CENSOR(ks, ii) * DRHOipw(ks, ii) / S0(ks, ii));
            }
          }
        }
        LAMBDAk0(ks, Itgrid) = TEMPB;
      }
    } 
  }
  
  return Rcpp::List::create(
    Rcpp::Named("BETA") = BETA,
    Rcpp::Named("F") = invF,
    Rcpp::Named("var") = var,
    Rcpp::Named("LAMBDA0") = LAMBDA0,
    Rcpp::Named("LAMBDAk0") = LAMBDAk0
    );
}

// [[Rcpp::export]]
Rcpp::List AsigInvTempaug(arma::cube SigmaInv, arma::cube G, arma::cube LAMBDAUG, arma::mat CLAMBDAUG,
                          arma::mat wght, arma::mat time, arma::mat censor, arma::mat markm,
                          int ncov, int nvgrid, int kk, arma::vec nsamp, double vstep, double hband,
                          int iskip, double tau) {
  
  arma::cube AsigInv(ncov * ncov, nvgrid, nvgrid, arma::fill::zeros);
  arma::cube tempaug(kk, arma::max(nsamp), nvgrid, arma::fill::zeros);
  
  for (int ispot = iskip; ispot <= nvgrid; ispot++) {
    double va = vstep * (ispot - 1.0);
    double vb = vstep * ispot;
    double vspot = ispot * vstep;
    arma::mat sigtemp(ncov, ncov, arma::fill::zeros);
    
    for (int mspot = iskip; mspot <= nvgrid; mspot++) {
      double vmspot = mspot * vstep;
      sigtemp += SigmaInv.slice(mspot - 1) * Epankercplusplus(vspot, vmspot, hband, 1.0) * vstep;
      AsigInv.subcube(0, ispot - 1, mspot - 1, ncov*ncov-1, ispot - 1, mspot - 1) = arma::vectorise(sigtemp);
    }
    
    arma::mat DRHOaug(kk, arma::max(nsamp), arma::fill::zeros);
    
    for (int ks = 0; ks < kk; ks++) {
      arma::uvec valid_indices = arma::find((time.row(ks) <= tau) && (censor.row(ks) >= 0.5) && (CLAMBDAUG.row(ks) >= 0.00000001));
      arma::rowvec Grow(arma::max(nsamp));
      arma::rowvec LAMBDAUGrow(arma::max(nsamp));
      arma::rowvec DRHOaugrow(arma::max(nsamp));
      DRHOaugrow.zeros();
      for (int ii = 0; ii < arma::max(nsamp); ii++){
        Grow(ii) = G(ks, ii, ispot - 1);
        LAMBDAUGrow(ii) = LAMBDAUG(ks, ii, ispot - 1);
      }
      arma:: rowvec CLAMBDAUGrow = CLAMBDAUG.row(ks);
      DRHOaugrow.elem(valid_indices) = Grow.elem(valid_indices) % LAMBDAUGrow.elem(valid_indices) / CLAMBDAUGrow.elem(valid_indices);
      
      arma::rowvec temp1 = (1.0 - wght.row(ks)) % DRHOaugrow * vstep;
      temp1.elem(arma::find((time.row(ks) > tau) || (censor.row(ks) < 0.5))).zeros();
      
      arma::rowvec temp2 = wght.row(ks);
      temp2.elem(arma::find((time.row(ks) > tau) || (markm.row(ks) < va) || (markm.row(ks) > vb) || (censor.row(ks) <= 0.5))).zeros();
      tempaug.subcube(ks, 0, ispot - 1, ks, max(nsamp)-1, ispot - 1) = temp1 + temp2;
      
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("AsigInv") = AsigInv,
                            Rcpp::Named("tempaug") = tempaug);
}


// [[Rcpp::export]]
Rcpp::List S0NS1N(arma::cube covart, arma::mat beta, arma::mat time, arma::uvec nsamp,
                  int kk, int nvgrid, int ncov, int iskip) {
  
  int mxn = arma::max(nsamp);
  arma::cube S0N(kk, mxn, nvgrid, arma::fill::zeros);
  arma::cube S1N(kk * ncov, mxn, nvgrid, arma::fill::zeros);
  arma::vec BZt(mxn);
  int s = 0;
  for (int ispot = iskip - 1; ispot < nvgrid; ispot++) {
    
    for (int ks = 0; ks < kk; ks++) {
      arma::mat covartks(ncov, mxn);
      covartks = covart.subcube(ks, 0, 0, ks, ncov - 1, mxn - 1);
      BZt = arma::trans(covartks) * beta.col(ispot);
      for (int i = 0; i < nsamp(ks); i++) {
        arma::uvec valid_indices = arma::find(time.row(ks) >= time(ks, i));
        S0N(ks, i, ispot) = arma::sum(arma::exp(BZt(valid_indices)));
        s = ks *ncov;
        for (int j = 0; j < ncov; j++) {
          arma::vec covartks = covart.subcube(ks, j, 0, ks, j, nsamp(ks) - 1);
          S1N(s, i, ispot) = arma::sum(covartks.elem(valid_indices) % arma::exp(BZt.elem(valid_indices)));
          s = s + 1;
        }
      }
      
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("S0N") = S0N,
                            Rcpp::Named("S1N") = S1N);
}

// [[Rcpp::export]]
arma::mat GDIST2Ncplusplus(int nvgrid, int iskip, arma::mat zdev, int KK, arma::ivec N, int NP,
                  arma::mat X, arma::cube ZT, arma::mat beta, arma::cube SigmaInv,
                  arma::cube S0N, arma::cube S1N, arma::cube tempaug, arma::cube AsigInv) {
  
  int mxn = max(N);
  arma::vec BZt(mxn, arma::fill::zeros);
  arma::vec Sx0(mxn, arma::fill::zeros);
  arma::mat Sx1(NP, mxn, arma::fill::zeros);
  arma::mat BU1(NP, nvgrid , arma::fill::zeros);
  arma::mat CUMBDIST(NP, nvgrid, arma::fill::zeros);
  
  for (int ispot = iskip - 1; ispot < nvgrid; ispot++) {
    for (int ks = 0; ks < KK; ks++) {
      arma::mat ZTks(NP, mxn);
      ZTks = ZT.subcube(ks, 0, 0, ks, NP - 1, mxn - 1);
      BZt.col(0) = arma::trans(ZTks) * beta.col(ispot);
      for (int i = 0; i < N(ks); i++) {
        Sx0.zeros();
        Sx1.zeros();
        arma::uvec risk_indices = arma::find(X.row(ks) >= X(ks, i));
        arma::rowvec zdevks = zdev.row(ks);
        Sx0(i) = arma::sum(arma::exp(BZt(risk_indices)) % zdevks.elem(risk_indices));
        for(int j = 0; j < NP; j++){
          arma::vec ZTjL = ZT.subcube(ks, j, 0, ks, j, mxn - 1);
          Sx1(j, i) = Sx1(j, i) + arma::sum(arma::exp(BZt(risk_indices)) % zdevks.elem(risk_indices) % ZTjL(risk_indices));
        }
        
        //Rcpp::Rcout << "2 " << std::endl;
        if (S0N(ks, i, ispot) > 0.00000001) {
          for (int j = 0; j < NP; j++) {
            arma::mat AsigInvk = AsigInv.subcube(j*NP, ispot, iskip - 1, (j+1)*NP - 1, ispot, nvgrid - 1);
            arma::vec S1Nk = S1N.subcube(ks*NP, i, ispot, (ks+1)*NP - 1, i, ispot);
            arma::vec ZTk = ZT.subcube(ks, 0, i, ks, NP - 1, i);
            arma::vec Sx1k = Sx1.col(i);
            arma::vec tempBU1 = arma::trans(AsigInvk) * (ZTk - S1Nk / S0N(ks, i, ispot));
            arma::vec tempXBU2 = arma::trans(AsigInvk) * (Sx1k - Sx0(i) * S1Nk / S0N(ks, i, ispot)) / S0N(ks, i, ispot);

            CUMBDIST.submat(j, iskip - 1, j, nvgrid - 1) += arma::trans((tempBU1 * zdev(ks, i) - tempXBU2) * tempaug(ks, i, ispot));
           
          }
        }
      }
    }//loop ks
  }//loop ispot
  return CUMBDIST;
}



