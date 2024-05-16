#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]

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
  
  for(int j = 0; j < NP; j++){
    BETA(j) = BETA0(j);
  }

  while (arma::sum(change) > 0.00001 && Kiter <= maxit) {
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
              TEMPB += Epankercplusplus(X(ks, ii), tvalue, TBAND, CENSOR(ks, ii)) * (DELTA(ks, ii)*WGHT(ks, ii)/ S0(ks,ii));
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
  while (arma::sum(change) > 0.00001 && Kiter <= maxit) {
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
arma::mat GDIST2Ncplusplus(int nvgrid, int iskip, arma::mat zdev, int KK, arma::ivec N, int NP,
                  arma::mat X, arma::cube ZT, arma::mat betaofv, arma::cube SigmaInv,
                  arma::cube S0N, arma::cube S1N, arma::cube tempaug, arma::cube AsigInv) {
  const int mxn = max(N);

  arma::mat BZt(mxn, 1, arma::fill::zeros);
  arma::vec Sx0(mxn, arma::fill::zeros);
  arma::mat Sx1(NP, mxn, arma::fill::zeros);
  arma::mat BU1(NP, nvgrid , arma::fill::zeros);
  arma::mat CUMBDIST(NP, nvgrid, arma::fill::zeros);

  for (int ispot = iskip - 1; ispot < nvgrid; ispot++) {
    for (int ks = 0; ks < KK; ks++) {
      for (int i = 0; i < N(ks); i++) {
        BZt(i,0) = 0.0;
        for(int j = 0; j < NP; j++){
          BZt(i, 0) = BZt(i, 0) + betaofv(j, ispot)* ZT (ks, j, i);
        }
      }
      //Rcpp::Rcout << BZt << std::endl;

      for (int i = 0; i < N(ks); i++) {
        Sx0(i) = 0.0;
        Sx1.col(i).zeros();

        for (int L = 0; L < N(ks); L++) {
          if (X(ks, L) >= X(ks, i)) {
            Sx0(i) += exp(BZt(L, 0)) * zdev(ks, L);
            for (int j = 0; j < NP; j++) {
              Sx1(j, i) += exp(BZt(L, 0)) * ZT(ks, j, L) * zdev(ks, L);
            }
          }
        }
       // Rcpp::Rcout << "S0(I) "<< Sx0(I) << std::endl;
        //Rcpp::Rcout << "S1(I) "<< Sx1(0, I) << std::endl;

        if (S0N(ks, i, ispot) > 0.00000001) {
          for (int j = 0; j < NP; j++) {
            for (int mspot = iskip - 1; mspot < nvgrid; mspot++) {
              double tempBU1 = 0.0;
              double tempXBU2 = 0.0;
              for(int k = 0; k < NP; k++){
                int index = j*NP + k;
                int indexS1N = ks*NP + k;
                tempBU1 = tempBU1 + AsigInv (index, ispot, mspot)*(ZT(ks, k, i) - S1N(indexS1N, i, ispot)/ S0N(ks, i, ispot));

                tempXBU2 = tempXBU2 + AsigInv (index, ispot, mspot)*(Sx1(k, i) - Sx0(i) * S1N(indexS1N, i, ispot) / S0N(ks, i, ispot)) / S0N(ks, i, ispot);
              }


              CUMBDIST(j, mspot) += (tempBU1 * zdev(ks, i) - tempXBU2) * tempaug(ks, i, ispot);
              //Rcpp::Rcout << "tempBU1 "<< tempBU1 << std::endl;


            }
          }
        }
      }
    }
  }


  //CUMBDIST = BU1;

  return CUMBDIST;
}
