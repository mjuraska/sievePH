#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double EpankerCplusplus(double tk, double tvalue, double hband, double delt) {
  double result = 0.0;
  
  if (std::abs(tk - tvalue) < hband) {
    result = 3.0 / (4.0 * hband) * delt * (1 - std::pow((tk - tvalue) / hband, 2));
  }
  
  return result;
}


// [[Rcpp::export]]
Rcpp::List estpcomcplusplus(double tau, int KK, arma::ivec N, int NP,
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
  for (int KITER = 1; KITER <= KL; KITER++) {
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
Rcpp::List estpipwcplusplus(double tau, int KK, arma::ivec N, double TBAND, int NP,
                   arma::mat X, arma::cube ZT, arma::mat CENSOR, arma::mat DELTA,
                   arma::mat WGHT, arma::vec BETA0) {
  const int mxn = max(N);
  
  arma::mat BZt(mxn, 1, arma::fill::zeros);
  arma::mat S0(mxn, 1, arma::fill::zeros);
  arma::mat S1(NP, mxn, arma::fill::zeros);
  arma::mat S2(NP, NP, arma::fill::zeros);
  arma::vec U(NP, arma::fill::zeros);
  arma::mat F(NP, NP, arma::fill::zeros);
  arma::mat invF(NP, NP, arma::fill::zeros);
  arma::mat F2(NP, NP, arma::fill::zeros);
  arma::mat var(NP, NP, arma::fill::zeros);
  arma::mat LAMBDA0(KK, mxn, arma::fill::zeros);
  arma::vec BETA (NP);
  
  for(int j = 0; j < NP; j++){
    BETA(j) = BETA0(j);
  }
  
  int KL = 6;
  
  for (int KITER = 1; KITER <= KL; KITER++) {
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
          S0(i, 0) = 0.0;
          S1.col(i).zeros();
          S2.zeros();
          
          for (int l = 0; l < N(ks); l++) {
            if (X(ks, l) >= X(ks, i)) {
              S0(i, 0) += exp(BZt(l, 0)) * WGHT(ks, l);
              
              for (int j = 0; j < NP; j++) {
                S1(j, i) += exp(BZt(l, 0)) * ZT(ks, j, l) * WGHT(ks, l);
                
                for (int k = 0; k < NP; k++) {
                  S2(j, k) += exp(BZt(l, 0)) * ZT(ks, j, l) * ZT(ks, k, l) * WGHT(ks, l);
                }
              }
            }
          }
          
          if (S0(i, 0) > 0.00000001) {
            for (int j = 0; j < NP; j++) {
              U(j) += DELTA(ks, i) * (ZT(ks, j, i) - S1(j, i) / S0(i, 0)) * WGHT(ks, i);
              
              for (int k = 0; k < NP; k++) {
                F(j, k) += DELTA(ks, i) * (S2(j, k) / S0(i, 0) - S1(j, i) * S1(k, i) / pow(S0(i, 0), 2)) * WGHT(ks, i);
                F2(j, k) += pow(DELTA(ks, i), 2) * (ZT(ks, j, i) - S1(j, i) / S0(i, 0)) * (ZT(ks, k, i) - S1(k, i) / S0(i, 0)) * pow(WGHT(ks, i), 2);
              }
            }
          }
        }
      }
      
      if (KITER == KL) {
        for (int i = 0; i < N(ks); i++) {
          double TEMPB = 0.0;
          for (int ii = 0; ii < N(ks); ii++) {
            if (X(ks, ii) <= tau && S0(ii, 0) > 0.00000001) {
              TEMPB += EpankerCplusplus(X(ks, ii), X(ks, i), TBAND, CENSOR(ks, ii)) * DELTA(ks, ii) * WGHT(ks, ii) / S0(ii, 0);
            }
          }
          LAMBDA0(ks, i) = TEMPB;
        }
      }
    }
    
    BETA += arma::solve(F, U);
  }
  
  var = arma::inv(F) * F2 * arma::inv(F);
  invF  = arma::inv(F);
  return Rcpp::List::create(
    Rcpp::Named("BETA") = BETA,
    Rcpp::Named("F") = invF,
    Rcpp::Named("var") = var,
    Rcpp::Named("LAMBDA0") = LAMBDA0
  );
}



// [[Rcpp::export]]
Rcpp::List estpaugcplusplus(double tau, double tstep, int ntgrid, double TBAND, int KK, arma::ivec N, int NP,
             arma::mat X, arma::cube ZT, arma::mat CENSOR, arma::mat DELTA,
             arma::mat WGHT, arma::mat DRHOipw, arma::vec BETA0) {
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
  arma::mat LAMBDAk0(KK,ntgrid, arma::fill::zeros);
  arma::mat F(NP, NP, arma::fill::zeros);
  arma::mat invF(NP, NP, arma::fill::zeros);
  arma::mat var(NP, NP, arma::fill::zeros);
  
  int KL = 6;
  for(int j = 0; j < NP; j++){
    BETA(j) = BETA0(j);
  }
  for (int KITER = 1; KITER <= KL; KITER++) {
    F.zeros();
    F2.zeros();
    S0.zeros();
    U.zeros();
    for (int ks = 0; ks < KK; ks++) {
      for (int I = 0; I < N(ks); I++) {
        BZt.row(I).zeros();
        for (int J = 0; J < NP; J++) {
          BZt(I, 0) += BETA(J) * ZT(ks, J, I);
        }
      }
      
      for (int I = 0; I < N(ks); I++) {
        if (X(ks, I) <= tau) {
          S0(ks, I) = 0.0;
          S1.col(I).zeros();
          S2.zeros();
          
          for (int L = 0; L < N(ks); L++) {
            if (X(ks, L) >= X(ks, I)) {
              S0(ks, I) += exp(BZt(L, 0));
              for (int J = 0; J < NP; J++) {
                S1(J, I) += exp(BZt(L, 0)) * ZT(ks, J, L);
                for (int K = 0; K < NP; K++) {
                  S2(J, K) += exp(BZt(L, 0)) * ZT(ks, J, L) * ZT(ks, K, L);
                }
              }
            }
          }
          
          if (S0(ks, I) > 0.00000001) {
            for (int J = 0; J < NP; J++) {
              U(J) += (ZT(ks, J, I) - S1(J, I) / S0(ks, I)) * (WGHT(ks, I) * DELTA(ks, I) + (1.0 - WGHT(ks, I)) * CENSOR(ks, I) * DRHOipw(ks, I));
             
              for (int K = 0; K < NP; K++) {
                F(J, K) += (S2(J, K) / S0(ks, I) - S1(J, I) * S1(K, I) / pow(S0(ks, I), 2)) * (WGHT(ks, I) * DELTA(ks, I) + (1.0 - WGHT(ks, I)) * CENSOR(ks, I) * DRHOipw(ks, I));
                F2(J, K) += (ZT(ks, J, I) - S1(J, I) / S0(ks, I)) * (ZT(ks, K, I) - S1(K, I) / S0(ks, I)) * pow((WGHT(ks, I) * DELTA(ks, I) + (1.0 - WGHT(ks, I)) * CENSOR(ks, I) * DRHOipw(ks, I)), 2);
              }
            }
          }
        }
      }
      
      if (KITER == KL) {
        for (int I = 0; I < N(ks); I++) {
          double TEMPB = 0.0;
          for (int II = 0; II < N(ks); II++) {
            if (X(ks, II) <= tau) {
              if (S0(ks, II) > 0.00000001) {
                TEMPB += EpankerCplusplus(X(ks, II), X(ks, I), TBAND, CENSOR(ks, II)) * (WGHT(ks, II) * DELTA(ks, II) + (1.0 - WGHT(ks, II)) * CENSOR(ks, II) * DRHOipw(ks, II)) / S0(ks, II);
              }
            }
          }
          LAMBDA0(ks, I) = TEMPB;
        }
      }
    }
    
    BETA += arma::solve(F, U);
  }
  
  invF = arma::inv(F);
  var = invF * F2 * invF;
  
  for (int ks = 0; ks < KK; ks++) {
    for (int Itgrid = 1; Itgrid <= ntgrid; Itgrid++) {
      double tvalue = tstep * Itgrid;
      double TEMPB = 0.0;
      for (int I = 0; I < N(ks); I++) {
        if (X(ks, I) <= tau) {
          if (S0(ks, I) > 0.00000001) {
            TEMPB += EpankerCplusplus(X(ks, I), tvalue, TBAND, CENSOR(ks, I)) * (WGHT(ks, I) * DELTA(ks, I) / S0(ks, I) + (1.0 - WGHT(ks, I)) * CENSOR(ks, I) * DRHOipw(ks, I) / S0(ks, I));
          }
        }
      }
      LAMBDAk0(ks, Itgrid - 1) = TEMPB;
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
      for (int I = 0; I < N(ks); I++) {
        BZt(I,0) = 0.0;
        for(int J = 0; J < NP; J++){
          BZt(I, 0) = BZt(I, 0) + betaofv(J, ispot)* ZT (ks, J, I);
        }
      }
      //Rcpp::Rcout << BZt << std::endl;
      
      for (int I = 0; I < N(ks); I++) {
        Sx0(I) = 0.0;
        Sx1.col(I).zeros();
        
        for (int L = 0; L < N(ks); L++) {
          if (X(ks, L) >= X(ks, I)) {
            Sx0(I) += exp(BZt(L, 0)) * zdev(ks, L);
            for (int J = 0; J < NP; J++) {
              Sx1(J, I) += exp(BZt(L, 0)) * ZT(ks, J, L) * zdev(ks, L);
            }
          }
        }
       // Rcpp::Rcout << "S0(I) "<< Sx0(I) << std::endl;
        //Rcpp::Rcout << "S1(I) "<< Sx1(0, I) << std::endl;
        
        if (S0N(ks, I, ispot) > 0.00000001) {
          for (int J = 0; J < NP; J++) {
            for (int mspot = iskip - 1; mspot < nvgrid; mspot++) {
              double tempBU1 = 0.0;
              double tempXBU2 = 0.0;
              for(int K = 0; K < NP; K++){
                int index = J*NP + K;
                int indexS1N = ks*NP+K;
                tempBU1 = tempBU1 + AsigInv (index, ispot, mspot)*(ZT(ks, K, I) - S1N(indexS1N, I, ispot)/ S0N(ks, I, ispot));
                
                tempXBU2 = tempXBU2 + AsigInv (index, ispot, mspot)*(Sx1(K, I) - Sx0(I) * S1N(indexS1N, I, ispot) / S0N(ks, I, ispot)) / S0N(ks, I, ispot);
              }
             
              
              CUMBDIST(J, mspot) += (tempBU1 * zdev(ks, I) - tempXBU2) * tempaug(ks, I, ispot);
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
