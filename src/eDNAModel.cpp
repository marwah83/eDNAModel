#define TMB_LIB_INIT R_init_occ
#include <TMB.hpp>
#include<math.h>


enum valid_link{
  log_link = 0,
  logit_link = 1,
  probit_link = 2,
  cloglog_link = 3
};

//modelled after glmmTMB's inverse_linkfun
template<class Type>
Type inverse_link(Type eta, int link){
  Type ans;
  switch (link){
  case log_link:
    eta = exp(eta);
    break;
  case logit_link:
    ans = invlogit(eta);
    break;
  case probit_link:
    ans = pnorm(eta);
    break;
  case cloglog_link:
    ans = 1-exp(-exp(eta));
    break;
  default:
    error("Invalid link.");
  }
  return(ans);
};

//piece of code stolen from gllvm package
template<class Type>
matrix<Type> constructL(const vector<Type> &theta){
  int n = (1 + sqrt(1 + 8 *  theta.size())) / 2;

  matrix<Type> L = Eigen::MatrixXd::Zero(n, n);  // Initialize Cholesky factor
  L.diagonal().fill(1.0);
  int covscounter = 0;
  for (int j = 0; j < n; ++j) {
    for (int i = j + 1; i < n; ++i) {
      L(i, j) = theta[covscounter];
      covscounter++;
    }
  }
  for (int i = 1; i < n; ++i) {
    L.row(i) /= L.row(i).norm();
  }
  return L;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_IMATRIX(Y);
  DATA_IMATRIX(Ysites);
  DATA_MATRIX(offset);
  DATA_MATRIX(Xa);//abundance covs
  DATA_MATRIX(Xo);//occupancy covs
  DATA_SPARSE_MATRIX(Za);//abundance re matrix
  DATA_SPARSE_MATRIX(Zo);//occupancy re matrix
  DATA_IMATRIX(csa);//cor structure for abundance REs
  DATA_IMATRIX(cso);//cor structure for occupancy REs
  DATA_INTEGER(family);
  DATA_MATRIX(NTrials);
  DATA_INTEGER(linka);
  DATA_INTEGER(linko);
  DATA_IVECTOR(sites);//should be ordered and zero-indexed
  PARAMETER_MATRIX(Ba);//abundance fix effs
  PARAMETER_MATRIX(Bo);//occupancy fix effs
  PARAMETER_MATRIX(Ua);//abundance re effs
  PARAMETER_VECTOR(logsda);//abundance re logsds
  PARAMETER_VECTOR(corsa);//abundance re cors
  PARAMETER_MATRIX(Uo);//occupancy re effs
  PARAMETER_VECTOR(logsdo);//occupancy re logsds
  PARAMETER_VECTOR(corso);//occupancy re cors
  PARAMETER_VECTOR(logphi);//for ZINB

  REPORT(Ba);
  REPORT(Bo);

  Type nll = 0;

  //log-mean of abundance process
  //replicate-specific
  matrix<Type>etaa(Xa.rows(),Y.cols());
  etaa.fill(0.0);

  etaa = Xa*Ba;

  if(offset.rows()>0){
  etaa += offset;
  }

  if(Za.rows()==Xa.rows()){
    etaa += Za*Ua;

    //re likelihood component
    //Set-up SDs of REs
    matrix<Type> sdsa = Eigen::MatrixXd::Zero(Za.cols(),Za.cols());
    sdsa.diagonal() = exp(logsda);

    //Set-up correlation parameters
    vector<Type>sigmaija((Za.cols()*Za.cols()-Za.cols())/2);
    sigmaija.fill(0.0);
    matrix<Type>SigmaaL(Za.cols(),Za.cols());//This will be the cholesky of the covariance matrix
    SigmaaL.setIdentity();
    matrix<Type>Sigmaa = SigmaaL;//Create an empty covariance matrix
    if(csa.cols()>1){
      //need a vector with covariances and zeros in the right places
      for(int i=0; i<csa.rows(); i++){
        sigmaija((csa(i,0) - 1) * (csa(i,0) - 2) / 2 + csa(i,1)-1) = corsa(i);
      }
      SigmaaL = sdsa*constructL(sigmaija);//scale the cholesky of the correlation matrix
    }else{
      SigmaaL = sdsa;
    }
      Sigmaa = SigmaaL*SigmaaL.transpose();
      REPORT(Sigmaa);
      REPORT(SigmaaL);

    density::MVNORM_t<Type> mvnorm(Sigmaa);
    for(int s = 0; s<Y.cols(); s++){
        nll += mvnorm(Ua.col(s));
    }
  }

  //mean of occupancy process
  //site-specific
  matrix<Type>etao(Ysites.rows(),Y.cols());
  etao.fill(0.0);
  etao = Xo*Bo;
  if(Zo.rows()==Ysites.rows()){
    etao += Zo*Uo;

    //re likelihood component
    //Set-up SDs of REs
    matrix<Type> sdso = Eigen::MatrixXd::Zero(Zo.cols(),Zo.cols());
    sdso.diagonal() = exp(logsdo);

    //Set-up correlation parameters
    vector<Type>sigmaijo((Zo.cols()*Zo.cols()-Zo.cols())/2);
    sigmaijo.fill(0.0);
    matrix<Type>SigmaoL(Zo.cols(),Zo.cols());//This will be the cholesky of the covariance matrix
    SigmaoL.setIdentity();
    matrix<Type>Sigmao = SigmaoL;//Create an empty covariance matrix
    if(cso.cols()>1){
      //need a vector with covariances and zeros in the right places
      for(int i=0; i<cso.rows(); i++){
        sigmaijo((cso(i,0) - 1) * (cso(i,0) - 2) / 2 + cso(i,1)-1) = corso(i);
      }
      SigmaoL = sdso*constructL(sigmaijo);//scale the cholesky of the correlation matrix
    }else{
      SigmaoL = sdso;
    }
    Sigmao = SigmaoL*SigmaoL.transpose();
    REPORT(Sigmao);

    density::MVNORM_t<Type> mvnorm(Sigmao);
    for(int s = 0; s<Y.cols(); s++){
      nll += mvnorm(Uo.col(s));
    }
  }

  REPORT(etaa);
  REPORT(etao);

  matrix<Type>logProbs(Ysites.rows(), Y.cols());
  logProbs.fill(0.0);

  if(family==0){ // ZIP
  for(int s = 0; s<Y.cols(); s++){
    for(int j = 0; j<Y.rows(); j++){
      logProbs(sites(j),s) += dpois(Type(Y(j,s)), inverse_link(etaa(j,s), linka), true);
    }
  }
  }else if(family==1){ //ZINB
    for(int s = 0; s<Y.cols(); s++){
      for(int j = 0; j<Y.rows(); j++){
        logProbs(sites(j),s) += dnbinom_robust(Type(Y(j,s)), inverse_link(etaa(j,s), linka), 2*etaa(j,s)-logphi(s), true);
      }
    }
    // for(int s = 0; s<Y.cols(); s++){
    //   int site = sites(0);
    //   Type res = 0;
    //   int yreps = 0;
    //   for(int i = 0; i<Xa.rows(); i++){
    //     //when we have switched to a new siten
    //     //we first add the likelihood calculation of the last site
    //     //as it means that we have accumulated the necessary product in the likelihood
    //     if(site!=sites(i)){
    //       if(yreps==0){
    //         //no detection
    //         nll -= log(invlogit(etao(site,s)) + (Type(1)-invlogit(etao(site,s)))*exp(res));
    //       }else{
    //         //detection
    //         nll -= log(Type(1)-invlogit(etao(site,s))) + res;
    //       }
    //       site = sites(i);
    //       yreps = 0;//reset summation
    //       res = 0;//reset product of replicate probabilities
    //     }
    //
    //     //as long as we  are in the same site
    //     // we need to take a product over replicates
    //     // and sum over Y
    //     yreps += Y(i,s);
    //     res += dnbinom_robust(Type(Y(i,s)), exp(etaa(i,s)), 2*etaa(i,s)-logphi(s), true);//product of probabilities
    //   }
    // }
  }else if(family == 2){ // Binomial
    for(int s = 0; s<Y.cols(); s++){
      for(int j = 0; j<Y.rows(); j++){
        logProbs(sites(j),s) += dbinom_robust(Type(Y(j,s)), NTrials(j,s), inverse_link(etaa(j,s), linka), true);
      }
    }
    // for(int s = 0; s<Y.cols(); s++){
    //   int site = sites(0);
    //   Type res = 0;
    //   int yreps = 0;
    //   for(int i = 0; i<Xa.rows(); i++){
    //     //when we have switched to a new siten
    //     //we first add the likelihood calculation of the last site
    //     //as it means that we have accumulated the necessary product in the likelihood
    //     if(site!=sites(i)){
    //       if(yreps==0){
    //         //no detection
    //         nll -= log(invlogit(etao(site,s)) + (Type(1)-invlogit(etao(site,s)))*exp(res));
    //       }else{
    //         //detection
    //         nll -= log(Type(1)-invlogit(etao(site,s))) + res;
    //       }
    //       site = sites(i);
    //       yreps = 0;//reset summation
    //       res = 0;//reset product of replicate probabilities
    //     }
    //
    //     //as long as we  are in the same site
    //     // we need to take a product over replicates
    //     // and sum over Y
    //     yreps += Y(i,s);
    //     res += dbinom_robust(Type(Y(i,s)), NTrials(i,s), invlogit(etaa(i,s)), true);//product of probabilities
    //   }
    // }
  }

  //finalize LL
  for(int s = 0; s<Y.cols(); s++){
    for(int i = 0; i<Ysites.rows(); i++){
      if(Ysites(i,s)==0){
        //no detection
        nll -= log(inverse_link(etao(i,s), linko) + (Type(1)-inverse_link(etao(i,s), linko))*exp(logProbs(i,s)));
      }else{
        //detection
        nll -= log(Type(1)-inverse_link(etao(i,s), linko)) + logProbs(i,s);
      }
    }
  }


  return nll;
}
