#include <RcppArmadillo.h>
  
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat schur(arma::mat a, arma::mat b, int dim) { 
  arma::mat res(dim, dim);
  res.fill(0);
  for(int i=0; i<dim;i++){
    for(int j=0; j<dim;j++){
      res(i,j) = a(i,j) * b(i,j);
    }
  }
  return(res);
}

// [[Rcpp::export]]
double mean_cpp(NumericVector x){
  int n = x.size();
  double sum = 0;
  for(int i = 0; i < n; i++){
    sum += x[i];
  }
  return sum/n;
}

// [[Rcpp::export]]
arma::mat covariance_cpp(NumericMatrix M){
  int dim = M.nrow();
  int N = M.ncol();
  arma::mat R(dim, dim);
  R.fill(0);
  NumericVector mean_vec(dim);
  for(int i=0; i<dim; i++){
    mean_vec[i] = mean_cpp(M(i,_));
  }
  for(int i=0; i<N; i++){
    R = R + as<arma::vec>(wrap(M(_,i)-mean_vec)) * (as<arma::vec>(wrap(M(_,i)-mean_vec))).t();
  }
  return (R / (N-1));
}

// [[Rcpp::export]]
NumericVector seq_cpp(int from, int to, int by, Function f) {
  NumericVector res = f(from, to, by);
  return res;
}

// [[Rcpp::export]]
NumericVector solve_ctrimat_cpp(NumericVector a, NumericVector b,
                                NumericVector c, NumericVector d) {
  int a_size = a.size();
  int c_size = c.size();
  int d_size = d.size();
  NumericVector solution(a_size+1);
  NumericVector c_pr(c_size);
  NumericVector d_pr(d_size);
  for(int i = 0; i < a_size; i++){
    if(i == 0){
      c_pr[i] = c[i]/b[i];
      d_pr[i] = d[i]/b[i];
    }else{
      c_pr[i] = c[i]/(b[i] - a[i-1]*c_pr[i-1]);
      d_pr[i] = (d[i] - a[i-1]*d_pr[i-1])/(b[i]-a[i-1]*c_pr[i-1]);
    }
  }
  
  int n = a_size;
  d_pr[n] = (d[n] - a[n-1]*d_pr[n-1])/(b[n]-a[n-1]*c_pr[n-1]);
  solution[n] = d_pr[n];
  for(int i = n-1; i>-1; i--){
    solution[i] = d_pr[i] - c_pr[i]*solution[i+1];
  }
  return solution;
}

// [[Rcpp::export]]
NumericVector solve_system_cpp(double a_0, NumericVector a, NumericVector b,
                               NumericVector c, double c_0, NumericVector d){
  
  int dim = d.size();
  NumericVector arg_1(dim-2);
  NumericVector arg_2(dim-1);
  NumericVector arg_3(dim-2);
  NumericVector arg_4(dim-1);
  for(int i=0; i<dim-1; i++){
    arg_2[i] = b[i];
    arg_4[i] = d[i];
  }
  for(int i=0; i<dim-2; i++){
    arg_1[i] = a[i];
    arg_3[i] = c[i];
  }
  NumericVector X_1 = solve_ctrimat_cpp(arg_1,arg_2,arg_3,arg_4);
  
  arg_4[0] = -c_0;
  arg_4[dim-2] = -c[dim-2];
  for(int i=1;i < dim-2; i++){
    arg_4[i] = 0;
  }
  NumericVector X_2 = solve_ctrimat_cpp(arg_1,arg_2,arg_3,arg_4);
  
  
  double x_last = (d[dim-1]-a_0*X_1[0]-a[dim-2]*X_1[dim-2])/(a_0*X_2[0]+a[dim-2]*X_2[dim-2]+b[dim-1]);
  NumericVector result(dim);
  for(int i=0;i<dim-1;i++){
    result[i] = X_1[i] + x_last*X_2[i];
  }
  result[dim-1] = x_last;
  return result;
}

// [[Rcpp::export]]
NumericVector for_1_cpp(NumericVector x){
  int size = x.size();
  double temp = x[0];
  for(int i=0; i<size-1; i++){
    x[i] = x[i+1];
  }
  x[size-1] = temp;
  return x;
}

// [[Rcpp::export]]
NumericVector rec_1_cpp(NumericVector x){
  int size = x.size();
  double temp = x[size-1];
  for(int i=size-1; i>0; i--){
    x[i] = x[i-1];
  }
  x[0] = temp;
  return x;
}


// [[Rcpp::export]]
NumericVector predict_cpp(int dim, double delta_x, double delta_t, NumericVector u,
                          NumericVector nu, NumericVector rho, 
                          NumericVector analysis, NumericVector noise) {
  
  double a_0; // low left corner
  if(u[dim-1] >= 0){
    a_0 = -nu[dim-1]/(delta_x*delta_x);
  }else{
    a_0 = u[dim-1]/delta_x-nu[dim-1]/(delta_x*delta_x);
  }
  
  NumericVector a(dim-1);
  for(int i=1;i<dim; i++){
    if(u[i] >= 0){
      a[i-1] = -u[i]/delta_x - nu[i]/(delta_x*delta_x);
    }else{
      a[i-1] = - nu[i]/(delta_x*delta_x);
    }
    
  }
  
  NumericVector b(dim);
  NumericVector d(dim);
  for(int i=0;i<dim; i++){
    if(u[i] >= 0){
      b[i] = 1/delta_t + rho[i] + 2*nu[i]/(delta_x*delta_x) + u[i]/delta_x;
    }else{
      b[i] = 1/delta_t + rho[i] + 2*nu[i]/(delta_x*delta_x) - u[i]/delta_x;
    }
  }
  
  NumericVector c(dim-1);
  for(int i=0;i<dim-1; i++){
    if(u[i] >= 0){
      c[i] = - nu[i]/(delta_x*delta_x);
    }else{
      c[i] = u[i]/delta_x - nu[i]/(delta_x*delta_x);
    }
  }
  double c_0; // upper rigth corner
  if(u[0] >= 0){
    c_0 = -u[0]/delta_x-nu[0]/(delta_x*delta_x);
  }else{
    c_0 = -nu[0]/(delta_x*delta_x);
  }
  
  
  
  for(int i=0;i<dim; i++){
    d[i] = analysis[i]/delta_t + noise[i];
  }
  
  NumericVector result = solve_system_cpp(a_0,a,b,c,c_0,d);
  return result;
}

// [[Rcpp::export]]
NumericVector noise_cpp(arma::mat X, arma::mat Y, arma::vec z) {

  arma::vec result = (X * Y * z);

  return wrap(result);
}

NumericVector double_mult_cpp(arma::mat X, arma::mat Y, arma::vec z) {

  arma::mat result = (X * Y * z);
  return wrap(result);
}


// [[Rcpp::export]]
NumericMatrix generate_field_cpp(int time, int dim, NumericVector start_X, double delta_t, 
                                 NumericMatrix U, NumericMatrix RHO, NumericMatrix Sigma, NumericMatrix Nu, 
                                 double L_C, Function create_cov_matrix){
  double Re = 6370000; 
  double delta_x = 2*PI*Re/dim;
  NumericMatrix C = create_cov_matrix(dim, rep(1,dim), L_C);
  NumericMatrix X(dim, time);
  X(_,0) = start_X;
  NumericMatrix S(dim,dim);
  NumericVector nor(dim);
  NumericVector noise(dim);
  for(int i=1; i<time; i++){
    for(int k=0; k<dim; k++){
      S(k,k) = Sigma(k,i);
    }
    
    nor = rnorm(dim,0, 1/(sqrt(delta_x)*sqrt(delta_t)));
    noise = noise_cpp(as<arma::mat>(S), as<arma::mat>(C), as<arma::vec>(nor));
    X(_,i) = predict_cpp(dim, delta_x, delta_t, U(_,i), Nu(_,i), RHO(_,i), X(_,i-1), noise);
  }
  return X;
}

// [[Rcpp::export]]
int generate(int time, int dim, NumericVector start_X, double delta_t, 
             NumericMatrix U, NumericMatrix RHO, NumericMatrix Sigma, NumericMatrix Nu, double L_C,
             Function create_cov_matrix, int TOTAL){
  NumericMatrix X(dim, time);
  for(int tot=0; tot < TOTAL; tot++){
    X = generate_field_cpp(time, dim, start_X, delta_t, U, RHO, Sigma, Nu, L_C, create_cov_matrix);
  }
  return 1;
}


// [[Rcpp::export]]
double inflation(double diag_1, double diag_2, NumericVector inflation_bounds, NumericVector inflation_coef){
  double result;
  double val1;
  double val2;
  val1 = inflation_coef[1];
  if(diag_1 < inflation_bounds[0]) val1 = inflation_coef[0];
  if(diag_1 > inflation_bounds[1]) val1 = inflation_coef[2];
  
  val2 = inflation_coef[1];
  if(diag_2 < inflation_bounds[0]) val2 = inflation_coef[0];
  if(diag_2 > inflation_bounds[1]) val2 = inflation_coef[2];
  result = val1 * val2;
  return result;
}

// [[Rcpp::export]]
List ensembl_kalman_cpp(int time, int dim, double delta_t, NumericMatrix U, NumericMatrix RHO,
                        NumericMatrix Sigma, NumericMatrix Nu, double L_C, NumericVector R, 
                        NumericVector ind_obs, NumericMatrix OBS, NumericVector start_X,
                        Function create_cov_matrix, int ens_size, NumericVector inflation_bounds,
                        NumericVector inflation_coef, arma::mat C_enkf, int thin_time_for_filer){
  
  NumericMatrix X_a(dim, time/thin_time_for_filer);
  NumericMatrix X_f(dim, time/thin_time_for_filer);
  NumericVector b_arr(time*dim*dim/thin_time_for_filer);
  arma::cube B_arr(b_arr.begin(), dim, dim, time/thin_time_for_filer, false);
  NumericVector s_nonloc_arr(time*dim*dim/thin_time_for_filer);
  arma::cube S_nonloc_arr(s_nonloc_arr.begin(), dim, dim, time/thin_time_for_filer, false);
  
  NumericMatrix x_me_arr(dim, time);
  NumericVector y;
  arma::mat BH(dim, ind_obs.size());
  arma::mat HBH(ind_obs.size(), ind_obs.size());
  arma::vec z;
  arma::vec v;
  X_f(_,0) = start_X;
  X_a(_,0) = start_X;
  int Re = 6370000;  
  double delta_x = 2*PI*Re/dim;
  NumericMatrix C = create_cov_matrix(dim, rep(1,dim), L_C);
  
  //A_arr.slice(0) = as<arma::mat>(start_A); 
  NumericMatrix PHI(dim,dim);
  NumericMatrix A_num(dim,dim);
  NumericVector noise(dim, 0.0);
  NumericMatrix B_1(dim,dim);
  NumericVector nor(dim);
  NumericVector eta(ind_obs.size());
  arma::mat B(dim,dim);
  arma::mat S(dim,dim);
  arma::mat R_mat(ind_obs.size(), ind_obs.size());
  arma::mat A;
  S.fill(0);
  B.fill(0);
  R_mat.fill(0);
  R_mat = arma::diagmat(as<arma::vec>(R));
  
  NumericMatrix X_en_f(dim, ens_size);  //filled with 0
  NumericMatrix temp_cov(dim, ens_size);
  NumericMatrix X_en_a(dim, ens_size);
  NumericMatrix temp1(dim, ens_size);
  NumericMatrix temp2(dim, ens_size);
  NumericMatrix X_me(dim, ens_size);
  NumericMatrix X_pe(dim, ens_size);
  NumericMatrix X_en_f_mean(dim);
  
  NumericMatrix diag(dim, time/thin_time_for_filer);
  NumericMatrix sup2diag(dim-2, time/thin_time_for_filer);
  NumericMatrix sup3diag(dim-3, time/thin_time_for_filer);
  NumericMatrix supdiag(dim-1, time/thin_time_for_filer);
  NumericMatrix Lambda(dim, time/thin_time_for_filer);
  NumericMatrix spectr(dim, time/thin_time_for_filer);
  for(int i=0; i<dim; i++){
    diag(i,0) = B(i,i);
  }

  arma::vec eigval;
  arma::mat eigvec;
  
  eig_sym(eigval, eigvec, B);
  spectr(_,0) = as<NumericVector>(wrap(eigval));
  
  for(int i=0; i<ens_size; i++){
    X_en_a(_,i) = start_X;
    X_en_f(_,i) = start_X;
  }
  X_a(_,0) = start_X;
  X_f(_,0) = start_X;
  
  for(int step = 1; step < time; step ++){
    for(int i=0; i<ens_size; i++){
      X_pe(_,i) = predict_cpp(dim, delta_x, delta_t, U(_,step), Nu(_,step), RHO(_,step), X_en_a(_,i), noise);
    }
    
    for(int k=0; k<dim; k++){
      S(k,k) = Sigma(k,step);
    }
    
    
    
    for(int i=0; i<ens_size; i++){  
      nor = rnorm(dim,0,sqrt(delta_t)/sqrt(delta_x));
      X_me(_,i) = noise_cpp(S, as<arma::mat>(C), as<arma::vec>(nor));
      X_me(_,i) = predict_cpp(dim, delta_x, delta_t, U(_,step), Nu(_,step), RHO(_,step), X_me(_,i), noise);
    }
    
    
    X_en_f = wrap(as<arma::mat>(X_pe) + as<arma::mat>(X_me));
    if((step % thin_time_for_filer)){
      X_en_a = X_en_f;
    }else{
      for(int i=0; i<dim; i++){
        X_en_f_mean[i] = mean_cpp(X_en_f(i,_)); 
      }
    
      X_f(_,step/thin_time_for_filer) = X_en_f_mean;
    
      
      temp1 = X_en_f;
      S_nonloc_arr.slice(step/thin_time_for_filer) = covariance_cpp(X_en_f);
      temp_cov = as<NumericMatrix>(wrap(covariance_cpp(X_en_f)));
      for(int i=0; i<dim; i++){
        for(int j=0; j<ens_size; j++){
          X_en_f(i,j) = X_en_f_mean[i] + (X_en_f(i,j) - X_en_f_mean[i]) * inflation(temp_cov(i,i), temp_cov(j,j), inflation_bounds, inflation_coef);
        }
      }
      temp2 = X_en_f;
    
      eig_sym(eigval, eigvec, covariance_cpp(X_en_f));
      spectr(_,step/thin_time_for_filer) = as<NumericVector>(wrap(eigval));
    
      B = schur(covariance_cpp(X_en_f), C_enkf, dim);
      //B = schur(covariance_cpp(X_en_f), C_enkf, dim) + inflation_beta * B_enkf_mean;
      B_arr.slice(step/thin_time_for_filer) = B;  
      y             = OBS(_,step);              
      BH            = B.cols(as<arma::uvec>(ind_obs));
      HBH           = B.submat(as<arma::uvec>(ind_obs), as<arma::uvec>(ind_obs));
   
      for(int i_en=0; i_en<ens_size;i_en++){
        v = X_en_f(_,i_en);
        v = v.elem(as<arma::uvec>(ind_obs));
        eta = rnorm(ind_obs.size(),0,1);
        for(int k=0; k<(ind_obs.size()); k++){
          eta[k] = eta[k] * sqrt(R[k]);
        }
        z             = arma::solve(HBH + R_mat, as<arma::vec>(y) + as<arma::vec>(eta) - v);
        X_en_a(_,i_en)      = X_en_f(_,i_en) + as<NumericVector>(wrap(BH * z));
      }
    
      for(int i=0; i<dim; i++){
        X_a(i,step/thin_time_for_filer) = mean_cpp(X_en_a(i,_)); 
      }
    
      for(int i=0; i<dim; i++){
        diag(i,step/thin_time_for_filer) = B(i,i);
      }

    
    }

    }
    
    return Rcpp::List::create(Rcpp::Named("X_a") = X_a,
                              Rcpp::Named("X_f") = X_f,
                              Rcpp::Named("B_arr") = B_arr,
                              Rcpp::Named("S_nonloc_arr") = S_nonloc_arr,
                              Rcpp::Named("B") = B,
                              Rcpp::Named("spectr") = spectr,
                              Rcpp::Named("diag") = diag);
}          


