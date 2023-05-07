#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);                                 // Data vector transmitted from R
  DATA_VECTOR(sex);
  
  DATA_FACTOR(subjectId_factor);                  // Data vector transmitted from R
  DATA_FACTOR(subjectId_day_factor);                  // Data vector transmitted from R
  DATA_FACTOR(subjectId_day_factor_gamma); 
  
  PARAMETER_VECTOR(u);                             // Random effects
  PARAMETER_VECTOR(v);                             // Random effects
  PARAMETER_VECTOR(gamma);                        // Random effects
  PARAMETER_VECTOR(beta);                          // Parameter value transmitted from R
  
  PARAMETER(alpha);
  PARAMETER(sigma2_u_log);
  PARAMETER(sigma2_v_log);
  PARAMETER(sigma2_G_log);
  PARAMETER(sigma2_log);
  
  int nobs = y.size();
  int nsubjects = u.size();
  int ndays = v.size();
  
  Type mean_random_subject = Type(0);
  Type mean_random_day = Type(0);
  Type mean_random_gamma = Type(0);
  
  int i;
  int j;
  
  Type f = 0;                     // Declare the "objective function" (neg. log. likelihood)
  
  for(int i=0; i < nsubjects; i++){
    f -= dnorm(u[i], mean_random_subject,
               sqrt(exp(sigma2_u_log) *exp(alpha * sex[index]) * exp(-gamma[i])),
              true);
  }
  
  for(int j=0; j < ndays; j++){
    i = subjectId_day_factor_gamma[j];
    f -= dnorm(v[j], mean_random_day,
               sqrt(exp(sigma2_v_log) *exp(alpha * sex[index]) * exp(-gamma[i])),
               true);
  }
  
  for(int i=0; i < nsubjects; i++){
    f -= dnorm(gamma[i], mean_random_gamma, sqrt(exp(sigma2_G_log)), true);
  }
  
  // is the indexing correctly here? C++ starts with 0 ...
  for(int index = 0; index < nobs; index++){
    i = subjectId_factor[index];
    j = subjectId_day_factor[index];
    // mu = beta[0] + beta[1]*sex[i] + u[j];
    f -= dnorm(y[index], (beta[0] + beta[1]*sex[index] + u[i] + v[j]),
               sqrt(exp(sigma2_log) *exp(alpha * sex[index]) * exp(-gamma[i])),
               true);
  }
  
  return f;
}