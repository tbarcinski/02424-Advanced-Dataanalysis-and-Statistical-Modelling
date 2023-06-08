#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);                                 // Data vector transmitted from R
  DATA_VECTOR(sex);
  
  DATA_FACTOR(subjectId_factor);                  // Data vector transmitted from R
  
  PARAMETER_VECTOR(u);                             // Random effects
  PARAMETER_VECTOR(beta);                          // Parameter value transmitted from R

  int nobs = y.size();
  int nsubjects = u.size();
  Type mean_random_effects = Type(0);
  
  int j;
  
  Type f = 0;                     // Declare the "objective function" (neg. log. likelihood)
  
  for(int j=0; j < nsubjects; j++){
    f -= dnorm(u[j], mean_random_effects, sqrt(exp(beta[2])), true);
  }
  
  for(int i =0; i < nobs; i++){
    j = subjectId_factor[i];
    // mu = beta[0] + beta[1]*sex[i] + u[j];
    f -= dnorm(y[i], (beta[0] + beta[1]*sex[i] + u[j]), sqrt(exp(beta[3])), true);
  }
  
  return f;
}