#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);                                 // Data vector transmitted from R
  DATA_VECTOR(sex);
  
  DATA_FACTOR(subjectId_frequencies);
  DATA_FACTOR(subjectId_factor);                  // Data vector transmitted from R
  
  PARAMETER_VECTOR(beta);                          // Parameter value transmitted from R
  
  PARAMETER(nsubjects_size);
  PARAMETER(sigma2_u_log);
  PARAMETER(sigma2_log);
  
  int nsubjects = nsubjects_size;
  
  using namespace density;
  int subject_length_i;
  int subject_original_indexing = 0;
  
  Type f = 0;                     // Declare the "objective function" (neg. log. likelihood)
  
  // for(int j=0; j < nsubjects; j++){
  //   f -= dnorm(u[j], mean_random_effects, sqrt(exp(beta[2])), true);
  // }
  // 
  // for(int i =0; i < nobs; i++){
  //   j = subjectId_factor[i];
  //   // mu = beta[0] + beta[1]*sex[i] + u[j];
  //   f -= dnorm(y[i], (beta[0] + beta[1]*sex[i] + u[j]), sqrt(exp(beta[3])), true);
  // }
  
  for (int index=0;index<nsubjects;index++){
    
    subject_length_i = subjectId_frequencies[index];
    
    matrix<Type> Sigma(subject_length_i,subject_length_i);
    vector<Type> values_subject_residuals_i(subject_length_i);
    
    for (int row=0; row < subject_length_i; row++){
      // on the diagonal there are sigma
      Sigma(row,row) += exp(sigma2_log);
      // j = subjectId_factor[subject_original_indexing];
      
      values_subject_residuals_i[row] = y[subject_original_indexing] - 
        (beta[0] + beta[1]*sex[subject_original_indexing]);
      
      subject_original_indexing += row;
      for (int column=0; column < subject_length_i; column++){
        Sigma(row,column) += exp(sigma2_u_log);
      }
    }
    // evaluate negative log likelihood

    f += MVNORM(Sigma)(values_subject_residuals_i);
  }  
  
  return f;
}