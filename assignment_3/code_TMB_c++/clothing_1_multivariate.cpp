#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);                                 // Data vector transmitted from R
  DATA_VECTOR(sex);
  DATA_INTEGER(nsubjects);
  
  DATA_FACTOR(subjectId_frequencies);

  PARAMETER_VECTOR(beta);                          // Parameter value transmitted from R
  PARAMETER(sigma2_u_log);
  PARAMETER(sigma2_log);
  
  using namespace density;
  int subject_length_i;
  int subject_original_indexing = 1;
  
  Type f = 0;                     // Declare the "objective function" (neg. log. likelihood)
  
  for (int index=0;index<nsubjects;index++){
    
    subject_length_i = subjectId_frequencies[index];
    
    matrix<Type> Sigma(subject_length_i,subject_length_i);
    Sigma.setZero();
    vector<Type> values_subject_residuals_i(subject_length_i);

    for (int row=0; row < 10; row++){

      Sigma(row, row) += exp(sigma2_log);
      
      values_subject_residuals_i[row] =
        y[subject_original_indexing] - 
        (beta[0] + beta[1]*sex[subject_original_indexing]);
      subject_original_indexing += 1;

      for (int column=0; column < 10; column++){
        Sigma(row, column) += exp(sigma2_u_log);
      }
    }
    f += MVNORM_t(Sigma)(values_subject_residuals_i);
  }  
  return f;
}
