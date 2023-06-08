#include <TMB.hpp>
using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_VECTOR(sex);
  DATA_INTEGER(nsubjects);
  DATA_FACTOR(subjectId_frequencies);
  
  PARAMETER_VECTOR(beta);
  PARAMETER(sigma2_u_log);
  PARAMETER(sigma2_log);
  
  using namespace density;
  int subject_length_i;
  int subject_original_indexing = 0;
  
  Type f = 0;
  
  for (int index = 0; index < nsubjects; index++) {
    subject_length_i = subjectId_frequencies[index];
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Sigma(subject_length_i, subject_length_i);
    Sigma.setZero();
    Eigen::Matrix<Type, Eigen::Dynamic, 1> values_subject_residuals_i(subject_length_i);
    
    for (int row = 0; row < subject_length_i; row++) {
      Sigma(row, row) += exp(sigma2_log);
      
      values_subject_residuals_i[row] =
        y[subject_original_indexing] -
        (beta[0] + beta[1] * sex[subject_original_indexing]);
      subject_original_indexing += 1;
      
      for (int column = 0; column < subject_length_i; column++) {
        Sigma(row, column) += exp(sigma2_u_log);
      }
    }
    f += MVNORM_t<Sigma_type>::pdf(values_subject_residuals_i, Sigma);
  }
  return f;
}
