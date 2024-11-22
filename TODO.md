# TODO:

In the evaluation of the e-step in line 27, instead of computing the determinant as "const double Ldet = L.determinant();", 
we can use the pre-computed cholesky factor from line 25 and use the formula: 

log_det == 2 * sum(log(diag(chol(M)))),

where M is a symmetric semidefinite positive matrix.


