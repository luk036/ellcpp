%
% [L,D]=ldlt(A)
%
% This function computes the square root free Cholesky factorization
%
%    A=L*D*L'
%
% where L is a lower triangular matrix with ones on the diagonal, and D
% is a diagonal matrix.  
%
% It is assumed that A is symmetric and postive definite.
%
% Reference: Golub and Van Loan, "Matrix Computations", second edition, 
%            p 137.  
% Author: Brian Borchers (borchers@nmt.edu)
%
function [L,D]=ldlt(A)
%
%  Figure out the size of A.
%
n=size(A,1);
%
% The main loop.  See Golub and Van Loan for details.  
%
L=zeros(n,n);
for j=1:n,
  if (j > 1),
    v(1:j-1)=L(j,1:j-1).*d(1:j-1);
    v(j)=A(j,j)-L(j,1:j-1)*v(1:j-1)';
    d(j)=v(j);
    if (j < n),
      L(j+1:n,j)=(A(j+1:n,j)-L(j+1:n,1:j-1)*v(1:j-1)')/v(j);
    end;
  else
    v(1)=A(1,1);
    d(1)=v(1);
    L(2:n,1)=A(2:n,1)/v(1);    
  end;
end;
%
%  Put d into a matrix.
%
D=diag(d);
%
%  Put ones on the diagonal of L.
%
L=L+eye(n);