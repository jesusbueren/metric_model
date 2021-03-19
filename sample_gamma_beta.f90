subroutine sample_gamma(type_i,y_star,gamma)
    use global_var
    implicit none
    integer,dimension(indv,1),intent(in)::type_i
    real(DP),dimension(indv,habits,generations),intent(in)::y_star
    real(DP),dimension(covariates_habits,habits,types),intent(out)::gamma
    integer::e_l,h_l,c_l
    real(DP),dimension(covariates_habits,1)::z
    real(DP),dimension(covariates_habits,covariates_habits)::Sigma,inv_Sigma,A
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    do e_l=1,types; do h_l=1,habits
        do c_l=1,covariates_habits
            z(c_l,1)=c4_normal_01(  )
        end do
        Sigma=matmul(transpose(big_X(1:counter_big_X(e_l,h_l),e_l,h_l,:)),big_X(1:counter_big_X(e_l,h_l),e_l,h_l,:))
        call inverse(Sigma,inv_Sigma,covariates_habits)
        A=inv_Sigma
        call choldc(A,covariates_habits)
        gamma(:,h_l,e_l)=matmul(inv_Sigma,matmul(transpose(big_X(1:counter_big_X(e_l,h_l),e_l,h_l,:)),big_Y(1:counter_big_X(e_l,h_l),e_l,h_l)))+matmul(A,z(:,1))
    end do; end do
        
end subroutine
    
subroutine sample_beta(beta)
    use global_var
    implicit none
    real(DP),dimension(covariates,clusters,clusters+1),intent(out)::beta
    integer::h_l2,h_l,c_l
    real(DP),dimension(covariates,1)::z
    real(DP),dimension(covariates,covariates)::Sigma,inv_Sigma,A
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface

    beta=0.0_dp
    do h_l=1,clusters;do h_l2=1,clusters
        do c_l=1,covariates
            z(c_l,1)=c4_normal_01(  )
        end do
        Sigma=matmul(transpose(big_X_h(1:counter_big_X_h(h_l),h_l,:)),big_X_h(1:counter_big_X_h(h_l),h_l,:))
        call inverse(Sigma,inv_Sigma,covariates)
        A=inv_Sigma
        call choldc(A,covariates)
        beta(:,h_l,h_l2)=matmul(inv_Sigma,matmul(transpose(big_X_h(1:counter_big_X_h(h_l),h_l,:)),big_Y_h(1:counter_big_X_h(h_l),h_l,h_l2)))+matmul(A,z(:,1))
    end do;end do
        
end subroutine    
    
subroutine inverse(a,c,n)
use nrtype
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
real(DP) :: a(n,n), c(n,n)
real(DP):: L(n,n), U(n,n), b(n), d(n), x(n)
real(DP) :: coeff
integer :: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0d0
U=0.0d0
b=0.0d0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0d0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0d0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0d0
end do
end subroutine inverse
    
