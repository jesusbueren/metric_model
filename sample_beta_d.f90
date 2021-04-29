  
subroutine sample_beta_d(beta_d,type_i,sample_k)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates,types,clusters),intent(inout)::beta_d
    integer,dimension(indv,1),intent(in)::type_i
    integer,dimension(indv,generations),intent(in)::sample_k
    integer::h_l,c_l,g_l,ge_l,age,ge_d,it,i_l,health_d,d_l,t_l
    real(DP)::d_star
    double precision,dimension(L_educ-1)::educ_d
    double precision,dimension(types-1)::type_d
    real(DP)::gender_d
    real(DP),dimension(covariates,1)::z
    real(DP),dimension(covariates,covariates)::Sigma,inv_Sigma,A
    real(DP),dimension(indv*g_max,clusters,types,covariates)::big_X_d
    real(DP),dimension(indv*g_max,clusters,types)::big_Y_d
    integer,dimension(clusters,types)::counter_big_X_d
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    real(DP),dimension(covariates,1)::x

     counter_big_X_d=0
    do i_l=1,indv;do g_l=first_age(i_l),last_age(i_l)-1
            x=-9.0d0
            age=initial_age+(g_l-1)*2-70
            gender_d=dble(gender(i_l)-1)
            x(1:4,1)=(/1.0_dp,dble(age),gender_d,dble(age)*gender_d/)
            educ_d=0.0d0
            if (educ(i_l)>1) then
                educ_d(educ(i_l)-1)=1.0d0
            end if
            x(5:6,1)=educ_d
            x(7:8,1)=educ_d*dble(age)
            if (sample_k(i_l,g_l)>=1 .and. sample_k(i_l,g_l+1)>=1) then
                counter_big_X_d(sample_k(i_l,g_l),type_i(i_l,1))=counter_big_X_d(sample_k(i_l,g_l),type_i(i_l,1))+1
                big_X_d(counter_big_X_d(sample_k(i_l,g_l),type_i(i_l,1)),sample_k(i_l,g_l),type_i(i_l,1),:)=x(:,1)
            end if
            if (sample_k(i_l,g_l)>=1 .and. sample_k(i_l,g_l+1)==clusters+1 ) then
                call TRUNCATED_NORMAL_A_SAMPLE(sum(x(:,1)*beta_d(:,type_i(i_l,1),sample_k(i_l,g_l))),1.0_dp,0.0_dp,d_star)
                big_Y_d(counter_big_X_d(sample_k(i_l,g_l),type_i(i_l,1)),sample_k(i_l,g_l),type_i(i_l,1))=d_star
            elseif (sample_k(i_l,g_l)>=1 .and. sample_k(i_l,g_l+1)<clusters+1  ) then
                call TRUNCATED_NORMAL_B_SAMPLE(sum(x(:,1)*beta_d(:,type_i(i_l,1),sample_k(i_l,g_l))),1.0_dp,0.0_dp,d_star)
                big_Y_d(counter_big_X_d(sample_k(i_l,g_l),type_i(i_l,1)),sample_k(i_l,g_l),type_i(i_l,1))=d_star
            end if
    end do; end do

    beta_d=0.0_dp
    do h_l=1,clusters;do t_l=1,types
        do c_l=1,covariates
            z(c_l,1)=c4_normal_01(  )
        end do
        Sigma=matmul(transpose(big_X_d(1:counter_big_X_d(h_l,t_l),h_l,t_l,:)),big_X_d(1:counter_big_X_d(h_l,t_l),h_l,t_l,:))
        call inverse(Sigma,inv_Sigma,covariates)
        A=inv_Sigma
        call choldc(A,covariates)
        beta_d(:,t_l,h_l)=matmul(inv_Sigma,matmul(transpose(big_X_d(1:counter_big_X_d(h_l,t_l),h_l,t_l,:)),big_Y_d(1:counter_big_X_d(h_l,t_l),h_l,t_l)))+matmul(A,z(:,1))
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
    
