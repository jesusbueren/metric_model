  
subroutine sample_beta_d(beta_d,type_i,sample_k)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates,clusters,L_gender,L_educ),intent(inout)::beta_d
    integer,dimension(indv,1),intent(in)::type_i
    integer,dimension(indv,generations),intent(in)::sample_k
    integer::h_l,c_l,g_l,ge_l,age,ge_d,it,i_l,health_d,d_l,e_l
    real(DP)::d_star
    real(DP)::gender_d
    real(DP),dimension(covariates,1)::z
    real(DP),dimension(covariates,covariates)::Sigma,inv_Sigma,A
    real(DP),dimension(indv*g_max,clusters,L_gender,L_educ,covariates)::big_X_d
    real(DP),dimension(indv*g_max,clusters,L_gender,L_educ)::big_Y_d
    integer,dimension(clusters,L_gender,L_educ)::counter_big_X_d
    real(DP),dimension(types-1)::dummy_type,dummy_type_x_age
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    real(DP),dimension(covariates,1)::x
    character::pause_k

     counter_big_X_d=0
    do i_l=1,indv_HRS;do g_l=first_age(i_l),last_age(i_l)-1
            x=-9.0d0
            age=initial_age+(g_l-1)*2-70
            dummy_type=0.0d0
            dummy_type_x_age=0.0d0
            if (type_i(i_l,1)>1)then
                dummy_type(type_i(i_l,1)-1)=1.0d0
                dummy_type_x_age(type_i(i_l,1)-1)=dble(age)
            end if
            x(:,1)=[(/1.0_dp,dble(age)/),dummy_type,dummy_type_x_age]!,dble(age)**2.0d0
            if (sample_k(i_l,g_l)>=1 .and. sample_k(i_l,g_l+1)>=1 .and. race(i_l)==1) then
                counter_big_X_d(sample_k(i_l,g_l),gender(i_l),educ(i_l))=counter_big_X_d(sample_k(i_l,g_l),gender(i_l),educ(i_l))+1
                big_X_d(counter_big_X_d(sample_k(i_l,g_l),gender(i_l),educ(i_l)),sample_k(i_l,g_l),gender(i_l),educ(i_l),:)=x(:,1)
                if (sample_k(i_l,g_l+1)==clusters+1 ) then
                    call TRUNCATED_NORMAL_A_SAMPLE(sum(x(:,1)*beta_d(:,sample_k(i_l,g_l),gender(i_l),educ(i_l))),1.0_dp,0.0_dp,d_star)
                elseif (sample_k(i_l,g_l+1)<clusters+1  ) then
                    call TRUNCATED_NORMAL_B_SAMPLE(sum(x(:,1)*beta_d(:,sample_k(i_l,g_l),gender(i_l),educ(i_l))),1.0_dp,0.0_dp,d_star)
                end if
                if (d_star>100 .or. d_star<-100) then
                    print*,d_star
                    print*,x(:,1)
                    print*,beta_d(:,sample_k(i_l,g_l),gender(i_l),educ(i_l))
                    read*,pause_k
                end if
                big_Y_d(counter_big_X_d(sample_k(i_l,g_l),gender(i_l),educ(i_l)),sample_k(i_l,g_l),gender(i_l),educ(i_l))=d_star
            end if
            
    end do; end do

    do h_l=1,clusters;do e_l=1,L_educ;do ge_l=1,L_gender
        do c_l=1,covariates
            z(c_l,1)=c4_normal_01(  )
        end do
        if (counter_big_X_d(h_l,ge_l,e_l)>1) then
            beta_d(:,h_l,ge_l,e_l)=0.0d0
            Sigma=matmul(transpose(big_X_d(1:counter_big_X_d(h_l,ge_l,e_l),h_l,ge_l,e_l,:)),big_X_d(1:counter_big_X_d(h_l,ge_l,e_l),h_l,ge_l,e_l,:))
            call inverse(Sigma,inv_Sigma,covariates)
            A=inv_Sigma
            call choldc(A,covariates)
            beta_d(:,h_l,ge_l,e_l)=matmul(inv_Sigma,matmul(transpose(big_X_d(1:counter_big_X_d(h_l,ge_l,e_l),h_l,ge_l,e_l,:)),big_Y_d(1:counter_big_X_d(h_l,ge_l,e_l),h_l,ge_l,e_l)))+matmul(A,z(:,1))
            if (isnan(sum(beta_d(:,h_l,ge_l,e_l)))) then
                print*,'pb beta_d'
            end if
        else
            print*,'strange sampled beta_d'
        end if
    end do;end do; end do
        
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
    
