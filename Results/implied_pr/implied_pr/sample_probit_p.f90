subroutine sample_probit_p(y,beta_w) 
    use global_var;use nrtype; use mixtures_vars
    implicit none
    integer,dimension(indv,1),intent(in)::y
    real(DP),dimension(covariates_mix,types,L_educ),intent(inout)::beta_w
    integer::i_l,g_l,y_l,e_l,c_l,age
    real(DP)::u,nw_star
    real(DP),dimension(covariates_mix,1)::x,z
    integer,dimension(types,L_educ)::counter_big_X
    real(DP),dimension(indv*10,types,L_educ,covariates_mix)::big_X
    real(DP),dimension(indv*10,types,L_educ)::big_Y
    real(DP),dimension(covariates_mix,covariates_mix)::Sigma,inv_Sigma,A
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    counter_big_X=0
    do i_l=1,indv;do g_l=first_age(i_l),last_age(i_l)
        age=initial_age+(g_l-1)*2-70
        x(1:4,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0/)     
        if (data_wealth(i_l,g_l)/=-9.0d0 .and. gender(i_l)==1 .and. initial_age+(g_l-1)*2<80) then
            counter_big_X(y(i_l,1),educ(i_l))=counter_big_X(y(i_l,1),educ(i_l))+1
            big_X(counter_big_X(y(i_l,1),educ(i_l)),y(i_l,1),educ(i_l),:)=x(:,1)
            if (data_wealth(i_l,g_l)<=0.0d0 .and.  initial_age+(g_l-1)*2<80) then
                call TRUNCATED_NORMAL_A_SAMPLE(sum(x(:,1)*beta_w(:,y(i_l,1),educ(i_l))),1.0_dp,0.0_dp,nw_star)
                big_Y(counter_big_X(y(i_l,1),educ(i_l)),y(i_l,1),educ(i_l))=nw_star
            else
                call TRUNCATED_NORMAL_B_SAMPLE(sum(x(:,1)*beta_w(:,y(i_l,1),educ(i_l))),1.0_dp,0.0_dp,nw_star)
                big_Y(counter_big_X(y(i_l,1),educ(i_l)),y(i_l,1),educ(i_l))=nw_star
            end if
        end if
    end do; end do
   
    beta_w=0.0_dp
    do y_l=1,types;do e_l=1,L_educ
        do c_l=1,covariates_mix
            z(c_l,1)=c4_normal_01(  )
        end do
        Sigma=matmul(transpose(big_X(1:counter_big_X(y_l,e_l),y_l,e_l,:)),big_X(1:counter_big_X(y_l,e_l),y_l,e_l,:))
        call inverse(Sigma,inv_Sigma,covariates_mix)
        A=inv_Sigma
        call choldc(A,covariates_mix)
        beta_w(:,y_l,e_l)=matmul(inv_Sigma,matmul(transpose(big_X(1:counter_big_X(y_l,e_l),y_l,e_l,:)),big_Y(1:counter_big_X(y_l,e_l),y_l,e_l)))+matmul(A,z(:,1))
    end do;end do
    
end subroutine
    
subroutine sample_probit_p_income(y,beta_i) 
    use global_var;use nrtype; use mixtures_vars_income
    implicit none
    integer,dimension(indv,1),intent(in)::y
    real(DP),dimension(covariates_mix,L_educ),intent(inout)::beta_i
    integer::i_l,g_l,y_l,e_l,c_l,age
    real(DP)::u,i_star
    real(DP),dimension(covariates_mix,1)::x,z
    integer,dimension(L_educ)::counter_big_X
    real(DP),dimension(indv*10,L_educ,covariates_mix)::big_X
    real(DP),dimension(indv*10,L_educ)::big_Y
    real(DP),dimension(covariates_mix,covariates_mix)::Sigma,inv_Sigma,A
    real(DP),dimension(types)::y_d
    real(DP),dimension(cohorts)::cohort_d
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    counter_big_X=0
    do i_l=indv_HRS+1,indv;do g_l=first_age(i_l),last_age(i_l)
        age=initial_age+(g_l-1)*2-70
        y_d=0.0d0
        y_d(y(i_l,1))=1.0d0
        cohort_d=0.0d0
        cohort_d(birth_cohort(i_l))=1.0d0
        x(1:covariates_mix,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(data_shlt(i_l,g_l)-1),cohort_d(4:5)/)     
        if (data_income(i_l,g_l)/=-9.0d0 .and. gender(i_l)==1 .and. initial_age+(g_l-1)*2<62) then
            counter_big_X(educ(i_l))=counter_big_X(educ(i_l))+1
            big_X(counter_big_X(educ(i_l)),educ(i_l),:)=x(:,1)
            if (data_income(i_l,g_l)<=520.0d0*7.25d0 ) then
                call TRUNCATED_NORMAL_A_SAMPLE(sum(x(:,1)*beta_i(:,educ(i_l))),1.0_dp,0.0_dp,i_star)
                big_Y(counter_big_X(educ(i_l)),educ(i_l))=i_star
            else
                call TRUNCATED_NORMAL_B_SAMPLE(sum(x(:,1)*beta_i(:,educ(i_l))),1.0_dp,0.0_dp,i_star)
                big_Y(counter_big_X(educ(i_l)),educ(i_l))=i_star
            end if
        end if
    end do; end do
   
    beta_i=0.0_dp
    do e_l=1,L_educ
        do c_l=1,covariates_mix
            z(c_l,1)=c4_normal_01(  )
        end do
        Sigma=matmul(transpose(big_X(1:counter_big_X(e_l),e_l,:)),big_X(1:counter_big_X(e_l),e_l,:))
        call inverse(Sigma,inv_Sigma,covariates_mix)
        A=inv_Sigma
        call choldc(A,covariates_mix)
        beta_i(:,e_l)=matmul(inv_Sigma,matmul(transpose(big_X(1:counter_big_X(e_l),e_l,:)),big_Y(1:counter_big_X(e_l),e_l)))+matmul(A,z(:,1))
    end do
    

    end subroutine    
    
    subroutine sample_probit_p_income_dynamic(y,beta_i) 
    use global_var;use nrtype; use mixtures_vars_income
    implicit none
    integer,dimension(indv,1),intent(in)::y
    real(DP),dimension(covariates_mix_d,L_educ),intent(inout)::beta_i
    integer::i_l,g_l,y_l,e_l,c_l,age,LF
    real(DP)::u,i_star
    real(DP),dimension(covariates_mix_d,1)::x,z
    integer,dimension(L_educ)::counter_big_X
    real(DP),dimension(indv*10,L_educ,covariates_mix_d)::big_X
    real(DP),dimension(indv*10,L_educ)::big_Y
    real(DP),dimension(covariates_mix_d,covariates_mix_d)::Sigma,inv_Sigma,A
    real(DP),dimension(cohorts)::cohort_d
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    counter_big_X=0
    do i_l=indv_HRS+1,indv;do g_l=first_age(i_l)+1,last_age(i_l)
        if (data_income(i_l,g_l-1)/=-9.0d0 .and.data_income(i_l,g_l)/=-9.0d0 .and. gender(i_l)==1 .and. initial_age+(g_l-1)*2<62) then
            age=initial_age+(g_l-1)*2-70
            if (data_income(i_l,g_l-1)<=520.0d0*7.25d0) then
                LF=0
            else
                LF=1
            end if
            cohort_d=0.0d0
            cohort_d(birth_cohort(i_l))=1.0d0
            x(1:covariates_mix_d,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(data_shlt(i_l,g_l)-1),dble(LF),cohort_d(4:5)/) 
            counter_big_X(educ(i_l))=counter_big_X(educ(i_l))+1
            big_X(counter_big_X(educ(i_l)),educ(i_l),:)=x(:,1)
            if (data_income(i_l,g_l)<=520.0d0*7.25d0 ) then
                call TRUNCATED_NORMAL_A_SAMPLE(sum(x(:,1)*beta_i(:,educ(i_l))),1.0_dp,0.0_dp,i_star)
                big_Y(counter_big_X(educ(i_l)),educ(i_l))=i_star
            else
                call TRUNCATED_NORMAL_B_SAMPLE(sum(x(:,1)*beta_i(:,educ(i_l))),1.0_dp,0.0_dp,i_star)
                big_Y(counter_big_X(educ(i_l)),educ(i_l))=i_star
            end if
        end if
    end do; end do
   
    beta_i=0.0_dp
    do e_l=1,L_educ
        do c_l=1,covariates_mix_d
            z(c_l,1)=c4_normal_01(  )
        end do
        Sigma=matmul(transpose(big_X(1:counter_big_X(e_l),e_l,:)),big_X(1:counter_big_X(e_l),e_l,:))
        call inverse(Sigma,inv_Sigma,covariates_mix_d)
        A=inv_Sigma
        call choldc(A,covariates_mix_d)
        beta_i(:,e_l)=matmul(inv_Sigma,matmul(transpose(big_X(1:counter_big_X(e_l),e_l,:)),big_Y(1:counter_big_X(e_l),e_l)))+matmul(A,z(:,1))
    end do
    

end subroutine    
    

    
    
double precision function c4_normal_01 (  )
!------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  double precision, parameter :: r4_pi=3.14159265358979323846264338327950288419716939937510
  double precision:: v1
  double precision:: v2
  double precision:: x_c
  double precision:: x_r 
  call random_number(v1)
  call random_number(v2)
  x_r =sqrt(-2.0d0*log(v1))*cos(2.0d0*r4_pi*v2)
  x_c =sqrt(-2.0d0*log(v1))*sin(2.0d0*r4_pi*v2)
  c4_normal_01=x_r
  return
end
    

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

SUBROUTINE choldc(a,n)
use nrtype
    IMPLICIT NONE
    integer,intent(in)::n
    real(DP), DIMENSION(n,n), INTENT(INOUT) :: a
    real(DP), DIMENSION(n) :: p
    INTEGER :: i,j
    real(DP) :: summ
    do i=1,n
	    summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
	    p(i)=sqrt(summ)
	    a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
    end do
    
    do i=1,n
        do j=1,n
            if (i==j) then
                a(i,i)=p(i)
            elseif (i<j) then
                a(i,j)=0.0d0
            end if
        end do
    end do
    
END SUBROUTINE choldc
    