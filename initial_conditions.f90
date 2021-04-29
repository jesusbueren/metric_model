module global_ini
    use global_var
    implicit none
    integer,dimension(indv,generations)::sample_k_ini
    integer,dimension(indv,1)::y_ini
end module
    
subroutine initial_conditions(beta_h,beta_d,gamma,y)
    use global_var; use nrtype; use global_ini
    implicit none
    real(DP),dimension(covariates,types,clusters,clusters),intent(out)::beta_h
    real(DP),dimension(covariates,types,clusters),intent(out)::beta_d
    real(DP),dimension(covariates_habits,habits,types),intent(out)::gamma
    integer,dimension(indv,1),intent(in)::y
    
    !Assign individuals to clusters
    call cluster_assign(sample_k_ini)
    !Initial conditions for y is set by a random number generator
    y_ini=y
    !Posterior distribution of transition parameters given the sampled health states & health types
    print*,'initial_conditions_survival'
    call initial_conditions_sur(beta_d) 
    print*,'initial_conditions_tr'
    call initial_conditions_tr(beta_h) 

    !Posterior distribution of the habits given the sampled health states
    print*,'initial_conditions_habits'
    call initial_conditions_habits(gamma)
      
end subroutine
    
subroutine cluster_assign(sample_k)
    use global_var; use nrtype
    implicit none
    integer,dimension(indv,generations),intent(out)::sample_k
    real(DP)::u
    integer::i_l,g_l,ind
    
    sample_k=data_shlt
    
end subroutine
    
subroutine initial_conditions_tr(beta_h)
use global_var; use nrtype; use global_ini
implicit none
real(DP),dimension(covariates,types,clusters,clusters),intent(out)::beta_h
real(DP),dimension(covariates*clusters**2,1)::c_tr
integer::it

beta_h=0.0_dp
do it=1,100
    print*,it
    call sample_beta_h(beta_h,y_ini,sample_k_ini) 
end do

end subroutine
    
subroutine initial_conditions_sur(beta_d)
use global_var; use nrtype; use global_ini
implicit none
real(DP),dimension(covariates,types,clusters),intent(out)::beta_d
real(DP),dimension(covariates*clusters,1)::c_tr
integer::it

beta_d=0.0_dp
do it=1,100
    print*,it
    !Sample beta given h_star
    call sample_beta_d(beta_d,y_ini,sample_k_ini)
end do

end subroutine    

subroutine initial_conditions_habits(gamma)
use global_var; use nrtype; use global_ini
implicit none
real(DP),dimension(covariates_habits,habits,types),intent(out)::gamma
integer::ind,it
real(DP),dimension(covariates_habits*habits*types,1)::c_ga
real(DP)::factor_A=1.0_dp

gamma=0.0_dp
do it=0,100
    print*,it
    !Sample y* given gamma
    call sample_gamma_y(gamma,y_ini,sample_k_ini)    
end do

end subroutine
    
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
    
 

