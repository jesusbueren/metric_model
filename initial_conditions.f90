module global_ini
    use global_var
    implicit none
    integer,dimension(indv,generations)::sample_k_ini
    integer,dimension(indv,1)::y_ini
end module
    
subroutine initial_conditions(beta,gamma,y)
    use global_var; use nrtype; use global_ini
    implicit none
    real(DP),dimension(covariates,types,clusters,clusters+1),intent(out)::beta
    real(DP),dimension(covariates_habits,habits,types),intent(out)::gamma
    integer,dimension(indv,1),intent(in)::y
    
    !Assign individuals to clusters
    call cluster_assign(sample_k_ini)
    !Initial conditions for y is set by a random number generator
    y_ini=y
    !Posterior distribution of transition parameters given the sampled health states & health types
    print*,'initial_conditions_tr'
    call initial_conditions_tr(beta)  
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
    
subroutine initial_conditions_tr(beta)
use global_var; use nrtype; use global_ini
implicit none
real(DP),dimension(covariates,types,clusters,clusters+1),intent(out)::beta
real(DP),dimension(covariates*clusters**2,1)::c_tr
real(DP),dimension(clusters,generations,types,L_gender)::dist_init
integer::it

beta=0.0_dp
do it=1,100
    print*,it
    !Sample h* given beta
    call sample_h_star(beta,y_ini,sample_k_ini) 
    !Sample beta given h_star
    call sample_beta(beta)
end do

end subroutine

subroutine initial_conditions_habits(gamma)
use global_var; use nrtype; use global_ini
implicit none
real(DP),dimension(covariates_habits,habits,types),intent(out)::gamma
integer::ind,it
real(DP),dimension(covariates_habits*habits*types,1)::c_ga
real(DP)::factor_A=1.0_dp
real(DP),dimension(indv,habits,generations)::y_star

gamma=0.0_dp
do it=0,100
    print*,it
    !Sample y* given gamma
    call sample_y_star(gamma,y_ini,sample_k_ini,y_star)
    !Sample gamma given y*
    call sample_gamma(y_ini,y_star,gamma)
    
    c_ga=reshape(gamma,(/covariates_habits*habits*types,1/))
    if (it==1) then
        open(unit=9,file=path_s//'c_habits.txt')
            write(9,'(F20.8)') c_ga
        close(9)
    else
        open(unit=9,file=path_s//'c_habits.txt',access='append')
            write(9,'(F20.8)') c_ga
        close(9)
    end if
    
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
    
 

