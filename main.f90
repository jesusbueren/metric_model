program main
    use global_var;use nrtype
    implicit none
    integer,dimension(1)::seed=254 
    real(DP),dimension(adls,clusters)::p
    real(DP),dimension(covariates,clusters,clusters+1)::beta
    real(DP),dimension(covariates_habits,habits,types)::gamma
    integer,dimension(indv,1)::y
    integer::i_l
    real(DP)::u
    
    call random_seed(PUT=seed)
    
    !Load original data
    call charge_data()
    
    !Original types are sampled at random
    do i_l=1,indv
        call RANDOM_NUMBER(u)
        if (u<0.5)then
            y(i_l,1)=1
        else
            y(i_l,1)=2
        end if
    end do
    
    !Initial conditions
    call initial_conditions(p,beta,gamma,y)
    open(unit=9,file=path_s//'initial_conditions.txt')
        write(9,'(F20.10)') p,beta,gamma
    close(9)
    
    !Full posterior
    open(unit=9,file=path_s//'initial_conditions.txt')
        read(9,'(F20.10)') p,beta,gamma
    close(9)
    call full_posterior(p,beta,gamma,y)
    
end program