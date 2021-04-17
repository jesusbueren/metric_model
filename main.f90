program main
    use global_var;use nrtype
    implicit none
    integer,dimension(1)::seed=254 
    real(DP),dimension(covariates,types,clusters,clusters+1)::beta
    real(DP),dimension(covariates_habits,habits,types)::gamma
    integer,dimension(indv,1)::y
    integer::i_l,ind
    real(DP)::u
    
    call random_seed(PUT=seed)
    
    !Load original data
    call charge_data()
    
    !Original types are sampled at random
    y=-1
    do i_l=1,indv
        call RANDOM_NUMBER(u)
        ind=1
        do while (y(i_l,1)==-1)
            if (u<1.0d0/dble(types)*dble(ind))then
                y(i_l,1)=ind
            else
                ind=ind+1
            end if
        end do
    end do
    
    !Initial conditions
    call initial_conditions(beta,gamma,y)
    open(unit=9,file=path_s//'initial_conditions.txt')
        write(9,'(F20.10)') beta,gamma
    close(9)
    
    !Full posterior
    open(unit=9,file=path_s//'initial_conditions.txt')
        read(9,'(F20.10)') beta,gamma
    close(9)
    call full_posterior(beta,gamma,y)
    
end program