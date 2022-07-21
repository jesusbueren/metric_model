program main
    use global_var;use nrtype
    implicit none
    integer,dimension(1)::seed=456
    real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_h
    real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_d
    real(DP),dimension(covariates_habits,habits,types)::gamma
    real(DP),dimension(covariates_mixture,L_gender,L_educ,types)::delta
    integer,dimension(indv,1)::y
    integer::i_l,ind
    real(DP)::u
    
    call random_seed(PUT=seed)
    
    !Load original data
    call charge_data()
    
    !Original types are sampled with higher probability of smokers are type 3, if high bmi type 2 else type 1
    y=-1
    do i_l=1,indv
        call random_number(u)
        if (data_habits(i_l,3,first_age(i_l))==1) then
            if (u<0.8d0)then
                y(i_l,1)=3
            end if
        end if
        if (data_habits(i_l,6,first_age(i_l))==1 .and. y(i_l,1)==-1) then
            if (u<0.8d0)then
                y(i_l,1)=2
            end if
        end if
        if (y(i_l,1)==-1) then
            y(i_l,1)=1
        end if
    end do
    

    
    
    !call simulate_data()
    
    !Initial conditions
    !call initial_conditions(beta_h,beta_d,gamma,y,delta)
    !open(unit=9,file=path_s//'initial_conditions.txt')
    !    write(9,'(F20.10)') beta_h,gamma,beta_d,delta
    !close(9)
    
    !Full posterior
    open(unit=9,file=path_s//'initial_conditions.txt')
        read(9,'(F20.10)') beta_h,gamma,beta_d,delta
    close(9)
    call full_posterior(beta_h,beta_d,gamma,y,delta)
    
end program