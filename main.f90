program main
    use global_var;use nrtype
    implicit none
    integer,dimension(1)::seed=123
    real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_h
    real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_d
    real(DP),dimension(covariates_habits,habits_nomed,types)::gamma
    real(DP),dimension(covariates_habits_med,habits_med,types)::gamma_med
    real(DP),dimension(covariates_mixture*(types-1),L_gender,L_educ)::delta
    integer,dimension(indv,1)::y
    integer::i_l,g_l,ind
    real(DP)::u
    real(DP),dimension(indv,types-1)::type_pr_1
    
    call random_seed(PUT=seed)
    
    !Load original data
    call charge_data()
    
    !Initial types
    y=-1
    !do i_l=1,indv 
    !    if (sample_selection(i_l)) then
    !        call RANDOM_NUMBER(u)
    !        ind=1
    !        do while (y(i_l,1)==-1)
    !            if (u<1.0d0/dble(types)*ind) then
    !                y(i_l,1)=ind
    !            else
    !                ind=ind+1
    !            end if
    !        end do
    !    end if
    !end do
        
    do i_l=1,indv; ;do g_l=first_age(i_l),last_age(i_l)        
        if (data_habits(i_l,3,g_l)==1) then !smoking 
                y(i_l,1)=types
                exit
        end if
        !if (data_habits(i_l,2,g_l)==0 .and. y(i_l,1)==-1 .and. types>2) then !no cancer test
        !    y(i_l,1)=types-1
        !end if
        if (data_habits(i_l,5,g_l)==1 .and. y(i_l,1)==-1) then !flu shot
            y(i_l,1)=1
        end if
    end do;
        if (y(i_l,1)==-1) then
            y(i_l,1)=2
            !if (types>2) then
            !    y(i_l,1)=2
            !end if
        end if
    end do      

    
    !Initial guess
    beta_h=0.0d0
    beta_d=0.0d0
    gamma=0.0d0
    delta=0.0d0
    gamma_med=0.0d0
    call full_posterior(beta_h,beta_d,gamma,gamma_med,y,delta)
    
   !call simulate_model()
    
    pause
    
end program