program main
    use global_var;use nrtype
    implicit none
    integer,dimension(1)::seed=456
    real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_h
    real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_d
    real(DP),dimension(covariates_habits,habits_nomed,types)::gamma
    real(DP),dimension(covariates_habits_med,habits_med,types)::gamma_med
    real(DP),dimension(covariates_mixture*(types-1),L_gender,L_educ)::delta
    integer,dimension(indv,1)::y
    integer::i_l
    real(DP)::u
    
    call random_seed(PUT=seed)
    
    !Load original data
    call charge_data()
    
    !Original types are sampled with higher probability of smokers are type 3, if no exercise type 2 else type 1
    ! This is irrelevant for the estimation results but it ensures that the type 1 is always the protective, Type 2 detrimental and type 3 harmful
    y=-1
    do i_l=1,indv
        
        if (data_habits(i_l,3,first_age(i_l))==1) then !smoking
                y(i_l,1)=types
        end if
        !if (data_habits(i_l,6,first_age(i_l))==1) then !Exercise
        !    call random_number(u)
        !    if (u<0.9d0)then
        !        y(i_l,1)=2
        !    else
        !        y(i_l,1)=1
        !    end if
        !end if
        !if (y(i_l,1)==-1 ) then
        !    call random_number(u)
        !    if (u<0.5d0)then
        !        y(i_l,1)=2
        !    else
        !        y(i_l,1)=1
        !    end if
        !end if
        !if (types>2 .and. y(i_l,1)==-1) then
        !    if (data_habits(i_l,6,first_age(i_l))==0) then !Exercise
        !        if (u<0.9d0)then
        !            y(i_l,1)=2
        !        end if
        !    else
        !        y(i_l,1)=types-1
        !    end if
        !end if
        !if (types==4.and. y(i_l,1)==-1) then
        !    if (data_habits(i_l,2,first_age(i_l))==1) then !Exercise
        !        if (u<0.9d0)then
        !            y(i_l,1)=3
        !        end if
        !    else
        !        y(i_l,1)=2
        !    end if
        !end if
        !
        !if (data_habits(i_l,1,first_age(i_l))==1 .and. y(i_l,1)==-1) then  !cancer test
        !    if (u<0.8d0)then
        !        y(i_l,1)=1
        !    end if
        !end if
        if (y(i_l,1)==-1 .and. types==2) then
            y(i_l,1)=1
        end if
        !if (y(i_l,1)==-1 .and. types==3) then
        !    y(i_l,1)=2
        !end if
        !if (y(i_l,1)==-1 .and. types==4) then
        !    y(i_l,1)=3
        !end if
    end do

    
    !Full posterior
    
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