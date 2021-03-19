program main
    use global_var;use nrtype
    implicit none
    integer,dimension(1)::seed=789 
    double precision,dimension(adls,clusters)::p
    real(DP),dimension(covariates,clusters,clusters+1)::beta
    real(DP),dimension(covariates_habits,habits,types)::gamma
    double precision,dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H
    real(DP),dimension(clusters,types,L_gender,L_educ)::init_cond 
    double precision,dimension(indv,clusters,generations)::likelihood
    double precision,dimension(indv,clusters,generations)::filtered_states
    double precision,dimension(indv,clusters+1,generations)::smoothed_states
    integer,dimension(indv,generations)::sample_k
    real(DP),dimension(types,L_gender,L_educ,clusters+1)::LE
    character::end_of_program
    integer::ns,burn,g_l,i_l,it,e_l,h_l,age
    double precision::counter,chg
    integer,dimension(indv,1)::y
    real(DP),dimension(types,L_gender,L_educ)::fraction_t
    double precision,dimension(types)::pr
    double precision,dimension(indv,types)::type_pr
    real(DP),dimension(covariates_habits,1)::x
    real(DP),dimension(habits,generations,types)::alphas
    
    
    call random_seed(PUT=seed) 
    call charge_data()
    
    !Load high density point from estimated posterior distribution 
    ! and take it as true value
    call load_high_density(p,beta,gamma)
    
    print*,'Compute smoothed probability/full info for econometritian'
    !Initial guess for the smoothed states
    smoothed_states=1.0d0/dble(clusters)
    call likelihood_all(p,likelihood)
    do it=1,100 !I run it 100 so that initial conditions for filter don't matter
        print*,it
        call transitions(beta,init_cond,H,LE) 
        call sample_y(gamma,y,fraction_t)
        call filtration(H,p,smoothed_states,likelihood,y,filtered_states)
        call smoothing(filtered_states,H,y,sample_k,smoothed_states,init_cond) 
    end do
    
     do e_l=1,types
        do h_l=1,habits; do g_l=1,generations
            age=initial_age+(g_l-1)*2-70
            x(:,1)=(/1.0_dp,dble(age),dble(age**2.0_dp-1.0_dp)/)
            alphas(h_l,g_l,e_l)=1.0_dp-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*gamma(:,h_l,e_l))/sqrt(2.0_dp)))
        end do; end do
    end do
    
    do i_l=1,indv;
        pr=1.0_dp
        do g_l=first_age(i_l),last_age(i_l);do h_l=1,habits
            do e_l=1,types
                if (data_habits(i_l,h_l,g_l)==1) then
                    pr(e_l)=pr(e_l)*alphas(h_l,g_l,e_l)
                end if
            end do
        end do; end do
        pr=pr/sum(pr,1)
        type_pr(i_l,:)=pr
    end do

    open(unit=9,file=path_s//'smoothed_pr_.txt',status='replace')
    do i_l=1,indv; do g_l=1,generations;
        write(9,'(I5,I3,F6.3,F6.3,I3,F6.3,F6.3,F6.3)'), i_l,g_l,smoothed_states(i_l,1,g_l),smoothed_states(i_l,2,g_l),maxloc(smoothed_states(i_l,:,g_l)),type_pr(i_l,1),type_pr(i_l,2)
    end do;end do
    close(9)
    
end