program main
    use global_var;use nrtype
    implicit none
    integer,dimension(1)::seed=789 
    real(DP),dimension(covariates,types,clusters,clusters)::beta_h
    real(DP),dimension(covariates,types,clusters)::beta_d
    real(DP),dimension(covariates_habits,habits,types)::gamma
    double precision,dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H
    real(DP),dimension(clusters,types,L_gender,L_educ)::init_cond 
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
    call load_high_density(beta_h,beta_d,gamma)
    

    call sample_y(gamma,y,fraction_t)

    
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
                elseif (data_habits(i_l,h_l,g_l)==0) then
                    pr(e_l)=pr(e_l)*(1.0d0-alphas(h_l,g_l,e_l))
                end if
            end do
        end do; end do
        pr=pr/sum(pr,1)
        type_pr(i_l,:)=pr
    end do
    
    open(unit=9,file=path_s//'pr_type_hrs.txt',status='replace')
    do i_l=1,indv_HRS; do g_l=1,generations;
        write(9,'(I5,I3,<types>F6.3)'), i_l,g_l,type_pr(i_l,:)
    end do;end do
    close(9)
    
    open(unit=9,file=path_s//'pr_type_psid.txt',status='replace')
    do i_l=1,indv_PSID; do g_l=1,generations;
        write(9,'(I5,I3,<types>F6.3)'), i_l+indv_HRS,g_l,type_pr(i_l+indv_HRS,:)
    end do;end do
    close(9)
    
end