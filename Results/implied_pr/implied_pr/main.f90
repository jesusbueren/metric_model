program main
    use global_var;use nrtype
    implicit none
    integer,dimension(1)::seed=789 
    real(DP),dimension(covariates,types,clusters,clusters)::beta_h
    real(DP),dimension(covariates,types,clusters)::beta_d
    real(DP),dimension(covariates_habits,habits,types)::gamma
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H
    real(DP),dimension(clusters,types,L_gender,L_educ)::init_cond 
    integer,dimension(indv,generations)::sample_k
    real(DP),dimension(types,L_gender,L_educ,clusters+1)::LE
    character::end_of_program
    integer::ns,burn,g_l,i_l,it,e_l,h_l,age,c_l
    real(DP)::counter,chg,health_d
    integer,dimension(indv,1)::y
    real(DP),dimension(types,L_gender,L_educ)::fraction_t
    real(DP),dimension(types)::pr
    real(DP),dimension(indv,types)::type_pr
    real(DP),dimension(covariates_habits,1)::x
    real(DP),dimension(habits,generations,types,clusters)::alphas
    
    
    call random_seed(PUT=seed) 
    call charge_data()
    
    
    sample_k=data_shlt
    
    open(unit=9,file=path_s//'implied_probilities.txt')
        read(9,'(F20.8)') type_pr
    close(9)
    
    open(unit=9,file=path_s//'pr_type_hrs2.txt')
    do i_l=1,indv_HRS; do g_l=1,generations;
        write(9,'(I5,I3,<types>F6.3)'), i_l,g_l,type_pr(i_l,:)
    end do;end do
    close(9)
    
    open(unit=9,file=path_s//'pr_type_psid.txt')
    do i_l=1,indv_PSID; do g_l=1,generations;
        write(9,'(I5,I3,<types>F6.3)'), i_l+indv_HRS,g_l,type_pr(i_l+indv_HRS,:)
    end do;end do
    close(9)
    
    !call estimate_mixture_wealth(type_pr)
    
    call estimate_mixture_income(type_pr)
    
end program