subroutine full_posterior(beta_h,beta_d,gamma,y,delta)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates_mixture,L_gender,L_educ,types),intent(inout)::delta
    real(DP),dimension(covariates,types,clusters,L_gender,L_educ),intent(inout)::beta_h
    real(DP),dimension(covariates_habits,habits,types),intent(inout)::gamma
    integer,dimension(indv,1),intent(inout)::y
    real(DP),dimension(covariates,types,clusters,L_gender,L_educ),intent(inout)::beta_d
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H,H_g 
 
    real(DP),dimension(types,L_gender,L_educ)::fraction_t
    integer::g_l,it3,it2,chg_mu,it,v_l,ind,ge_l,e_l,i_l,burn
    real(DP)::prior_H
    integer,dimension(indv,generations)::sample_k
    real(DP)::u
    real(DP),dimension(indv,habits,generations)::y_star
    real(DP),dimension(types,L_gender,L_educ,clusters+1)::LE
    real(DP),dimension(generations,clusters,L_gender,L_educ,types)::weights,joint_yh
    real(DP),dimension(clusters,L_gender,L_educ)::share_h
    !Timer
    integer::calc
    real::calctime
    interface
        real function tock(t)
            integer, intent(in) :: t
        end function tock
    end interface
    character::end_k

    character::continue_program
    
    !Compute share of indv in good and bad health in the initial period across education and gender
    call fraction_h_e_g(share_h)
    
    sample_k=data_shlt
    H=1/dble(clusters+1)
    burn=1
    
    
    it2=-1
    it3=-1
    chg_mu=0

    !call tick(calc)
    do it=1,2000+burn
        print*,it
        it2=it2+1 
        it3=it3+1 
        !Sample h* given beta 
        !print*,'2'
        call sample_beta_h(beta_h,y,sample_k)
        !Sample d* given beta 
        !print*,'3'
        call sample_beta_d(beta_d,y,sample_k)
        !Sample y* given gamma
        !print*,'4'
        call sample_gamma_y(gamma,y,sample_k) 
        !Compute life-expectancy
        !print*,'6'
        call transitions(beta_h,beta_d,H,LE,joint_yh)
        !sample weights
        call sample_delta(delta,y,sample_k,share_h)
        !compute weights
        call compute_weights(delta,H,share_h,weights,joint_yh) 
        !print*,'1'
        call sample_y(gamma,y,sample_k,H,weights)
        !Save results
        if (it>burn) then
            call save_results(beta_h,beta_d,gamma,LE,sum(joint_yh,2),it-burn)
        end if
    end do
    
end subroutine
    
subroutine tick(t)
    integer, intent(OUT) :: t
    call system_clock(t)
end subroutine tick

! returns time in seconds from now to time described by t
real function tock(t)
    integer, intent(in) :: t
    integer :: now, clock_rate
    call system_clock(now,clock_rate)
    tock = real(now - t)/real(clock_rate)
end function tock