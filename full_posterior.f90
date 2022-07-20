subroutine full_posterior(beta_h,beta_d,gamma,y,delta)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(clusters,L_gender,L_educ,types),intent(inout)::delta
    real(DP),dimension(covariates,types,clusters,L_gender,L_educ),intent(inout)::beta_h
    real(DP),dimension(covariates_habits,habits,types),intent(inout)::gamma
    integer,dimension(indv,1),intent(inout)::y
    real(DP),dimension(covariates,types,clusters,L_gender,L_educ),intent(inout)::beta_d
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H,H_g 
    integer::it,burn
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
    
    sample_k=data_shlt
    
    !Compute share of indv in good and bad health in the initial period across education and gender
    call fraction_h_e_g(sample_k,share_h)
    
    H=1/dble(clusters+1)
    burn=0
    
    

    beta_h=0.0d0
    beta_d=0.0d0
    gamma=0.0d0
    delta=1.0d0/dble(types)
    

    !call tick(calc)
    do it=1,4000+burn
        print*,it
        !Sample health transitions parameters
        call sample_beta_h(beta_h,y,sample_k)
        !Sample survival parameters
        call sample_beta_d(beta_d,y,sample_k)
        !Sample health behavior parameters
        call sample_gamma_y(gamma,y,sample_k) 
        !Compute transitions and life expectancies
        call transitions(beta_h,beta_d,H,LE,joint_yh) 
        !Sample pr of type at initial age
        call sample_delta(delta,H,share_h,y,sample_k,weights,joint_yh)
        !sample type
        call sample_y(gamma,y,sample_k,H,weights)

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