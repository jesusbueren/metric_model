subroutine full_posterior(beta_h,beta_d,gamma,gamma_med,y,delta)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates_mixture*(types-1),L_gender,L_educ),intent(inout)::delta
    real(DP),dimension(covariates,clusters,L_gender,L_educ),intent(inout)::beta_h
    real(DP),dimension(covariates_habits,habits_nomed,types),intent(inout)::gamma
    real(DP),dimension(covariates_habits_med,habits_med,types),intent(inout)::gamma_med
    integer,dimension(indv,1),intent(inout)::y
    real(DP),dimension(covariates,clusters,L_gender,L_educ),intent(inout)::beta_d
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H,H_g 
    integer::it,burn,it2,h_l
    integer,dimension(indv,generations)::sample_k
    real(DP)::u
    real(DP),dimension(indv,habits,generations)::y_star
    real(DP),dimension(types,L_gender,L_educ,clusters+1)::LE
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts)::weights,joint_yh,aux
    real(DP),dimension(clusters,L_gender,L_educ)::share_h
    real(DP),dimension(indv,types)::type_pr,type_pr_av
    real(DP),dimension(covariates,covariates,clusters,L_gender,L_educ)::sigma_h,sigma_d
    real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_h_mean,beta_d_mean
    real(DP),dimension(covariates_mixture*(types-1),L_gender,L_educ)::mean_delta
    integer,dimension(L_gender,L_educ)::acc_delta
    real(DP),dimension(covariates_mixture*(types-1),covariates_mixture*(types-1),L_gender,L_educ)::cov_delta
    real(DP),dimension(L_gender,L_educ)::shrinkage
    real(DP),dimension(types,clusters,L_gender,L_educ)::shrinkage_h,shrinkage_d
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
    
    joint_yh=1.0d0/dble(clusters*types)
    weights=1.0d0/dble(types)
    mean_delta=0.0d0
    cov_delta=0.0d0
    acc_delta=0
    shrinkage=1.0d0
    shrinkage_h=1.0d0
    shrinkage_d=1.0d0
    !Burn iterations (avoid saving results before iteration)
    burn=0    

    type_pr_av=0.0d0
    it2=0
    beta_h_mean=-9.0d0
    sigma_h=-9.0d0
    beta_d_mean=-9.0d0
    sigma_d=-9.0d0
    do it=1,500000+burn
        print*,it

        !Sample health transitions parameters: gh/bh
        if (it<=2000) then
            call sample_beta_h(beta_h,y,sample_k)
        else
            call sample_beta_h_MH(beta_h,beta_d,share_h,H,y,sample_k,weights,joint_yh,beta_h_mean,sigma_h,it-2000,shrinkage_h)
        end if
        
        !Sample health transitions parameters: survival
        if (it<=2000) then
            call sample_beta_d(beta_d,y,sample_k)
        else
            call sample_beta_d_MH(beta_d,beta_h,share_h,H,y,sample_k,weights,joint_yh,beta_d_mean,sigma_d,it-2000,shrinkage_h)
        end if
        

        !Sample health behavior parameters
        !a) no med
        call sample_gamma_y(gamma,y,sample_k) 
        !b) med
        call sample_gamma_y_med(gamma_med,y,sample_k) 
        
        !Compute transitions and life expectancies        
        compute_LE=1
        call transitions(beta_h,beta_d,H,LE,joint_yh) 
            

        !Sample pr of type at initial age
        call sample_delta(delta,mean_delta,cov_delta,H,share_h,y,sample_k,weights,joint_yh,it,acc_delta,shrinkage)
        !weights(1,2,1,1,1,1) 
        !sample type
        
        

        call sample_y(gamma,gamma_med,y,sample_k,H,weights,type_pr)




        
        if (it>burn) then
            if (mod(it,10) == 0) then
                do h_l=1,clusters
                    aux(:,h_l,:,:,:,:)=sum(joint_yh,2) 
                end do
                call save_results(beta_h,beta_d,gamma,gamma_med,delta,LE,sum(joint_yh,2),joint_yh/aux,H,it2) !H(2,3,:,2,1,3)
                it2=it2+1
            end if
            type_pr_av=dble(it-burn-1)/dble(it-burn)*type_pr_av+1.0d0/dble(it-burn)*type_pr
        end if
    end do
    
    open(unit=9,file=path_s//'implied_probilities_'//types_s//'.txt')
        write(9,'(F20.8)') type_pr_av
    close(9)
    
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