subroutine full_posterior(p,beta,gamma,y)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(adls,clusters),intent(inout)::p
    real(DP),dimension(covariates,clusters,clusters+1),intent(inout)::beta
    real(DP),dimension(covariates_habits,habits,types),intent(inout)::gamma
    integer,dimension(indv,1),intent(inout)::y
    real(DP),dimension(adls,clusters)::p_g
    real(DP),dimension(adls*clusters,1)::c_p,c_p_g,mean_p,new_mean_p
    real(DP),dimension(adls*clusters,adls*clusters)::C_rnd_W_p,C_A_p
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H,H_g 
    real(DP),dimension(clusters,types,L_gender,L_educ)::init_cond 
    integer::c_l,c_l2,it3,it2,chg_mu,it,burn=100,v_l,ind
    real(DP),dimension(indv,clusters,generations)::likelihood
    real(DP)::factor_p,posterior_mu,prior_H,s_d2=2.4**2/(adls*clusters),posterior_mu_g
    real(DP),dimension(indv,clusters,generations)::filtered_states
    real(DP),dimension(indv,clusters+1,generations)::smoothed_states
    integer,dimension(indv,generations)::sample_k
    real(DP)::u,factor_mu
    real(DP),dimension(indv,habits,generations)::y_star
    real(DP),dimension(types,L_gender,L_educ,clusters+1)::LE
    real(DP),dimension(types,L_gender,L_educ)::fraction_t
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
    
    factor_mu=1.0_dp
    factor_p=1.0_dp/10000.0_dp
    C_rnd_W_p=0.0d0
    do c_l=1,adls*clusters; do c_l2=1,adls*clusters
        if (c_l==c_l2) then
            C_rnd_W_p(c_l,c_l2)=factor_p
        end if
    end do;end do
    
    !Initial guess for the smoothed states
    smoothed_states=0.0d0
    smoothed_states(:,1:clusters,:)=1.0d0/dble(clusters)
    
    call likelihood_all(p,likelihood)
    
    posterior_mu=-1.0/0d0
    it2=-1
    it3=-1
    chg_mu=0
    !call tick(calc)
    do it=0,100000000
        print*,it
        it2=it2+1 
        it3=it3+1 
        !sample health status
        call transitions(beta,init_cond,H,LE)  
        call filtration(H,p,smoothed_states,likelihood,y,filtered_states)
        call smoothing(filtered_states,H,y,sample_k,smoothed_states,init_cond)
        !Sample health behavior type
        call sample_y(gamma,y,fraction_t)
        !Sample p by cluster
        c_p=reshape(p,(/adls*clusters,1/))
        call compute_mean(mean_p,c_p,it,clusters*adls,new_mean_p)
        call compute_cov(C_A_p,it,s_d2,mean_p,new_mean_p,c_p,clusters*adls)
        mean_p=new_mean_p
        if (it>burn) then
            C_rnd_W_p=C_A_p*factor_mu
            call choldc(C_rnd_W_p,clusters*adls)
        end if 
        if (it==burn) then
            factor_mu=1.0d0
        end if 
        call proposal(c_p,c_p_g,C_rnd_W_p,adls*clusters)
        do c_l=1,clusters;do v_l=1,adls
            ind=v_l+(c_l-1)*adls
            p_g(v_l,c_l)=min(max(c_p_g(ind,1),1.0d-6),0.9999d0)
        end do;end do
        !Evaluate posterior at new guess
        call posterior_fct_mu(p,sample_k,posterior_mu)
        call posterior_fct_mu(p_g,sample_k,posterior_mu_g)
        !Acceptance or rejection of new guess
        call RANDOM_NUMBER(u)
        if (u<min(exp(posterior_mu_g-posterior_mu),1.0d0)) then
            p=p_g
            call likelihood_all(p,likelihood) 
            chg_mu=chg_mu+1
        end if      
        !Print acceptance rate
        if ((it3==50 .and. it<1000) .or. (it3==500 .and. it>=1000) ) then
            print*,'completed ',it,' iterations'
            print*,'acceptance rate:'
            print*,'mu:',real(chg_mu)/real(it3)*100,'%'
            print*,'posterior:',real(posterior_mu)
            print*,' '
            if (real(chg_mu)/real(it3)*100<=15 ) then
                factor_mu=factor_mu/10.0_dp
                print*,'factor_mu',factor_mu
            elseif (real(chg_mu)/real(it3)*100>50 ) then
                factor_mu=factor_mu*5.0_dp
                print*,'factor_mu',factor_mu
            end if
            it3=0
            chg_mu=0
        end if
        
        !Sample h* given beta
        call sample_h_star(beta,y,sample_k)
        !Sample beta given h_star
        call sample_beta(beta)
        
        !Sample y* given gamma
        call sample_y_star(gamma,y,sample_k,y_star)
        !Sample gamma given y*
        call sample_gamma(y,y_star,gamma)
        !Save results
        call save_results(beta,gamma,p,posterior_mu,LE,fraction_t,it)
        
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