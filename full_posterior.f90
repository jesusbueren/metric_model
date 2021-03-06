subroutine full_posterior(beta,gamma,y)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates,types,clusters,clusters+1),intent(inout)::beta
    real(DP),dimension(covariates_habits,habits,types),intent(inout)::gamma
    integer,dimension(indv,1),intent(inout)::y
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H,H_g 
    real(DP),dimension(clusters,types,L_gender,L_educ)::init_cond 
    integer::c_l,c_l2,it3,it2,chg_mu,it,burn=100,v_l,ind
    real(DP)::prior_H
    integer,dimension(indv,generations)::sample_k
    real(DP)::u
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
    
    sample_k=data_shlt
    
    it2=-1
    it3=-1
    chg_mu=0
    !call tick(calc)
    do it=0,100000000
        print*,it
        it2=it2+1 
        it3=it3+1 
        !Sample health behavior type
        call sample_y(gamma,y,fraction_t)        
        !Sample h* given beta
        call sample_h_star(beta,y,sample_k)
        !Sample beta given h_star
        call sample_beta(beta)
        !Sample y* given gamma
        call sample_y_star(gamma,y,sample_k,y_star)
        !Sample gamma given y*
        call sample_gamma(y,y_star,gamma)
        !Compute share of indv in each group
        call compute_init(sample_k,y,init_cond)        
        !Compute life-expectancy
        call transitions(beta,init_cond,H,LE)
        !Save results
        call save_results(beta,gamma,LE,fraction_t,it)
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