subroutine income_process(type_pr)
    use global_var;use nrtype; use mixtures_vars_income
    implicit none
    double precision,dimension(indv,types),intent(in)::type_pr
    integer,dimension(indv,1)::y
    integer,parameter::iterations=500
    real(DP),dimension(covariates_mix_mean,types,L_educ,iterations)::beta_mean
    real(DP),dimension(types,L_educ)::beta_var
    real(DP),dimension(iterations)::rho,s2_nu, s2_w
    integer::it,e_l,y_l,g_l
    real(DP),dimension(generations,types,L_educ)::median_i
    real(DP),dimension(indv,generations)::u_draw
    
    

    beta_mean(:,:,:,1)=0.0d0
    s2_w(1)=0.25d0
    rho(1)=0.9d0
    s2_nu(1)=0.25d0*(1.0d0-rho(1)**2)    
    u_draw=0.0d0
    
    
    do it=1,iterations-1    
        print*,it,rho(it)
        !Sample health behavior type
        call sample_health_behavior(type_pr,y) 
    
        !Sample parameters of mean wealth
        beta_var=s2_w(it)
        call sample_mean_p_income(y,beta_var,u_draw,beta_mean(:,:,:,it+1))
        
        !sample shocks
        call kalman_FS(rho(it),s2_nu(it),s2_w(it),beta_mean(:,:,:,it+1),y,u_draw)
        
        !sample rho and s2_nu
        call sample_rho_nu(u_draw,s2_nu(it),rho(it+1),s2_nu(it+1))
        
        !sample s2_w
        call sample_s2w(y,u_draw,beta_mean(:,:,:,it+1),s2_w(it+1)) 

    end do
    
    open(unit=9,file=path//'metric_model\Results\parameters_income.txt')
            write(9,*) sum(beta_mean,4)/dble(iterations),sum(rho)/dble(iterations),sum(s2_nu)/dble(iterations),sum(s2_w)/dble(iterations)
    close(9)

    open(unit=10,file=path//"metric_model\Results\median_income.txt")
    do e_l=1,L_educ; do y_l=1,types;do g_l=1,generations
            beta_var=sum(s2_w+s2_nu/(1-rho**2))/dble(iterations)
            call quantile_income_hat(y_l,e_l,g_l,sum(beta_mean,4)/dble(iterations),beta_var,0.5d0,median_i(g_l,y_l,e_l))
           write(10,'(I3,I3,I3,F15.2)') y_l,e_l,g_l,median_i(g_l,y_l,e_l)
    end do; end do; end do
    close(10)
    

end subroutine
    
