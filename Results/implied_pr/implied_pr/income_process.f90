subroutine income_process(type_pr)
    use global_var;use nrtype; use mixtures_vars_income
    implicit none
    double precision,dimension(indv,types),intent(in)::type_pr
    integer,dimension(indv,1)::y
    integer,parameter::burn=100,iterations=600!700
    real(DP),dimension(covariates_mix_mean,L_educ,iterations)::beta_mean
    real(DP),dimension(L_educ,cohorts,iterations)::s2_nu
    real(DP),dimension(L_educ,cohorts,iterations)::s2_w
    real(DP),dimension(L_educ,cohorts,iterations)::s2_0
    real(DP),dimension(L_educ,cohorts,iterations)::rho
    integer::it,e_l,y_l,g_l,c_l,s_l,g_l2
    real(DP),dimension(generations,types,L_educ,cohorts)::median_i,mean_i
    real(DP),dimension(indv,generations)::u_draw
    real(DP)::beta_var
    
    real(DP),dimension(covariates_mix_mean,L_educ)::beta_true
    real(DP),dimension(L_educ,cohorts)::s2_nu_true
    real(DP),dimension(L_educ,cohorts)::s2_w_true
    real(DP),dimension(L_educ,cohorts)::s2_0_true
    real(DP),dimension(L_educ,cohorts)::rho_true
    integer,parameter::sims=1!00
    
    real(DP),dimension(covariates_mix_mean,L_educ,sims)::beta_hat_s
    real(DP),dimension(L_educ,cohorts,sims)::s2_nu_hat_s,rho_hat_s,s2_w_hat_s,s2_0_hat_s

    
    
    do s_l=1,sims
        print*,s_l,'out of ',sims
        u_draw=-1.0d0 
        call sample_health_behavior(type_pr,y) 
    
        beta_true(:,1)=(/8.84794,	-0.0323242,	-0.000509,	7.53E-06,	-0.2686783,	-0.1,	-0.2,	0.0719227,	0.155986/) 
        beta_true(:,2)=(/9.84794,	-0.0323242,	-0.000509,	7.53E-06,	-0.2686783,	-0.1,	-0.2,	0.0719227,	0.155986/)
        beta_true(:,3)=(/10.84794,	-0.0323242,	-0.000509,	7.53E-06,	-0.2686783,	-0.1,	-0.2,	0.0719227,	0.155986/) 
        s2_w_true=0.01d0 
        s2_0_true=0.25d0 
        s2_nu_true=0.02d0 
        rho_true=0.9d0
    
       ! call simulate_income(y,beta_true,s2_nu_true,s2_w_true,s2_0_true,rho_true,u_draw)
    
        beta_mean(:,:,1)=0.0
        s2_w(:,:,1)=s2_w_true/2
        s2_0(:,:,1)=s2_0_true/4
        s2_nu(:,:,1)=s2_nu_true/3
        rho(:,:,1)=0.8d0!rho_true/2
        
       
        do it=1,iterations-1    
            !print*,'************************'
            !Sample health behavior type
            call sample_health_behavior(type_pr,y)  
    
            !Sample parameters of mean wealth
            call sample_mean_p_income(y,s2_w(:,:,it),u_draw,beta_mean(:,:,it+1)) 
        
            !sample shocks
            call kalman_FS(s2_0(:,:,it),s2_nu(:,:,it),s2_w(:,:,it),beta_mean(:,:,it+1),y,rho(:,:,it),u_draw) 
        
            !sample s2_0 and s2_nu
            call sample_rho_0_nu(u_draw,s2_nu(:,:,it),s2_0(:,:,it),rho(:,:,it),s2_0(:,:,it+1),s2_nu(:,:,it+1),rho(:,:,it+1)) !s2_nu(1,1,100:600)

            !sample s2_w    
            call sample_s2w(y,u_draw,beta_mean(:,:,it+1),s2_w(:,:,it+1)) 
            
        end do
    
        beta_hat_s(:,:,s_l)=sum(beta_mean(:,:,burn:iterations),3)/dble(iterations-burn) !beta_mean(1,1,100:200)
        s2_nu_hat_s(:,:,s_l)=sum(s2_nu(:,:,burn:iterations),3)/dble(iterations-burn)
        s2_w_hat_s(:,:,s_l)=sum(s2_w(:,:,burn:iterations),3)/dble(iterations-burn) !s2_w(1,1,:)
        s2_0_hat_s(:,:,s_l)=sum(s2_0(:,:,burn:iterations),3)/dble(iterations-burn)
        rho_hat_s(:,:,s_l)=sum(rho(:,:,burn:iterations),3)/dble(iterations-burn) !rho(:,:,600)
        print*,'beta',beta_hat_s(1,1,s_l)
        print*,'nu',s2_nu_hat_s(1,1,s_l)
        print*,'w',s2_w_hat_s(1,1,s_l)
        print*,'0',s2_0_hat_s(1,1,s_l)
        print*,'rho',rho_hat_s(1,1,s_l)
        
        if (s_l==1) then    
            open(unit=9,file=path//'metric_model\Results\montecarlo_beta.txt')
                    write(9,*) beta_hat_s(1,1,s_l)
            close(9)
            open(unit=9,file=path//'metric_model\Results\montecarlo_nu.txt')
                    write(9,*) s2_nu_hat_s(1,1,s_l)
            close(9)
            open(unit=9,file=path//'metric_model\Results\montecarlo_w.txt')
                    write(9,*) s2_w_hat_s(1,1,s_l)
            close(9)
            open(unit=9,file=path//'metric_model\Results\montecarlo_0.txt')
                    write(9,*) s2_0_hat_s(1,1,s_l)
            close(9)
            open(unit=9,file=path//'metric_model\Results\montecarlo_rho.txt')
                    write(9,*) rho_hat_s(1,1,s_l)
            close(9)
        else
            open(unit=9,file=path//'metric_model\Results\montecarlo_beta.txt',access="append")
                write(9,*) beta_hat_s(1,1,s_l)
            close(9)
            open(unit=9,file=path//'metric_model\Results\montecarlo_nu.txt',access="append")
                    write(9,*) s2_nu_hat_s(1,1,s_l)
            close(9)
            open(unit=9,file=path//'metric_model\Results\montecarlo_w.txt',access="append")
                    write(9,*) s2_w_hat_s(1,1,s_l)
            close(9)
            open(unit=9,file=path//'metric_model\Results\montecarlo_0.txt',access="append")
                    write(9,*) s2_0_hat_s(1,1,s_l)
            close(9)
            open(unit=9,file=path//'metric_model\Results\montecarlo_rho.txt',access="append")
                    write(9,*) rho_hat_s(1,1,s_l)
            close(9)
        end if    
    end do
    
    
    open(unit=9,file=path//'metric_model\Results\parameters_income.txt') !beta_mean(9,3,600)
            write(9,*) sum(beta_mean(:,:,burn:iterations),3)/dble(iterations-burn+1),sum(s2_nu(:,:,burn:iterations),3)/dble(iterations-burn+1),sum(s2_w(:,:,burn:iterations),3)/dble(iterations-burn+1),sum(rho(:,:,burn:iterations),3)/dble(iterations-burn+1),sum(s2_0(:,:,burn:iterations),3)/dble(iterations-burn+1)
    close(9)
    
    open(unit=10,file=path//"metric_model\Results\median_income.txt")
    do c_l=1,cohorts;do e_l=1,L_educ; do y_l=1,types;do g_l=1,generations
        beta_var=sum(s2_0(e_l,c_l,burn:iterations))/dble(iterations-burn+1)
        if (g_l>1) then
            do g_l2=2,g_l
                beta_var=(sum(rho(e_l,c_l,burn:iterations))/dble(iterations-burn+1))**2.0d0*beta_var+sum(s2_nu(e_l,c_l,burn:iterations))/dble(iterations-burn+1)
            end do
        end if
        beta_var=beta_var+sum(s2_w(e_l,c_l,burn:iterations))/dble(iterations-burn+1)
        print*,g_l,beta_var,(sum(rho(e_l,c_l,burn:iterations))/dble(iterations-burn+1))**2.0d0 !rho(1,1,100:it)
        call quantile_income_hat(y_l,e_l,c_l,g_l,sum(beta_mean(:,:,burn:iterations),3)/dble(iterations-burn+1),beta_var,0.5d0,median_i(g_l,y_l,e_l,c_l),mean_i(g_l,y_l,e_l,c_l))
        write(10,'(I3,I3,I3,I3,F15.2,F15.4,F15.4)') y_l,e_l,c_l,g_l,median_i(g_l,y_l,e_l,c_l),beta_var,mean_i(g_l,y_l,e_l,c_l)        
    end do; end do; end do;end do
    close(10)
    
    
    pause
end subroutine
    
