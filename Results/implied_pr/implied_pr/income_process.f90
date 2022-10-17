subroutine income_process(type_pr)
    use global_var;use nrtype; use mixtures_vars_income
    implicit none
    double precision,dimension(indv,types),intent(in)::type_pr
    integer,dimension(indv,1)::y
    integer,parameter::iterations=5000,burn=1000
    real(DP),dimension(covariates_mix_mean,L_educ,iterations)::beta_mean
    real(DP),dimension(L_educ,iterations)::s2_nu,s2_w
    real(DP),dimension(cohorts,L_educ,iterations)::s2_0
    integer::it,e_l,y_l,g_l,c_l
    real(DP),dimension(generations,types,L_educ,cohorts)::median_i
    real(DP),dimension(indv,generations)::u_draw
    real(DP)::beta_var
    
    

    beta_mean(:,:,1)=0.0d0
    s2_w(:,1)=0.25d0
    s2_0(:,:,1)=0.25d0
    s2_nu(:,1)=0.02d0
    u_draw=0.0d0
    
    
    do it=1,iterations-1    
        print*,'it',it
        print*,'dp',s2_0(3:5,1,it)
        print*,'hs',s2_0(3:5,2,it)
        print*,'cl',s2_0(3:5,3,it)
        print*,'************************'
        !Sample health behavior type
        call sample_health_behavior(type_pr,y) 
    
        !Sample parameters of mean wealth
        call sample_mean_p_income(y,s2_w(:,it),u_draw,beta_mean(:,:,it+1))
        
        !sample shocks
        call kalman_FS(s2_0(:,:,it),s2_nu(:,it),s2_w(:,it),beta_mean(:,:,it+1),y,u_draw)
        
        !sample s2_0 and s2_nu
        call sample_0_nu(u_draw,s2_nu(:,it),s2_0(:,:,it),s2_0(:,:,it+1),s2_nu(:,it+1)) 
        
        !sample s2_w    
        call sample_s2w(y,u_draw,beta_mean(:,:,it+1),s2_w(:,it+1)) 

    end do
    
    open(unit=9,file=path//'metric_model\Results\parameters_income.txt')
            write(9,*) sum(beta_mean(:,:,burn:iterations),3)/dble(iterations-burn+1),sum(s2_nu)/dble(iterations-burn+1),sum(s2_nu)/dble(iterations-burn+1),sum(s2_w)/dble(iterations-burn+1)
    close(9)

    open(unit=10,file=path//"metric_model\Results\median_income.txt")
    do c_l=1,cohorts;do e_l=1,L_educ; do y_l=1,types;do g_l=1,generations
            beta_var=s2_0(c_l,e_l,iterations)+(g_l-1)*s2_nu(e_l,iterations)+s2_w(e_l,iterations)
            call quantile_income_hat(y_l,e_l,c_l,g_l,sum(beta_mean,3)/dble(iterations),beta_var,0.5d0,median_i(g_l,y_l,e_l,c_l))
           write(10,'(I3,I3,I3,I3,F15.2,F10.4)') y_l,e_l,c_l,g_l,median_i(g_l,y_l,e_l,c_l),beta_var
    end do; end do; end do;end do
    close(10)
    
    pause
end subroutine
    
