subroutine quantile_wealth_hat(y_l,e_l,g_l,beta_mean,beta_var,quantile,q_wealth,mean_w,variance_w)
    use global_var; use nrtype; use mixtures_vars
    implicit none
    integer,intent(in)::y_l,e_l,g_l
    real(DP),dimension(covariates_mix,types,L_educ),intent(in)::beta_mean
    real(DP),dimension(types,L_educ),intent(in)::beta_var
    real(DP),intent(in)::quantile
    real(DP),intent(out)::q_wealth,mean_w,variance_w
    real(DP),dimension(covariates_mix,1)::x
    integer::age
    real(DP)::x2


        age=initial_age+(g_l-1)*2-70
        x(1:4,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0/) 
        x2=2.0d0*quantile-1.0d0
        q_wealth=exp(sum(x(:,1)*beta_mean(:,y_l,e_l))+sqrt(beta_var(y_l,e_l))*sqrt(2.0d0)*&
            sqrt(pi_d)*(0.5d0*x2+1.0d0/24.0d0*pi_d*x2**3.0d0+7.0d0/960.0d0*pi_d**2.0d0*x2**5.0d0+127.0d0/80640.0d0*pi_d**3.0d0*x2**7.0d0))
        mean_w=exp(sum(x(:,1)*beta_mean(:,y_l,e_l))+beta_var(y_l,e_l)/2.0d0)
        variance_w=(exp(beta_var(y_l,e_l))-1.0d0)*exp(2.0d0*sum(x(:,1)*beta_mean(:,y_l,e_l))+beta_var(y_l,e_l))
        
    
end subroutine
    
    
subroutine quantile_income_hat(y_l,e_l,g_l,beta_mean,beta_var,quantile,q_income)
    use global_var; use nrtype; use mixtures_vars_income
    implicit none
    integer,intent(in)::y_l,e_l,g_l
    real(DP),dimension(covariates_mix_mean,types,L_educ),intent(in)::beta_mean
    real(DP),dimension(types,L_educ),intent(in)::beta_var
    real(DP),intent(in)::quantile
    real(DP),intent(out)::q_income
    real(DP),dimension(covariates_mix_mean,1)::x
    integer::age
    real(DP)::x2


        age=initial_age+(g_l-1)*2-70
        x(1:covariates_mix_mean,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0,0.0d0/) 
        x2=2.0d0*quantile-1.0d0
        q_income=exp(sum(x(:,1)*beta_mean(:,y_l,e_l))+sqrt(beta_var(y_l,e_l))*sqrt(2.0d0)*&
            sqrt(pi_d)*(0.5d0*x2+1.0d0/24.0d0*pi_d*x2**3.0d0+7.0d0/960.0d0*pi_d**2.0d0*x2**5.0d0+127.0d0/80640.0d0*pi_d**3.0d0*x2**7.0d0))


    
end subroutine