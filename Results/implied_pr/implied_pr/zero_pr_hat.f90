subroutine zero_pr_hat(beta_w,pr_zero)
    use global_var; use nrtype; use mixtures_vars
    implicit none
    real(DP),dimension(covariates_mix,types,L_educ),intent(in)::beta_w
    real(DP),dimension(generations,types,L_educ,clusters,cohorts),intent(out)::pr_zero
    real(DP),dimension(covariates_mix,1)::x
    real(DP),dimension(cohorts)::cohort_d
    integer::y_l,e_l,g_l,age,h_l,c_l
    real(DP)::h_d
    
    do y_l=1,types;do e_l=1,L_educ;do h_l=1,clusters; do g_l=1,generations;do c_l=1,cohorts
        h_d=dble(h_l-1)
        age=initial_age+(g_l-1)*2-70
        cohort_d=0.0d0
        cohort_d(c_l)=1.0d0
        x(1:covariates_mix,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0,cohort_d(2:cohorts)/) 
        pr_zero(g_l,y_l,e_l,h_l,c_l)=1.0_dp-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*beta_w(:,y_l,e_l))/sqrt(2.0_dp)))
    end do;end do;end do;end do;end do
    
end subroutine
    
subroutine zero_pr_hat_income(beta_i,pr_zero)
    use global_var; use nrtype; use mixtures_vars_income
    implicit none
    real(DP),dimension(covariates_mix,L_educ),intent(in)::beta_i
    real(DP),dimension(generations,types,L_educ,clusters),intent(out)::pr_zero
    real(DP),dimension(covariates_mix,1)::x
    integer::h_l,e_l,g_l,age,c_l,y_l
    real(DP),dimension(types)::y_d
    real(DP),dimension(cohorts)::cohort_d
    
    do h_l=1,clusters;do e_l=1,L_educ; do g_l=1,generations;do y_l=1,types
        age=initial_age+(g_l-1)*2-70
        y_d=0.0d0
        y_d(y_l)=1.0d0
        x(1:covariates_mix,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0,dble(h_l-1),dble(h_l-1)*dble(age)/)    
        pr_zero(g_l,y_l,e_l,h_l)=1.0_dp-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*beta_i(:,e_l))/sqrt(2.0_dp)))
    end do;end do;end do;end do
    
end subroutine  
    
subroutine zero_pr_hat_income_dynamic(beta_i,pr_zero)
    use global_var; use nrtype; use mixtures_vars_income
    implicit none
    real(DP),dimension(covariates_mix_d,L_educ),intent(in)::beta_i
    real(DP),dimension(generations,types,L_educ,clusters,2),intent(out)::pr_zero
    real(DP),dimension(covariates_mix_d,1)::x
    integer::h_l,e_l,g_l,age,c_l,y_l,f_l
    real(DP),dimension(types)::y_d
    real(DP),dimension(cohorts)::cohort_d
    
    do h_l=1,clusters;do e_l=1,L_educ; do g_l=1,generations;do y_l=1,types;do f_l=1,2
        age=initial_age+(g_l-1)*2-70
        x(1:covariates_mix_d,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0,dble(h_l-1),dble(h_l-1)*dble(age),dble(f_l-1),dble(f_l-1)*dble(age)/)    
        pr_zero(g_l,y_l,e_l,h_l,f_l)=1.0_dp-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*beta_i(:,e_l))/sqrt(2.0_dp)))
    end do;end do;end do;end do;end do
    
end subroutine  
    