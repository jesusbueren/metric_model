subroutine zero_pr_hat(beta_w,pr_zero)
    use global_var; use nrtype; use mixtures_vars
    implicit none
    real(DP),dimension(covariates_mix,types,L_educ),intent(in)::beta_w
    real(DP),dimension(generations,types,L_educ),intent(out)::pr_zero
    real(DP),dimension(covariates_mix,1)::x
    integer::y_l,e_l,g_l,age
    
    do y_l=1,types;do e_l=1,L_educ; do g_l=1,generations
        age=initial_age+(g_l-1)*2-70
        x(1:4,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0/) 
        pr_zero(g_l,y_l,e_l)=1.0_dp-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*beta_w(:,y_l,e_l))/sqrt(2.0_dp)))
    end do;end do;end do
    
end subroutine
    
subroutine zero_pr_hat_income(beta_i,pr_zero)
    use global_var; use nrtype; use mixtures_vars_income
    implicit none
    real(DP),dimension(covariates_mix,L_educ),intent(in)::beta_i
    real(DP),dimension(generations,clusters,L_educ,cohorts),intent(out)::pr_zero
    real(DP),dimension(covariates_mix,1)::x
    integer::h_l,e_l,g_l,age,c_l,y_l
    real(DP),dimension(types)::y_d
    real(DP),dimension(cohorts)::cohort_d
    
    do h_l=1,clusters;do e_l=1,L_educ; do g_l=1,generations;do c_l=1,cohorts
        age=initial_age+(g_l-1)*2-70

        cohort_d=0.0d0
        cohort_d(c_l)=1.0d0
        x(1:covariates_mix,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(h_l-1),cohort_d(4:5)/)    
        pr_zero(g_l,h_l,e_l,c_l)=1.0_dp-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*beta_i(:,e_l))/sqrt(2.0_dp)))
    end do;end do;end do;end do
    
end subroutine  
    
subroutine zero_pr_hat_income_dynamic(beta_i,pr_zero)
    use global_var; use nrtype; use mixtures_vars_income
    implicit none
    real(DP),dimension(covariates_mix_d,L_educ),intent(in)::beta_i
    real(DP),dimension(generations,clusters,L_educ,cohorts,2),intent(out)::pr_zero
    real(DP),dimension(covariates_mix_d,1)::x
    integer::h_l,e_l,g_l,age,lf_l,c_l
    real(DP),dimension(cohorts)::cohort_d
    
    do h_l=1,clusters;do e_l=1,L_educ; do g_l=1,generations; do lf_l=1,2;do c_l=1,cohorts
        age=initial_age+(g_l-1)*2-70
        cohort_d=0.0d0
        cohort_d(c_l)=1.0d0
        x(1:covariates_mix_d,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(h_l-1),dble(lf_l-1),cohort_d(4:5)/) 
        pr_zero(g_l,h_l,e_l,c_l,lf_l)=1.0_dp-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*beta_i(:,e_l))/sqrt(2.0_dp)))
    end do;end do;end do;end do;end do
    
end subroutine  
    