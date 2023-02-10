subroutine pr_of_zero_income_dynamic()
    use global_var;use nrtype; use mixtures_vars_income
    implicit none
    integer,dimension(indv,1)::y
    real(DP),dimension(generations,clusters,L_educ,cohorts,2)::pr_zero
    real(DP),dimension(covariates_mix_d,clusters,L_educ)::beta_i
    integer,parameter::iterations=2000,burn=1000
    real(DP),dimension(generations,clusters,L_educ,cohorts,2,iterations)::pr_zero_it
    integer::it,e_l,h_l,g_l,lf_l,c_l
    

    beta_i=0.0d0
    
    do it=1,iterations    
        print*,it
        !Sample probit parameters
        call sample_probit_p_income_dynamic(y,beta_i)
        !Compute posterior for implied pr of zero
        call zero_pr_hat_income_dynamic(beta_i,pr_zero_it(:,:,:,:,:,it)) 
    end do
    
    pr_zero=sum(pr_zero_it(:,:,:,:,:,burn:iterations),6)/dble(iterations-burn+1)
    
    open(unit=10,file=path//"metric_model\Results\labor_force_participation_dynamic.txt")
         do lf_l=1,2; do c_l=1,cohorts;do e_l=1,L_educ; do h_l=1,clusters;do g_l=1,generations
            write(10,'(I3,I3,I3,I3,F10.5)') h_l,e_l,g_l,c_l,1.0d0-pr_zero(g_l,h_l,e_l,c_l,lf_l)
        end do; end do; end do; end do; end do
    close(10)
    

end subroutine
    
