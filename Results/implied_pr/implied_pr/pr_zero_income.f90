subroutine pr_of_zero_income(type_pr)
    use global_var;use nrtype; use mixtures_vars_income
    implicit none
    double precision,dimension(indv,types),intent(in)::type_pr
    integer,dimension(indv,1)::y
    real(DP),dimension(generations,clusters,L_educ,cohorts)::pr_zero
    real(DP),dimension(covariates_mix,clusters,L_educ)::beta_i
    integer,parameter::iterations=500
    real(DP),dimension(generations,clusters,L_educ,cohorts,iterations)::pr_zero_it
    integer::it,e_l,h_l,g_l,c_l,y_l
    

    beta_i=0.0d0
    
    do it=1,iterations    
        print*,it
        !Sample health behavior type
        call sample_health_behavior(type_pr,y) 
        !Sample probit parameters
        call sample_probit_p_income(y,beta_i)
        !Compute posterior for implied pr of zero
        call zero_pr_hat_income(beta_i,pr_zero_it(:,:,:,:,it)) 
    end do
    
    pr_zero=sum(pr_zero_it,5)/dble(iterations)
    
    open(unit=10,file=path//"metric_model\Results\labor_force_participation.txt")
        do c_l=1,cohorts;do e_l=1,L_educ; do h_l=1,clusters;do g_l=1,generations  
            write(10,'(<4>I3,F10.5)') h_l,e_l,g_l,c_l,1.0d0-pr_zero(g_l,h_l,e_l,c_l)
        end do; end do; end do;end do
    close(10)

end subroutine
    
