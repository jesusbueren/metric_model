subroutine pr_of_zero_wealth(type_pr,pr_zero)
    use global_var;use nrtype; use mixtures_vars
    double precision,dimension(indv,types),intent(in)::type_pr
    integer,dimension(indv,1)::y
    real(DP),dimension(generations,types,L_educ),intent(out)::pr_zero
    real(DP),dimension(covariates_mix,types,L_educ)::beta_w
    integer,parameter::iterations=500
    real(DP),dimension(generations,types,L_educ,iterations)::pr_zero_it
    integer::it
    

    beta_w=0.0d0
    
    do it=1,iterations    
        print*,it
        !Sample health behavior type
        call sample_health_behavior(type_pr,y) 
    
        !Sample probit parameters
        call sample_probit_p(y,beta_w)
        
        !Compute posterior for implied pr of zero
        call zero_pr_hat(beta_w,pr_zero_it(:,:,:,it))
    
    end do
    
    pr_zero=sum(pr_zero_it,4)/dble(iterations)
    

end subroutine
    
