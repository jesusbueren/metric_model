subroutine log_normal_dist(type_pr,beta_mean,beta_var)
    use global_var;use nrtype; use mixtures_vars
    implicit none
    double precision,dimension(indv,types),intent(in)::type_pr
    integer,dimension(indv,1)::y
    real(DP),dimension(covariates_mix,types,L_educ),intent(out)::beta_mean
    real(DP),dimension(types,L_educ),intent(out)::beta_var
    integer,parameter::iterations=500
    real(DP),dimension(covariates_mix,types,L_educ,iterations)::beta_mean_it
    real(DP),dimension(types,L_educ,iterations)::beta_var_it
    integer::it
    

    beta_mean_it(:,:,:,1)=0.0d0
    beta_var_it(:,:,1)=5.0d0
    
    do it=1,iterations-1    
        print*,it
        !Sample health behavior type
        call sample_health_behavior(type_pr,y) 
    
        !Sample parameters of mean wealth
        call sample_mean_p(y,beta_var_it(:,:,it),beta_mean_it(:,:,:,it+1))
        
        !Sample parameters of var of wealth
        call sample_var_p(y,beta_mean_it(:,:,:,it+1),beta_var_it(:,:,it+1))

    end do

    beta_mean=sum(beta_mean_it,4)/dble(iterations)
    beta_var=sum(beta_var_it,3)/dble(iterations)
    

end subroutine
    
