module mixtures_vars
    use global_var
    implicit none
    integer,parameter::covariates_mix=4
    real(DP),dimension(indv,generations)::data_wealth
end module
    

subroutine estimate_mixture_wealth(type_pr)
    use global_var;use nrtype; use mixtures_vars
    implicit none
    double precision,dimension(indv,types),intent(in)::type_pr
    real(DP),dimension(generations,indv_HRS)::data_wealth_hrs
    real(DP),dimension(generations,indv_psid)::data_wealth_psid
    real(DP),dimension(generations,types,L_educ)::pr_zero
    real(DP),dimension(generations,types,L_educ,3)::quantile_w
    real(DP),dimension(3)::quantiles=(/0.25d0,0.5d0,0.75d0/)
    real(DP),dimension(covariates_mix,types,L_educ)::beta_mean
    real(DP),dimension(types,L_educ)::beta_var
    real(DP)::q,p
    integer::e_l,g_l,y_l,p_l
    character::pause_k

    
    open(unit=10,file=path//"Data\savings.csv")
        read(10,*) data_wealth_hrs
    close(10)
    
    open(unit=10,file=path//"Data\savings_psid.csv")
        read(10,*) data_wealth_psid
    close(10)
    
    data_wealth(1:indv_HRS,:)=reshape(data_wealth_hrs,(/indv_HRS,generations/),order=(/2,1/))
    data_wealth(indv_HRS+1:indv,:)=reshape(data_wealth_psid,(/indv_psid,generations/),order=(/2,1/))  

    pr_zero=-9.0d0
    quantile_w=-9.0d0
    call pr_of_zero_wealth(type_pr,pr_zero)
    call log_normal_dist(type_pr,beta_mean,beta_var)
    
    do p_l=1,3
        p=quantiles(p_l)
        do e_l=1,L_educ; do y_l=1,types;do g_l=1,generations
            if (pr_zero(g_l,y_l,e_l)>p) then 
                quantile_w(g_l,y_l,e_l,p_l)=0.0d0
            else
                q=p-pr_zero(g_l,y_l,e_l)
                call quantile_wealth_hat(y_l,e_l,g_l,beta_mean,beta_var,q,quantile_w(g_l,y_l,e_l,p_l))
            end if    
        end do; end do; end do
    end do
    
    open(unit=10,file=path//"metric_model\Results\wealth_moments_data.txt")
        do e_l=1,L_educ; do y_l=1,types;do g_l=1,generations
            write(10,'(I3,I3,I3,<3>F14.1)') y_l,e_l,g_l,quantile_w(g_l,y_l,e_l,1),quantile_w(g_l,y_l,e_l,2),quantile_w(g_l,y_l,e_l,3)
        end do; end do; end do
    close(10)

    
    open(unit=10,file=path//"metric_model\Results\wealth_moments_data.txt")
        do e_l=1,L_educ; do y_l=1,types;do g_l=1,generations
            write(10,'(I3,I3,I3,<3>F14.1)') y_l,e_l,g_l,quantile_w(g_l,y_l,e_l,1),quantile_w(g_l,y_l,e_l,2),quantile_w(g_l,y_l,e_l,3) 
        end do; end do; end do
    close(10)    

    

end subroutine