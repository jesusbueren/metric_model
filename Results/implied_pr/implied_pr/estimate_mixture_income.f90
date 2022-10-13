module mixtures_vars_income
    use global_var
    implicit none
    integer,parameter::covariates_mix=4,covariates_mix_mean=9,covariates_mix_d=5
    real(DP),dimension(indv,generations)::data_income
end module
    
    
subroutine estimate_mixture_income(type_pr)
    use global_var;use nrtype; use mixtures_vars_income
    implicit none
    double precision,dimension(indv,types)::type_pr
    real(DP),dimension(generations,indv_psid)::data_income_psid
    
    
    open(unit=10,file=path//"Data\income_psid.csv")
        read(10,*) data_income_psid
    close(10)
    
    data_income=-9.0d0
    data_income(indv_HRS+1:indv,:)=reshape(data_income_psid,(/indv_psid,generations/),order=(/2,1/)) 

    !call pr_of_zero_income()
    !call pr_of_zero_income_dynamic()
    
    call income_process(type_pr)

end subroutine