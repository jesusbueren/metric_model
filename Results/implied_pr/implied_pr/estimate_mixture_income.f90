module mixtures_vars_income
    use global_var
    implicit none
    integer,parameter::covariates_mix=6,covariates_mix_mean=6,covariates_mix_d=8
    real(DP),dimension(indv,generations)::data_income
end module
    
    
subroutine estimate_mixture_income(type_pr,sample_k)
    use global_var;use nrtype; use mixtures_vars_income; use mixtures_vars
    implicit none
    double precision,dimension(indv,types),intent(in)::type_pr
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(generations,indv_psid)::data_income_psid
    
    
    open(unit=10,file=path//"Data\income_psid.csv")
        read(10,*) data_income_psid
    close(10)

    data_income=-9.0d0
    data_income(indv_HRS+1:indv,:)=reshape(data_income_psid,(/indv_psid,generations/),order=(/2,1/)) 

    call pr_of_zero_income(type_pr,sample_k)
    call pr_of_zero_income_dynamic(type_pr,sample_k)
    call income_process(type_pr)

end subroutine