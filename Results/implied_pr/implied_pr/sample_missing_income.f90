subroutine sample_missing_income(y,s2_w,u_draw,beta_mean) 
    use global_var;use nrtype; use mixtures_vars_income
    implicit none
    real(DP),dimension(covariates_mix_mean,L_educ),intent(in)::beta_mean
    real(DP),dimension(L_educ,cohorts),intent(in)::s2_w
    real(DP),dimension(types)::y_d
    real(DP),dimension(cohorts)::cohort_d
    real(DP),dimension(indv,generations),intent(in)::u_draw
    integer::i_l,g_l,ind
    real(DP)::age
    real(DP),dimension(covariates_mix_mean,1)::x
    double precision,dimension(indv*generations,1)::w_draw
    integer,dimension(indv,1),intent(in)::y
     interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
     end interface
    
     first_age=1
     last_age=17
     !data_income=0.0d0

    ind=0
    do i_l=indv_HRS+1,indv;do g_l=first_age(i_l),last_age(i_l)
        if (data_income(i_l,g_l)<520.0d0*7.25d0 ) then
            age=initial_age+(g_l-1)*2-70
            y_d=0.0d0
            y_d(y(i_l,1))=1.0d0
            cohort_d=0.0d0
            cohort_d(birth_cohort(i_l))=1.0d0
            x(1:covariates_mix_mean,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0,dble(data_shlt(i_l,g_l)-1),y_d(2:types),cohort_d(4:5)/)
            if ( gender(i_l)==1 .and. initial_age+(g_l-1)*2<63) then
                ind=ind+1
    1            w_draw(ind,1)=c4_normal_01( )*sqrt(s2_w(educ(i_l),birth_cohort(i_l)))
                if (exp(sum(beta_mean(:,educ(i_l))*x(:,1))+u_draw(i_l,g_l)+w_draw(ind,1))<520.0d0*7.25d0 ) then
                    data_income(i_l,g_l)=exp(sum(beta_mean(:,educ(i_l))*x(:,1))+u_draw(i_l,g_l)+w_draw(ind,1))
                else
                    go to 1
                end if
            end if
        end if
    end do; end do

end subroutine 