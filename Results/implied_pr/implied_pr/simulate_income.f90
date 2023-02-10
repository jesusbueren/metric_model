subroutine simulate_income(y,beta_mean,s2_nu,s2_w,s2_0,rho,u_draw)
    use global_var;use nrtype; use mixtures_vars_income
    implicit none
    real(DP),dimension(covariates_mix_mean,L_educ),intent(in)::beta_mean
    real(DP),dimension(L_educ,cohorts),intent(in)::s2_nu
    real(DP),dimension(L_educ,cohorts),intent(in)::s2_w
    real(DP),dimension(L_educ,cohorts),intent(in)::s2_0
    real(DP),dimension(L_educ,cohorts),intent(in)::rho
    real(DP),dimension(types)::y_d
    real(DP),dimension(cohorts)::cohort_d
    real(DP),dimension(indv,generations),intent(inout)::u_draw
    integer::i_l,g_l,ind,g_l2
    real(DP)::age,u_m,var_ind
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
        age=initial_age+(g_l-1)*2-70
        y_d=0.0d0
        y_d(y(i_l,1))=1.0d0
        cohort_d=0.0d0
        cohort_d(birth_cohort(i_l))=1.0d0
        x(1:covariates_mix_mean,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0,dble(data_shlt(i_l,g_l)-1),y_d(2:types),cohort_d(4:5)/)
        if (g_l==first_age(i_l)) then
            var_ind=s2_0(educ(i_l),birth_cohort(i_l))
            if (g_l>1) then
                do g_l2=2,first_age(i_l)
                    var_ind=rho(educ(i_l),birth_cohort(i_l))**2*var_ind+s2_nu(educ(i_l),birth_cohort(i_l))
                end do
            end if
            u_draw(i_l,g_l)= c4_normal_01( )*sqrt(var_ind)
        else
            u_draw(i_l,g_l)=rho(educ(i_l),birth_cohort(i_l))*u_draw(i_l,g_l-1)+c4_normal_01( )*sqrt(s2_nu(educ(i_l),birth_cohort(i_l)))
        end if
        if (gender(i_l)==1 .and. initial_age+(g_l-1)*2<63) then
            ind=ind+1
            w_draw(ind,1)=c4_normal_01( )*sqrt(s2_w(educ(i_l),birth_cohort(i_l))) 
            data_income(i_l,g_l)=exp(sum(beta_mean(:,educ(i_l))*x(:,1))+u_draw(i_l,g_l)+w_draw(ind,1))
        end if
        if (isnan(data_income(i_l,g_l))) then
            print*,sum(beta_mean(:,educ(i_l))*x(:,1)),exp(sum(beta_mean(:,educ(i_l))*x(:,1))),u_draw(i_l,g_l),c4_normal_01( )**sqrt(s2_w(educ(i_l),birth_cohort(i_l)))
            pause
        end if
    end do; end do
!data_income(27091:27200,:)
end subroutine 