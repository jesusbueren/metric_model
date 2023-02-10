subroutine kalman_FS(s2_0_e,s2_nu_e,s2_w_e,beta_mean,y,rho_e,u_draw)
    use global_var;use nrtype; use mixtures_vars_income
    implicit none
    double precision,dimension(indv,generations),intent(out)::u_draw
    double precision,dimension(L_educ,cohorts),intent(in)::s2_nu_e
    double precision,dimension(L_educ,cohorts),intent(in)::s2_w_e
    double precision,dimension(L_educ,cohorts),intent(in)::s2_0_e
    double precision,dimension(L_educ,cohorts),intent(in)::rho_e
    real(DP),dimension(covariates_mix_mean,L_educ),intent(in)::beta_mean
    double precision::rho,s2_nu,s2_w,s2_0
    real(DP),dimension(covariates_mix_mean,1)::x
    integer,dimension(indv,1),intent(in)::y
    integer::t_l,i_l,t_l2
    real(DP),dimension(types)::y_d
    real(DP),dimension(cohorts)::cohort_d
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    real(DP)::K,J,age
    real(DP),dimension(generations+1)::xi_t1_t1,P_t1_t1,xi_t1_t0,P_t1_t0,xi_t1_T, &
                                               xi_t1_T_p,P_t1_T
    
    
    u_draw=-1.0d0
    do i_l=indv_HRS+1,indv
        rho=rho_e(educ(i_l),birth_cohort(i_l))
        s2_0=s2_0_e(educ(i_l),birth_cohort(i_l))
        s2_nu=s2_nu_e(educ(i_l),birth_cohort(i_l))
        s2_w=s2_w_e(educ(i_l),birth_cohort(i_l))
        !Filter
        xi_t1_t1=0.0d0;P_t1_t1=0.0d0;xi_t1_t0=0.0d0;P_t1_t0=0.0d0;xi_t1_T=0.0d0
        do t_l=first_age(i_l),last_age(i_l)
            if (t_l==first_age(i_l)) then
                xi_t1_t0(t_l)=0
                P_t1_t0(t_l)=s2_0
                if (t_l>1) then
                    do t_l2=2,first_age(i_l)
                        P_t1_t0(t_l)=rho**2*P_t1_t0(t_l)+s2_nu
                    end do
                end if
            end if
            age=initial_age+(t_l-1)*2-70
            y_d=0.0d0
            y_d(y(i_l,1))=1.0d0
            cohort_d=0.0d0
            cohort_d(birth_cohort(i_l))=1.0d0
            x(1:covariates_mix_mean,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0,dble(data_shlt(i_l,t_l)-1),y_d(2:types),cohort_d(4:5)/)

            !Gain & Updating equations
            if ( gender(i_l)==1 .and. initial_age+(t_l-1)*2<63 .and. data_income(i_l,t_l)>520.0d0*7.25d0 ) then !
                xi_t1_t1(t_l)=xi_t1_t0(t_l)+(P_t1_t0(t_l)*((P_t1_t0(t_l)+s2_w)**-1.0d0))*(log(data_income(i_l,t_l))-sum(x(1:covariates_mix_mean,1)*beta_mean(:,educ(i_l)))-xi_t1_t0(t_l)) 
            else
                K=0.0d0
                xi_t1_t1(t_l)=xi_t1_t0(t_l)
            end if

            P_t1_t1(t_l)=P_t1_t0(t_l)-P_t1_t0(t_l)*((P_t1_t0(t_l)+s2_w)**-1.0d0)*P_t1_t0(t_l)
            !Forecasting equations
            xi_t1_t0(t_l+1)=rho*xi_t1_t1(t_l)
            P_t1_t0(t_l+1)=rho**2*P_t1_t1(t_l)+s2_nu            
        end do
        !Smoother
        do t_l=last_age(i_l),first_age(i_l),-1
            if (t_l==last_age(i_l)) then
                xi_t1_T(t_l)=xi_t1_t1(t_l)
                u_draw(i_l,t_l)=xi_t1_t1(t_l)+c4_normal_01( )*sqrt(P_t1_t1(t_l))
            else
                J=P_t1_t1(t_l)*rho*((rho*P_t1_t1(t_l)*rho+s2_nu)**-1.0d0)
                xi_t1_T(t_l)=xi_t1_t1(t_l)+J*(u_draw(i_l,t_l+1)-xi_t1_t0(t_l+1))
                P_t1_T(t_l)=P_t1_t1(t_l)-J*rho*P_t1_t1(t_l)
                u_draw(i_l,t_l)=xi_t1_T(t_l)+c4_normal_01( )*sqrt(P_t1_T(t_l))
            end if                       
            
            if (isnan(u_draw(i_l,t_l))) then 
                print*,''
            end if
        end do  
    end do
    
    !u_draw(27091,1:17)
end subroutine