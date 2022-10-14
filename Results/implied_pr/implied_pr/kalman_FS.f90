subroutine kalman_FS(rho_e,s2_nu_e,s2_w_e,beta_mean,y,u_draw)
    use global_var;use nrtype; use mixtures_vars_income
    implicit none
    double precision,dimension(indv,generations),intent(out)::u_draw
    double precision,dimension(L_educ),intent(in)::rho_e,s2_nu_e,s2_w_e
    real(DP),dimension(covariates_mix_mean,L_educ),intent(in)::beta_mean
    double precision::rho,s2_nu,s2_w
    real(DP),dimension(covariates_mix_mean,1)::x
    integer,dimension(indv,1),intent(in)::y
    integer::t_l,i_l
    real(DP),dimension(types)::y_d
    real(DP),dimension(cohorts)::cohort_d
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    real(DP)::K,J,age
    real(DP),dimension(generations+1)::xi_t1_t1,P_t1_t1,xi_t1_t0,P_t1_t0,xi_t1_T, &
                                               xi_t1_T_p,xi_t1_t1_p,xi_t1_t0_p,xi_p,y_p
    
    
    u_draw=-1.0d0
    do i_l=indv_HRS+1,indv
        rho=rho_e(educ(i_l))
        s2_nu=s2_nu_e(educ(i_l))
        s2_w=s2_w_e(educ(i_l))
        !Filter
        xi_t1_t1=0.0d0;P_t1_t1=0.0d0;xi_t1_t0=0.0d0;P_t1_t0=0.0d0;xi_t1_T=0.0d0
        xi_t1_T_p=0.0d0;xi_t1_t1_p=0.0d0;xi_t1_t0_p=0.0d0;xi_p=0.0d0;y_p=0.0d0
        do t_l=first_age(i_l),last_age(i_l)
            if (t_l==first_age(i_l)) then
                xi_t1_t0(t_l)=0
                P_t1_t0(t_l)=s2_nu/(1-rho**2)
                xi_t1_t0_p(t_l)=0
                xi_p(t_l)=sqrt(s2_nu/(1-rho**2))*c4_normal_01( ) !simulated draw
            else
                xi_p(t_l)=rho*xi_p(t_l-1)+sqrt(s2_nu)*c4_normal_01( ) !simulated draw
            end if
            age=initial_age+(t_l-1)*2-70
            y_d=0.0d0
            y_d(y(i_l,1))=1.0d0
            cohort_d=0.0d0
            cohort_d(birth_cohort(i_l))=1.0d0
            x(1:covariates_mix_mean,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0,dble(data_shlt(i_l,t_l)-1),y_d(2:types),cohort_d(4:5)/)
            y_p(t_l)=sum(x(1:covariates_mix_mean,1)*beta_mean(:,educ(i_l)))+xi_p(t_l)+c4_normal_01( )*sqrt(s2_w)
            !Gain & Updating equations
            if (data_income(i_l,t_l)>520.0d0*7.25d0 .and. gender(i_l)==1 .and. initial_age+(t_l-1)*2<60 ) then
                K=P_t1_t0(t_l)*(P_t1_t0(t_l)+s2_w)**-1.0d0
                xi_t1_t1(t_l)=xi_t1_t0(t_l)+K*(log(data_income(i_l,t_l))-sum(x(1:covariates_mix_mean,1)*beta_mean(:,educ(i_l)))-xi_t1_t0(t_l)) 
            else
                K=0.0d0
                xi_t1_t1(t_l)=xi_t1_t0(t_l)
            end if
            xi_t1_t1_p(t_l)=xi_t1_t0_p(t_l)+K*(y_p(t_l)-sum(x(1:covariates_mix_mean,1)*beta_mean(:,educ(i_l)))-xi_t1_t0_p(t_l))
            P_t1_t1(t_l)=P_t1_t0(t_l)-K*P_t1_t0(t_l)
            !Forecasting equations
            xi_t1_t0(t_l+1)=rho*xi_t1_t1(t_l)
            xi_t1_t0_p(t_l+1)=rho*xi_t1_t1_p(t_l)
            P_t1_t0(t_l+1)=rho*P_t1_t1(t_l)*rho+s2_nu            
        end do
        !Smoother
        do t_l=last_age(i_l),first_age(i_l),-1
            if (t_l==last_age(i_l)) then
                xi_t1_T(t_l)=xi_t1_t1(t_l)
                xi_t1_T_p(t_l)=xi_t1_t1_p(t_l)
            else
                J=P_t1_t1(t_l)*rho*(P_t1_t0(t_l+1)**-1)
                xi_t1_T(t_l)=xi_t1_t1(t_l)+J*(xi_t1_T(t_l+1)-xi_t1_t0(t_l+1))
                xi_t1_T_p(t_l)=xi_t1_t1_p(t_l)+J*(xi_t1_T_p(t_l+1)-xi_t1_t0_p(t_l+1))
            end if                       
            u_draw(i_l,t_l)=xi_t1_T(t_l)-xi_t1_T_p(t_l)+xi_p(t_l) 
            if (isnan(u_draw(i_l,t_l))) then
                print*,''
            end if
        end do 
    end do

end subroutine