subroutine sample_s2w(y,u_draw,beta_mean,s2_w)
    use global_var;use nrtype; use mixtures_vars_income
    implicit none
    integer,dimension(indv,1),intent(in)::y
    double precision,dimension(indv,generations),intent(in)::u_draw
    double precision,dimension(covariates_mix_mean,L_educ),intent(in)::beta_mean
    double precision,dimension(L_educ,cohorts),intent(out)::s2_w
    real(DP),dimension(covariates_mix_mean,1)::x
    real(DP)::age
    integer::ind,t_l,i_l,e_l,c_l
    double precision,dimension(indv*generations,1)::w_draw
    double precision::v,s2,shape,scale
    real(DP),dimension(types)::y_d
    real(DP),dimension(cohorts)::cohort_d
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    interface
        double precision function r8_gamma_01_sample ( shape )
            implicit none
            double precision,intent(in)::shape
        end function r8_gamma_01_sample
    end interface
    
    !Sample s2w
    !do c_l=3,cohorts;do e_l=1,L_educ
        c_l=3
        e_l=1
        w_draw=-1.0d0
        ind=0
        do i_l=indv_HRS+1,indv
            do t_l=first_age(i_l),last_age(i_l)
                if (u_draw(i_l,t_l)/=-1.0d0 .and. gender(i_l)==1 .and. initial_age+(t_l-1)*2<63    .and. birth_cohort(i_l)==c_l .and. data_income(i_l,t_l)>520.0d0*7.25d0) then   ! .and. educ(i_l)==e_l
                    ind=ind+1
                    age=initial_age+(t_l-1)*2-70
                    y_d=0.0d0
                    y_d(y(i_l,1))=1.0d0
                    cohort_d=0.0d0
                    cohort_d(birth_cohort(i_l))=1.0d0
                    x(1:covariates_mix_mean,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0,dble(data_shlt(i_l,t_l)-1),y_d(2:types),cohort_d(4:5)/)
                    w_draw(ind,1)=log(data_income(i_l,t_l))-u_draw(i_l,t_l)-sum(x(1:covariates_mix_mean,1)*beta_mean(:,educ(i_l))) !data_income(i_l,:)
                end if
            end do
        end do
        v=dble(ind-1)
        s2=sum(w_draw(1:ind,1)**2)/v 
        shape=v/2
        scale=1/(v*s2/2)
        s2_w(:,:)=1/(r8_gamma_01_sample(shape)*scale)
    !end do;end do
    
    do c_l=1,2;do e_l=1,3
        s2_w(e_l,c_l)=s2_w(e_l,3)
    end do;end do
        




    
    
end subroutine