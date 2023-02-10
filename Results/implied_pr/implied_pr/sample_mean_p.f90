subroutine sample_mean_p(y,beta_var,beta_mean)
    use global_var; use mixtures_vars; use nrtype
    implicit none
    integer,dimension(indv,1),intent(in)::y
    real(DP),dimension(types,L_educ),intent(in)::beta_var
    real(DP),dimension(covariates_mix,types,L_educ),intent(out)::beta_mean
    integer::i_l,g_l,e_l,y_l,c_l
    real(DP)::age
    real(DP),dimension(covariates_mix,1)::x,z
    integer,dimension(types,L_educ)::counter_big_X
    real(DP),dimension(indv*10,types,L_educ,covariates_mix)::big_X
    real(DP),dimension(indv*10,types,L_educ)::big_Y
    real(DP),dimension(covariates_mix,covariates_mix)::Sigma,inv_Sigma,A
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    counter_big_X=0
    do i_l=1,indv;do g_l=first_age(i_l),last_age(i_l)
        age=initial_age+(g_l-1)*2-70
        x(1:4,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0/)     
        if (data_wealth(i_l,g_l)>0.0d0 .and. gender(i_l)==1 .and. initial_age+(g_l-1)*2<80) then
            counter_big_X(y(i_l,1),educ(i_l))=counter_big_X(y(i_l,1),educ(i_l))+1
            big_X(counter_big_X(y(i_l,1),educ(i_l)),y(i_l,1),educ(i_l),:)=x(:,1)
            big_Y(counter_big_X(y(i_l,1),educ(i_l)),y(i_l,1),educ(i_l))=log(data_wealth(i_l,g_l))
        end if
    end do; end do
    
    beta_mean=0.0_dp
    do y_l=1,types;do e_l=1,L_educ
        do c_l=1,covariates_mix
            z(c_l,1)=c4_normal_01(  )*sqrt(beta_var(y_l,e_l))
        end do
        Sigma=matmul(transpose(big_X(1:counter_big_X(y_l,e_l),y_l,e_l,:)),big_X(1:counter_big_X(y_l,e_l),y_l,e_l,:))
        call inverse(Sigma,inv_Sigma,covariates_mix)
        A=inv_Sigma
        call choldc(A,covariates_mix)
        beta_mean(:,y_l,e_l)=matmul(inv_Sigma,matmul(transpose(big_X(1:counter_big_X(y_l,e_l),y_l,e_l,:)),big_Y(1:counter_big_X(y_l,e_l),y_l,e_l)))+matmul(A,z(:,1))
    end do;end do
    
    
end subroutine
    
subroutine sample_mean_p_income(y,s2_w,u_draw,beta_mean)
    use global_var; use mixtures_vars_income; use nrtype
    implicit none
    integer,dimension(indv,1),intent(in)::y
    real(DP),dimension(L_educ,cohorts),intent(in)::s2_w
    real(DP),dimension(indv,generations),intent(in)::u_draw
    real(DP),dimension(covariates_mix_mean,L_educ),intent(out)::beta_mean
    integer::i_l,g_l,e_l,y_l,c_l
    real(DP)::age
    real(DP),dimension(covariates_mix_mean,1)::x,z
    integer,dimension(L_educ)::counter_big_X
    real(DP),dimension(indv*10,L_educ,covariates_mix_mean)::big_X
    real(DP),dimension(indv*10,L_educ)::big_Y
    real(DP),dimension(covariates_mix_mean,covariates_mix_mean)::Sigma,inv_Sigma,A
    real(DP),dimension(types)::y_d
    real(DP),dimension(cohorts)::cohort_d
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    counter_big_X=0
    do i_l=indv_HRS+1,indv;do g_l=first_age(i_l),last_age(i_l)
        age=initial_age+(g_l-1)*2-70
        y_d=0.0d0
        y_d(y(i_l,1))=1.0d0
        cohort_d=0.0d0
        cohort_d(birth_cohort(i_l))=1.0d0
        x(1:covariates_mix_mean,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0,dble(data_shlt(i_l,g_l)-1),y_d(2:types),cohort_d(4:5)/)     
        if ( gender(i_l)==1 .and. initial_age+(g_l-1)*2<63  .and. data_income(i_l,g_l)>520.0d0*7.25d0  ) then ! 
            counter_big_X(educ(i_l))=counter_big_X(educ(i_l))+1
            big_X(counter_big_X(educ(i_l)),educ(i_l),:)=x(:,1)/sqrt(s2_w(educ(i_l),birth_cohort(i_l)))
            big_Y(counter_big_X(educ(i_l)),educ(i_l))=(log(data_income(i_l,g_l))-u_draw(i_l,g_l))/sqrt(s2_w(educ(i_l),birth_cohort(i_l)))
            if (isnan((log(data_income(i_l,g_l))-u_draw(i_l,g_l))/sqrt(s2_w(educ(i_l),birth_cohort(i_l))))) then
                print*,''
            end if
        end if
    end do; end do
    
    beta_mean=0.0_dp
    do e_l=1,L_educ
        do c_l=1,covariates_mix_mean
            z(c_l,1)=c4_normal_01(  )
        end do
        Sigma=matmul(transpose(big_X(1:counter_big_X(e_l),e_l,:)),big_X(1:counter_big_X(e_l),e_l,:))
        call inverse(Sigma,inv_Sigma,covariates_mix_mean)
        A=inv_Sigma
        call choldc(A,covariates_mix_mean)
        beta_mean(:,e_l)=matmul(inv_Sigma,matmul(transpose(big_X(1:counter_big_X(e_l),e_l,:)),big_Y(1:counter_big_X(e_l),e_l)))+matmul(A,z(:,1))
    end do
    
    
end subroutine
    
    
    