subroutine sample_var_p(y,beta_mean,beta_var)
    use global_var; use mixtures_vars; use nrtype
    implicit none
    integer,dimension(indv,1),intent(in)::y
    real(DP),dimension(types,L_educ),intent(out)::beta_var
    real(DP),dimension(covariates_mix,types,L_educ),intent(in)::beta_mean
    integer::i_l,g_l,e_l,y_l,c_l
    real(DP)::age,n0=10,S0=1
    real(DP),dimension(covariates_mix,1)::x,z
    integer,dimension(types,L_educ)::counter_big_X
    real(DP),dimension(indv*10,types,L_educ)::big_u2
    real(SP)::shape,scale
    interface
        function gengam(a,r )
            implicit none
              real ( kind = 4 ) a
              real ( kind = 4 ) gengam
              real ( kind = 4 ) r
        end function gengam
    end interface
    
    counter_big_X=0
    do i_l=1,indv;do g_l=first_age(i_l),last_age(i_l)  
       age=initial_age+(g_l-1)*2-70
        x(1:4,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(age)**3.0d0/)     
        if (data_wealth(i_l,g_l)>0.0d0 .and. gender(i_l)==1) then
            counter_big_X(y(i_l,1),educ(i_l))=counter_big_X(y(i_l,1),educ(i_l))+1
            big_u2(counter_big_X(y(i_l,1),educ(i_l)),y(i_l,1),educ(i_l))=(log(data_wealth(i_l,g_l))-sum(x(:,1)*beta_mean(:,y(i_l,1),educ(i_l))))**2.0d0 !big_u2(1:500,1,1)
        end if
    end do; end do
    
    beta_var=0.0_dp
    do y_l=1,types;do e_l=1,L_educ
        shape=real(n0+counter_big_X(y_l,e_l))/2.0_sp
        scale=real(S0+sum(big_u2(1:counter_big_X(y_l,e_l),y_l,e_l)))/2.0_sp 
        beta_var(y_l,e_l)=1.0d0/gengam(scale,shape)
    end do;end do
    
    
    end subroutine
    
    
subroutine sample_var_p_income(y,beta_mean,beta_var)
    use global_var; use mixtures_vars_income; use nrtype
    implicit none
    integer,dimension(indv,1),intent(in)::y
    real(DP),dimension(types,L_educ),intent(out)::beta_var
    real(DP),dimension(covariates_mix,types,L_educ),intent(in)::beta_mean
    integer::i_l,g_l,e_l,y_l,c_l
    real(DP)::age,n0=10,S0=1
    real(DP),dimension(covariates_mix,1)::x,z
    integer,dimension(types,L_educ)::counter_big_X
    real(DP),dimension(indv*10,types,L_educ)::big_u2
    real(SP)::shape,scale
    interface
        function gengam(a,r )
            implicit none
              real ( kind = 4 ) a
              real ( kind = 4 ) gengam
              real ( kind = 4 ) r
        end function gengam
    end interface
    
    counter_big_X=0
    do i_l=indv_HRS+1,indv;do g_l=first_age(i_l),last_age(i_l)  
       age=initial_age+(g_l-1)*2-70
        x(1:4,1)=(/1.0_dp,dble(age),dble(age)**2.0d0,dble(data_shlt(i_l,g_l)-1)/)     
        if (data_income(i_l,g_l)>520.0d0*7.25d0 .and. gender(i_l)==1 .and. initial_age+(g_l-1)*2<62) then
            counter_big_X(y(i_l,1),educ(i_l))=counter_big_X(y(i_l,1),educ(i_l))+1
            big_u2(counter_big_X(y(i_l,1),educ(i_l)),y(i_l,1),educ(i_l))=(log(data_income(i_l,g_l))-sum(x(:,1)*beta_mean(:,y(i_l,1),educ(i_l))))**2.0d0 !big_u2(1:500,1,1)
        end if
    end do; end do
    
    beta_var=0.0_dp
    do y_l=1,types;do e_l=1,L_educ
        shape=real(n0+counter_big_X(y_l,e_l))/2.0_sp
        scale=real(S0+sum(big_u2(1:counter_big_X(y_l,e_l),y_l,e_l)))/2.0_sp 
        beta_var(y_l,e_l)=1.0d0/gengam(scale,shape)
    end do;end do
    
    
end subroutine
    
    
    