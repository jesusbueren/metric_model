subroutine sample_gamma_y(gamma,type_i,sample_k)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates_habits,habits,types),intent(inout)::gamma
    integer,dimension(indv,1),intent(in)::type_i
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(indv,habits,generations)::y_star
    real(DP),dimension(covariates_habits,1)::x
    integer::h_l,c_l,g_l,e_d,ge_l,age,it,i_l,health_d,e_l
    real(DP),dimension(indv*g_max,types,habits,covariates_habits)::big_X
    real(DP),dimension(indv*g_max,types,habits)::big_Y
    integer,dimension(types,habits)::counter_big_X
    real(DP),dimension(covariates_habits,1)::z
    real(DP),dimension(covariates_habits,covariates_habits)::Sigma,inv_Sigma,A
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    y_star=-9.0_sp
    counter_big_X=0
    do i_l=1,indv;
        do g_l=first_age(i_l),last_age(i_l);do h_l=1,habits
            health_d=sample_k(i_l,g_l)-1
            age=initial_age+(g_l-1)*2-70
            x(:,1)=(/1.0_dp,dble(age),dble(age**2.0_dp-1.0_dp)/)
            if (data_habits(i_l,h_l,g_l)==1 .or. data_habits(i_l,h_l,g_l)==0 ) then
                counter_big_X(type_i(i_l,1),h_l)=counter_big_X(type_i(i_l,1),h_l)+1
                big_X(counter_big_X(type_i(i_l,1),h_l),type_i(i_l,1),h_l,:)=x(:,1)
            end if
            if (data_habits(i_l,h_l,g_l)==1 ) then
                call TRUNCATED_NORMAL_A_SAMPLE(sum(x(:,1)*gamma(:,h_l,type_i(i_l,1))),1.0_dp,0.0_dp,y_star(i_l,h_l,g_l))
                big_Y(counter_big_X(type_i(i_l,1),h_l),type_i(i_l,1),h_l)=y_star(i_l,h_l,g_l)
            elseif (data_habits(i_l,h_l,g_l)==0 ) then
                call TRUNCATED_NORMAL_B_SAMPLE(sum(x(:,1)*gamma(:,h_l,type_i(i_l,1))),1.0_dp,0.0_dp,y_star(i_l,h_l,g_l))
                big_Y(counter_big_X(type_i(i_l,1),h_l),type_i(i_l,1),h_l)=y_star(i_l,h_l,g_l)
            end if
        end do; end do
    end do
    
    
    do e_l=1,types; do h_l=1,habits
        do c_l=1,covariates_habits
            z(c_l,1)=c4_normal_01(  )
        end do
        Sigma=matmul(transpose(big_X(1:counter_big_X(e_l,h_l),e_l,h_l,:)),big_X(1:counter_big_X(e_l,h_l),e_l,h_l,:))
        call inverse(Sigma,inv_Sigma,covariates_habits)
        A=inv_Sigma
        call choldc(A,covariates_habits)
        gamma(:,h_l,e_l)=matmul(inv_Sigma,matmul(transpose(big_X(1:counter_big_X(e_l,h_l),e_l,h_l,:)),big_Y(1:counter_big_X(e_l,h_l),e_l,h_l)))+matmul(A,z(:,1))
    end do; end do
    
end subroutine