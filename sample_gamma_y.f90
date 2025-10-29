subroutine sample_gamma_y(gamma,type_i,sample_k)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates_habits,habits_nomed,types),intent(inout)::gamma
    integer,dimension(indv,1),intent(in)::type_i
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(indv,habits_nomed,generations)::y_star
    real(DP),dimension(covariates_habits,1)::x
    integer::h_l,c_l,g_l,e_d,ge_l,age,it,i_l,e_l
    real(DP)::health_d
    real(DP),dimension(indv*g_max,types,habits_nomed,covariates_habits)::big_X
    real(DP),dimension(indv*g_max,types,habits_nomed)::big_Y
    integer,dimension(types,habits)::counter_big_X
    real(DP),dimension(covariates_habits,1)::z
    real(DP),dimension(covariates_habits,covariates_habits)::Sigma,inv_Sigma,A,B_0
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    y_star=-9.0_sp
    counter_big_X=0
    do i_l=1,indv;
        if (sample_selection(i_l)) then
            do g_l=first_age(i_l),last_age(i_l);do h_l=1,habits_nomed
                if (sample_k(i_l,g_l)/=-1) then
                    health_d=dble(sample_k(i_l,g_l)-1)
                    age=initial_age+(g_l-1)*2
                    x(:,1)=(/1.0_dp,dble(age),dble(age**2.0_dp-1.0_dp),health_d/)
                    if ((data_habits(i_l,habits_vec(h_l),g_l)==1 .or. data_habits(i_l,habits_vec(h_l),g_l)==0) ) then
                        counter_big_X(type_i(i_l,1),h_l)=counter_big_X(type_i(i_l,1),h_l)+1
                        big_X(counter_big_X(type_i(i_l,1),h_l),type_i(i_l,1),h_l,:)=x(:,1)
                    end if
                    if (data_habits(i_l,habits_vec(h_l),g_l)==1 ) then
                        call TRUNCATED_NORMAL_A_SAMPLE(sum(x(:,1)*gamma(:,h_l,type_i(i_l,1))),1.0_dp,0.0_dp,y_star(i_l,h_l,g_l))
                        big_Y(counter_big_X(type_i(i_l,1),h_l),type_i(i_l,1),h_l)=y_star(i_l,h_l,g_l)
                    elseif (data_habits(i_l,habits_vec(h_l),g_l)==0 ) then
                        call TRUNCATED_NORMAL_B_SAMPLE(sum(x(:,1)*gamma(:,h_l,type_i(i_l,1))),1.0_dp,0.0_dp,y_star(i_l,h_l,g_l))
                        big_Y(counter_big_X(type_i(i_l,1),h_l),type_i(i_l,1),h_l)=y_star(i_l,h_l,g_l)
                    end if
                end if
            end do; end do
        end if
    end do
    
    
    do e_l=1,types; do h_l=1,habits_nomed
        if (counter_big_X(e_l,h_l)>1)then
            B_0=0.0d0
            do c_l=1,covariates_habits
                z(c_l,1)=c4_normal_01(  )
                B_0(c_l,c_l)=0.0d0
            end do
            Sigma=B_0+matmul(transpose(big_X(1:counter_big_X(e_l,h_l),e_l,h_l,:)),big_X(1:counter_big_X(e_l,h_l),e_l,h_l,:))
            call inverse(Sigma,inv_Sigma,covariates_habits)
            A=inv_Sigma

            call choldc(A,covariates_habits)
            gamma(:,h_l,e_l)=matmul(inv_Sigma,matmul(transpose(big_X(1:counter_big_X(e_l,h_l),e_l,h_l,:)),big_Y(1:counter_big_X(e_l,h_l),e_l,h_l)))+matmul(A,z(:,1))
            if (isnan(sum(gamma(:,h_l,e_l)))) then
                print*,'pr gamma y'
                print*,counter_big_X(e_l,h_l),sigma
                pause
            end if
        end if
    end do; end do
    
end subroutine
    
subroutine sample_gamma_y_med(gamma_med,type_i,sample_k)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates_habits_med,habits_med,types),intent(inout)::gamma_med
    integer,dimension(indv,1),intent(in)::type_i
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(indv,habits_med,generations)::y_star
    real(DP),dimension(covariates_habits_med,1)::x
    integer::h_l,c_l,g_l,e_d,ge_l,age,it,i_l,e_l
    real(DP)::health_d,ins_d
    real(DP),dimension(indv*g_max,types,habits_med,covariates_habits_med)::big_X
    real(DP),dimension(indv*g_max,types,habits_med)::big_Y
    integer,dimension(types,habits)::counter_big_X
    real(DP),dimension(covariates_habits_med,1)::z
    real(DP),dimension(covariates_habits_med,covariates_habits_med)::Sigma,inv_Sigma,A,B_0
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    y_star=-9.0_sp
    counter_big_X=0
    do i_l=1,indv_HRS;
        if (sample_selection(i_l)) then
            do g_l=first_age(i_l),last_age(i_l);do h_l=1,habits_med
                if (sample_k(i_l,g_l)/=-1 .and. data_ins_hrs(i_l,g_l)/=-1) then
                    health_d=dble(sample_k(i_l,g_l)-1)
                    ins_d=dble(data_ins_hrs(i_l,g_l)) !data_ins_hrs(i_l,:)
                    age=initial_age+(g_l-1)*2
                    x(:,1)=(/1.0_dp,dble(age),dble(age**2.0_dp-1.0_dp),health_d,ins_d/)
                    !x(:,1)=(/1.0_dp,dble(age),dble(age**2.0_dp-1.0_dp),health_d/)
                    if ((data_habits(i_l,habits_med_vec(h_l),g_l)==1 .or. data_habits(i_l,habits_med_vec(h_l),g_l)==0) ) then
                        counter_big_X(type_i(i_l,1),h_l)=counter_big_X(type_i(i_l,1),h_l)+1
                        big_X(counter_big_X(type_i(i_l,1),h_l),type_i(i_l,1),h_l,:)=x(:,1)
                    end if
                    if (data_habits(i_l,habits_med_vec(h_l),g_l)==1 ) then
                        call TRUNCATED_NORMAL_A_SAMPLE(sum(x(:,1)*gamma_med(:,h_l,type_i(i_l,1))),1.0_dp,0.0_dp,y_star(i_l,h_l,g_l))
                        big_Y(counter_big_X(type_i(i_l,1),h_l),type_i(i_l,1),h_l)=y_star(i_l,h_l,g_l)
                    elseif (data_habits(i_l,habits_med_vec(h_l),g_l)==0 ) then
                        call TRUNCATED_NORMAL_B_SAMPLE(sum(x(:,1)*gamma_med(:,h_l,type_i(i_l,1))),1.0_dp,0.0_dp,y_star(i_l,h_l,g_l))
                        big_Y(counter_big_X(type_i(i_l,1),h_l),type_i(i_l,1),h_l)=y_star(i_l,h_l,g_l)
                    end if
                end if
            end do; end do
        end if
    end do
    
    
    do e_l=1,types; do h_l=1,habits_med
        if (counter_big_X(e_l,h_l)>1)then
            B_0=0.0d0
            do c_l=1,covariates_habits
                z(c_l,1)=c4_normal_01(  )
                B_0(c_l,c_l)=0.0d0
            end do
            Sigma=B_0+matmul(transpose(big_X(1:counter_big_X(e_l,h_l),e_l,h_l,:)),big_X(1:counter_big_X(e_l,h_l),e_l,h_l,:)) !big_X(1:counter_big_X(e_l,h_l),e_l,h_l,5)
            call inverse(Sigma,inv_Sigma,covariates_habits_med)
            A=inv_Sigma
            call choldc(A,covariates_habits_med)
            gamma_med(:,h_l,e_l)=matmul(inv_Sigma,matmul(transpose(big_X(1:counter_big_X(e_l,h_l),e_l,h_l,:)),big_Y(1:counter_big_X(e_l,h_l),e_l,h_l)))+matmul(A,z(:,1))
            if (isnan(sum(gamma_med(:,h_l,e_l)))) then
                print*,'pr gamma y med'
                print*,counter_big_X(e_l,h_l),sigma
                pause
            end if
        end if
    end do; end do
    
end subroutine    