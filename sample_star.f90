subroutine sample_y_star(gamma,type_i,sample_k,y_star)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates_habits,habits,types),intent(in)::gamma
    integer,dimension(indv,1),intent(in)::type_i
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(indv,habits,generations),intent(out)::y_star
    real(DP),dimension(covariates_habits,1)::x
    integer::h_l,c_l,g_l,e_d,ge_l,age,it,i_l,health_d
    
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
    
end subroutine
    
subroutine sample_h_star(beta,type_i,sample_k)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates,types,clusters,clusters+1),intent(in)::beta
    integer,dimension(indv,1),intent(in)::type_i
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(covariates,1)::x
    integer::h_l,c_l,g_l,ge_l,age,ge_d,it,i_l,health_d,d_l
    real(DP),dimension(clusters+1)::h_star1
    double precision,dimension(L_educ-1)::educ_d
    double precision,dimension(types-1)::type_d
    real(DP)::gender_d
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    counter_big_X_h=0
    do i_l=1,indv;do g_l=first_age(i_l),last_age(i_l)-1
        x=-9.0d0
        age=initial_age+(g_l-1)*2-70
        gender_d=dble(gender(i_l)-1)
        x(1:5,1)=(/1.0_dp,dble(age),dble(age**2.0_dp-1.0_dp),gender_d,dble(age)*gender_d/)
        educ_d=0.0d0
        if (educ(i_l)>1) then
            educ_d(educ(i_l)-1)=1.0d0
        end if
        x(6:7,1)=educ_d
        x(8:9,1)=educ_d*dble(age)

        if (sample_k(i_l,g_l)>=1 .and. sample_k(i_l,g_l+1)>=1) then
            counter_big_X_h(sample_k(i_l,g_l),type_i(i_l,1))=counter_big_X_h(sample_k(i_l,g_l),type_i(i_l,1))+1
            big_X_h(counter_big_X_h(sample_k(i_l,g_l),type_i(i_l,1)),sample_k(i_l,g_l),type_i(i_l,1),:)=x(:,1)
            !Sample latent h
1           do h_l=1,clusters+1
                h_star1(h_l)=c4_normal_01()+sum(x(:,1)*beta(:,type_i(i_l,1),sample_k(i_l,g_l),h_l))
            end do
            if (maxloc(h_star1,1)==sample_k(i_l,g_l+1)) then
                big_Y_h(counter_big_X_h(sample_k(i_l,g_l),type_i(i_l,1)),sample_k(i_l,g_l),type_i(i_l,1),:)=h_star1
            else
                go to 1
            end if
        end if
    end do;end do

end subroutine