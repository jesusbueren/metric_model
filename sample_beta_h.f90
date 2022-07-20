subroutine sample_beta_h(beta_h,type_i,sample_k)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates,types,clusters,L_gender,L_educ),intent(inout)::beta_h
    integer,dimension(indv,1),intent(in)::type_i
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(covariates,1)::x
    integer::h_l,c_l,g_l,ge_l,age,ge_d,it,i_l,health_d,d_l,t_l,e_l
    real(DP)::h_star1
    real(DP),dimension(covariates,1)::z
    real(DP),dimension(covariates,covariates)::Sigma,inv_Sigma,A,Sigma_aux
    real(DP),dimension(indv*g_max,clusters,types,L_gender,L_educ,covariates)::big_X_h
    real(DP),dimension(indv*g_max,clusters,types,L_gender,L_educ)::big_Y_h
    integer,dimension(clusters,types,L_gender,L_educ)::counter_big_X_h
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    counter_big_X_h=0
    do i_l=1,indv;do g_l=first_age(i_l),last_age(i_l)-1
        x=-9.0d0
        age=initial_age+(g_l-1)*2-70
        x(1:3,1)=(/1.0_dp,dble(age),dble(age)**2.0d0/)
        if (sample_k(i_l,g_l)>=1 .and. sample_k(i_l,g_l+1)>=1 .and. sample_k(i_l,g_l+1)<clusters+1 .and. race(i_l)==1) then
            counter_big_X_h(sample_k(i_l,g_l),type_i(i_l,1),gender(i_l),educ(i_l))=counter_big_X_h(sample_k(i_l,g_l),type_i(i_l,1),gender(i_l),educ(i_l))+1
            big_X_h(counter_big_X_h(sample_k(i_l,g_l),type_i(i_l,1),gender(i_l),educ(i_l)),sample_k(i_l,g_l),type_i(i_l,1),gender(i_l),educ(i_l),:)=x(:,1)
            !Sample latent h
            if (sample_k(i_l,g_l+1)==1) then
                call TRUNCATED_NORMAL_A_SAMPLE(sum(x(:,1)*beta_h(:,type_i(i_l,1),sample_k(i_l,g_l),gender(i_l),educ(i_l))),1.0_dp,0.0_dp,h_star1)
            elseif ( sample_k(i_l,g_l+1)==2) then
                call TRUNCATED_NORMAL_B_SAMPLE(sum(x(:,1)*beta_h(:,type_i(i_l,1),sample_k(i_l,g_l),gender(i_l),educ(i_l))),1.0_dp,0.0_dp,h_star1)
            end if
            if (h_star1>100 .or. h_star1<-100) then
                    print*,i_l,sample_k(i_l,g_l),sample_k(i_l,g_l+1),type_i(i_l,1),g_l,gender(i_l),educ(i_l)
                    pause
            end if
            big_Y_h(counter_big_X_h(sample_k(i_l,g_l),type_i(i_l,1),gender(i_l),educ(i_l)),sample_k(i_l,g_l),type_i(i_l,1),gender(i_l),educ(i_l))=h_star1
        end if
    end do;end do
    
    
    do h_l=1,clusters;do t_l=1,types;do e_l=1,L_educ;do ge_l=1,L_gender
        if (counter_big_X_h(h_l,t_l,ge_l,e_l)>5) then
            beta_h(:,t_l,h_l,ge_l,e_l)=0.0d0
            do c_l=1,covariates
                z(c_l,1)=c4_normal_01(  )
            end do
            Sigma=matmul(transpose(big_X_h(1:counter_big_X_h(h_l,t_l,ge_l,e_l),h_l,t_l,ge_l,e_l,:)),big_X_h(1:counter_big_X_h(h_l,t_l,ge_l,e_l),h_l,t_l,ge_l,e_l,:))
            Sigma_aux=Sigma
            call inverse(Sigma,inv_Sigma,covariates)
            A=inv_Sigma
            call choldc(A,covariates)
            beta_h(:,t_l,h_l,ge_l,e_l)=matmul(inv_Sigma,matmul(transpose(big_X_h(1:counter_big_X_h(h_l,t_l,ge_l,e_l),h_l,t_l,ge_l,e_l,:)),big_Y_h(1:counter_big_X_h(h_l,t_l,ge_l,e_l),h_l,t_l,ge_l,e_l)))+matmul(A,z(:,1))
            if (isnan(sum(beta_h(:,t_l,h_l,ge_l,e_l)))) then
                print*,'error beta_h'
            end if
        else
            print*,'strange beta_h'
        end if
    end do;end do;end do;end do

end subroutine