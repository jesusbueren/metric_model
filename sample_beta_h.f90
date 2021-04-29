subroutine sample_beta_h(beta_h,type_i,sample_k)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates,types,clusters,clusters),intent(inout)::beta_h
    integer,dimension(indv,1),intent(in)::type_i
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(covariates,1)::x
    integer::h_l,c_l,g_l,ge_l,age,ge_d,it,i_l,health_d,d_l,h_l2,t_l
    real(DP),dimension(clusters)::h_star1
    double precision,dimension(L_educ-1)::educ_d
    double precision,dimension(types-1)::type_d
    real(DP)::gender_d
    real(DP),dimension(covariates,1)::z
    real(DP),dimension(covariates,covariates)::Sigma,inv_Sigma,A
    real(DP),dimension(indv*g_max,clusters,types,covariates)::big_X_h
    real(DP),dimension(indv*g_max,clusters,types,clusters)::big_Y_h
    integer,dimension(clusters,types)::counter_big_X_h
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
        x(1:4,1)=(/1.0_dp,dble(age),gender_d,dble(age)*gender_d/)
        educ_d=0.0d0
        if (educ(i_l)>1) then
            educ_d(educ(i_l)-1)=1.0d0
        end if
        x(5:6,1)=educ_d
        x(7:8,1)=educ_d*dble(age)

        if (sample_k(i_l,g_l)>=1 .and. sample_k(i_l,g_l+1)>=1 .and. sample_k(i_l,g_l+1)<clusters+1) then
            counter_big_X_h(sample_k(i_l,g_l),type_i(i_l,1))=counter_big_X_h(sample_k(i_l,g_l),type_i(i_l,1))+1
            big_X_h(counter_big_X_h(sample_k(i_l,g_l),type_i(i_l,1)),sample_k(i_l,g_l),type_i(i_l,1),:)=x(:,1)
            !Sample latent h
1           do h_l=1,clusters
                h_star1(h_l)=c4_normal_01()+sum(x(:,1)*beta_h(:,type_i(i_l,1),sample_k(i_l,g_l),h_l))
            end do
            if (maxloc(h_star1,1)==sample_k(i_l,g_l+1)) then
                big_Y_h(counter_big_X_h(sample_k(i_l,g_l),type_i(i_l,1)),sample_k(i_l,g_l),type_i(i_l,1),:)=h_star1
            else
                go to 1
            end if
        end if
    end do;end do
    
    beta_h=0.0_dp
    do h_l=1,clusters;do h_l2=1,clusters-1;do t_l=1,types
        do c_l=1,covariates
            z(c_l,1)=c4_normal_01(  )
        end do
        Sigma=matmul(transpose(big_X_h(1:counter_big_X_h(h_l,t_l),h_l,t_l,:)),big_X_h(1:counter_big_X_h(h_l,t_l),h_l,t_l,:))
        call inverse(Sigma,inv_Sigma,covariates)
        A=inv_Sigma
        call choldc(A,covariates)
        beta_h(:,t_l,h_l,h_l2)=matmul(inv_Sigma,matmul(transpose(big_X_h(1:counter_big_X_h(h_l,t_l),h_l,t_l,:)),big_Y_h(1:counter_big_X_h(h_l,t_l),h_l,t_l,h_l2)))+matmul(A,z(:,1))
    end do;end do;end do

end subroutine