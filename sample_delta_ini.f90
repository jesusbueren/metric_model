subroutine sample_delta_ini(delta,type_i,sample_k) 
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates_mixture,L_gender,L_educ,types),intent(inout)::delta
    integer,dimension(indv,1),intent(in)::type_i
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(covariates_mixture,1)::x
    integer::h_l,g_l,ge_l,age,ge_d,it,i_l,d_l,h_l2,t_l,e_l,y_l,c_l
    real(DP),dimension(types)::y_star,epsilon,y_guess
    real(DP),dimension(covariates_mixture,1)::z
    real(DP),dimension(covariates_mixture,covariates_mixture)::Sigma,inv_Sigma,A
    real(DP),dimension(indv,L_gender,L_educ,covariates_mixture)::big_X
    real(DP),dimension(indv,types,L_gender,L_educ)::big_Y
    integer,dimension(L_gender,L_educ)::counter_big_X
    integer::ind
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    counter_big_X=0
    do i_l=1,indv;
        g_l=first_age(i_l)
        x=-9.0d0
        x(1:covariates_mixture,1)=(/1.0_dp,h_bar(i_l)/)
        if (sample_k(i_l,g_l)>=1 .and. race(i_l)==1 .and. first_age(i_l)<5) then
            counter_big_X(gender(i_l),educ(i_l))=counter_big_X(gender(i_l),educ(i_l))+1
            big_X(counter_big_X(gender(i_l),educ(i_l)),gender(i_l),educ(i_l),:)=x(:,1) 
            y_star=-1
            ind=0
            do while (y_star(1)==-1 .and. ind<10000)
                epsilon=-1
                do y_l=1,types
                    epsilon(y_l)=c4_normal_01() 
                    y_guess(y_l)=sum(x(1:covariates_mixture,1)*delta(1:covariates_mixture,gender(i_l),educ(i_l),y_l))+epsilon(y_l)
                end do
                if (maxloc(y_guess,1)==type_i(i_l,1)) then
                    y_star=y_guess
                end if
                ind=ind+1
            end do
            if (ind>10000) then
                print*,'problem in sample delta',i_l,sample_k(i_l,g_l),gender(i_l),educ(i_l)
            end if
            big_Y(counter_big_X(gender(i_l),educ(i_l)),:,gender(i_l),educ(i_l))=y_star
        end if
    end do
    
    delta=0.0_dp
    do y_l=1,types-1;do e_l=1,L_educ;do ge_l=1,L_gender
        do c_l=1,covariates_mixture
            z(c_l,1)=c4_normal_01(  )
        end do
        Sigma=matmul(transpose(big_X(1:counter_big_X(ge_l,e_l),ge_l,e_l,:)),big_X(1:counter_big_X(ge_l,e_l),ge_l,e_l,:))
        call inverse(Sigma,inv_Sigma,covariates_mixture)
        A=inv_Sigma
        call choldc(A,covariates_mixture)
        delta(:,ge_l,e_l,y_l)=matmul(inv_Sigma,matmul(transpose(big_X(1:counter_big_X(ge_l,e_l),ge_l,e_l,:)),big_Y(1:counter_big_X(ge_l,e_l),y_l,ge_l,e_l)))+matmul(A,z(:,1))
    end do;end do;end do
    

end subroutine