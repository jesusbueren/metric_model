subroutine sample_y(gamma,y,fraction_t)
    use nrtype; use global_var
    implicit none
    integer,dimension(indv,1),intent(out)::y
    real(DP),dimension(covariates_habits,habits,types),intent(in)::gamma
    real(DP),dimension(covariates_habits,1)::x
    integer::h_l,c_l,g_l,e_d,age,ge_d,it,i_l,e_l,health_d,ind
    real(dp)::d,p,u
    real(DP),dimension(habits,generations,types)::alphas
    real(DP),dimension(types)::pr
    real(DP),dimension(types,L_gender,L_educ),intent(out)::fraction_t
    real(DP),dimension(L_gender,L_educ)::counter
    
    do e_l=1,types
        do h_l=1,habits;do c_l=1,clusters; do g_l=1,generations
            age=initial_age+(g_l-1)*2-70
            health_d=c_l-1
            x(:,1)=(/1.0_dp,dble(age),dble(age**2.0_dp-1.0_dp)/)
            alphas(h_l,g_l,e_l)=1.0_dp-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*gamma(:,h_l,e_l))/sqrt(2.0_dp)))
        end do;end do; end do
    end do
    
    counter=0.0d0
    fraction_t=0.0d0
    do i_l=1,indv;
        pr=1.0_dp
        do g_l=first_age(i_l),last_age(i_l);do h_l=1,habits
            do e_l=1,types
                if (data_habits(i_l,h_l,g_l)==1) then
                    pr(e_l)=pr(e_l)*alphas(h_l,g_l,e_l)
                end if
            end do
        end do; end do
        pr=pr/sum(pr)
        counter(gender(i_l),educ(i_l))=counter(gender(i_l),educ(i_l))+1.0d0
        fraction_t(:,gender(i_l),educ(i_l))=(counter(gender(i_l),educ(i_l))-1.0d0)/counter(gender(i_l),educ(i_l))*fraction_t(:,gender(i_l),educ(i_l)) + &
                                             1.0d0/counter(gender(i_l),educ(i_l))*pr
        y(i_l,1)=-9
        call RANDOM_NUMBER(u)
        ind=1
        do while (y(i_l,1)==-9)
            if (u<sum(pr(1:ind),1) .or. ind==types) then
                y(i_l,1)=ind
            else
                ind=ind+1
            end if 
        end do
    end do
    
end subroutine