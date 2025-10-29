subroutine sample_y(gamma,gamma_med,y,sample_k,H,weights,type_pr)
    use nrtype; use global_var
    implicit none
    integer,dimension(indv,1),intent(inout)::y
    integer,dimension(indv,1)::y_new
    real(DP),dimension(covariates_habits,habits_nomed,types),intent(in)::gamma
    real(DP),dimension(covariates_habits_med,habits_med,types),intent(in)::gamma_med
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
    real(DP),dimension(covariates_habits,1)::x
    real(DP),dimension(covariates_habits_med,1)::x_med
    integer::h_l,c_l,g_l,e_d,age,ge_d,it,i_l,e_l,ind,changes,ins_l
    real(dp)::health_d,ins_d,maxlog
    real(dp)::d,p,u,log_likeli
    real(DP),dimension(habits_nomed,generations,types,clusters)::alphas
    real(DP),dimension(habits_med,generations,types,clusters,2)::alphas_med
    real(DP),dimension(types)::pr,filtered_pr,selection
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts),intent(in)::weights
    real(DP),dimension(indv,types),intent(out)::type_pr

    
    do e_l=1,types; do h_l=1,habits_nomed;do c_l=1,clusters; do g_l=1,generations
            age=initial_age+(g_l-1)*2
            health_d=dble(c_l-1)
            x(:,1)=(/1.0_dp,dble(age),dble(age**2.0_dp-1.0_dp),health_d/)
            alphas(h_l,g_l,e_l,c_l)=1.0-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*gamma(:,h_l,e_l))/sqrt(2.0_dp))) !alphas(2,:,2,1)
    end do;end do; end do;end do
    do e_l=1,types; do h_l=1,habits_med;do c_l=1,clusters; do g_l=1,generations; do ins_l=1,2
            age=initial_age+(g_l-1)*2
            health_d=dble(c_l-1)
            ins_d=dble(ins_l-1)
            x_med(:,1)=(/1.0_dp,dble(age),dble(age**2.0_dp-1.0_dp),health_d,ins_d/)
            !x_med(:,1)=(/1.0_dp,dble(age),dble(age**2.0_dp-1.0_dp),health_d/)
            alphas_med(h_l,g_l,e_l,c_l,ins_l)=1.0d0-0.5_dp*(1.0_dp+erf(-sum(x_med(:,1)*gamma_med(:,h_l,e_l))/sqrt(2.0_dp))) !alphas_med(1,1,:,1,2)
    end do;end do; end do;end do;end do

    changes=0
    log_likeli=0
    
    do i_l=1,indv;
        if (sample_selection(i_l)) then
            pr=0.0d0
            if (sample_k(i_l,first_age(i_l))/=-1) then
                filtered_pr=log(weights(first_age(i_l),sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l),:,birth_cohort(i_l)))
            else
                !just one observation either way
                filtered_pr=log(weights(first_age(i_l),1,gender(i_l),educ(i_l),:,birth_cohort(i_l)))
            end if
                
            do g_l=first_age(i_l),last_age(i_l)-1
                do h_l=1,habits_nomed; do e_l=1,types
                    if (data_habits(i_l,habits_vec(h_l),g_l)==1 .and. sample_k(i_l,g_l)/=-1 ) then !data_habits(i_l,habits_vec(2),:)
                        pr(e_l)=pr(e_l)+log(alphas(h_l,g_l,e_l,sample_k(i_l,g_l)))
                    elseif (data_habits(i_l,habits_vec(h_l),g_l)==0 .and. sample_k(i_l,g_l)/=-1) then  
                        pr(e_l)=pr(e_l)+log(1.0d0-alphas(h_l,g_l,e_l,sample_k(i_l,g_l)))
                    end if
                end do; end do !data_habits(i_l,habits_med_vec(3),:)
                if (i_l<=indv_HRS) then
                    do h_l=1,3; do e_l=1,types
                        if (data_habits(i_l,habits_med_vec(h_l),g_l)==1 .and. sample_k(i_l,g_l)/=-1 .and. data_ins_hrs(i_l,g_l)/=-1) then  !alphas_med(h_l,g_l,:,sample_k(i_l,g_l),data_ins_hrs(i_l,g_l)+1)
                            pr(e_l)=pr(e_l)+log(alphas_med(h_l,g_l,e_l,sample_k(i_l,g_l),data_ins_hrs(i_l,g_l)+1))
                        elseif (data_habits(i_l,habits_med_vec(h_l),g_l)==0 .and. sample_k(i_l,g_l)/=-1 .and. data_ins_hrs(i_l,g_l)/=-1) then 
                            pr(e_l)=pr(e_l)+log(1.0d0-alphas_med(h_l,g_l,e_l,sample_k(i_l,g_l),data_ins_hrs(i_l,g_l)+1))
                        end if
                    end do; end do
                end if
                do e_l=1,types
                    if (sample_k(i_l,g_l)>=1 .and. sample_k(i_l,g_l+1)>=1) then !sample_k(i_l,31) data_habits(i_l,:,29:31) sample_k(i_l,:) H(1,3,:,e_l,gender(i_l),educ(i_l))
                        if (i_l<=indv_HRS) then
                            filtered_pr(e_l)=filtered_pr(e_l)+log(H(sample_k(i_l,g_l),sample_k(i_l,g_l+1),g_l,e_l,gender(i_l),educ(i_l))) !H(sample_k(i_l,g_l),sample_k(i_l,g_l+1),g_l,e_l,gender(i_l),educ(i_l))
                        else
                            !filtered_pr(e_l)=filtered_pr(e_l)*min(H(sample_k(i_l,g_l),sample_k(i_l,g_l+1),g_l,e_l,gender(i_l),educ(i_l)), 1.0d-8)/(1.0d0-H(sample_k(i_l,g_l),clusters+1,g_l,e_l,gender(i_l),educ(i_l))) 
                            if (H(sample_k(i_l,g_l),clusters+1,g_l,e_l,gender(i_l),educ(i_l))==1.0d0) then
                                filtered_pr(e_l)=filtered_pr(e_l)+log(1.0d-15)
                                !print*,'got here caution in sample y'
                            else
                                filtered_pr(e_l)=filtered_pr(e_l)+log(H(sample_k(i_l,g_l),sample_k(i_l,g_l+1),g_l,e_l,gender(i_l),educ(i_l))/(1.0d0-H(sample_k(i_l,g_l),clusters+1,g_l,e_l,gender(i_l),educ(i_l))))
                            end if
                        end if
                    end if
                    if (isnan(sum(filtered_pr))) then
                        print*,'pb sample_y nan'
                    elseif (sum(filtered_pr)==0.0d0) then
                        !print*,'pb sample_y zero'
                    end if
                end do
            end do

            
            log_likeli=log_likeli+pr(y(i_l,1))
            pr=pr+filtered_pr
            maxlog = maxval(pr)
            pr = exp(pr - maxlog)
            pr=pr/sum(pr)
            type_pr(i_l,:)=pr
            log_likeli=log_likeli+filtered_pr(y(i_l,1))
            if (isnan(log_likeli))  then
                print*,'problem sample y'
            end if
                
            y_new(i_l,1)=-9
            call RANDOM_NUMBER(u)
            ind=1
            do while (y_new(i_l,1)==-9)
                if (u<sum(pr(1:ind),1) .or. ind==types) then
                    y_new(i_l,1)=ind
                    if (y_new(i_l,1)/=y(i_l,1)) then
                        changes=changes+1
                    end if
                else
                    ind=ind+1
                end if 
                if (ind>types) then
                    print*,'here'
                end if
            end do
        end if
    end do
    y=y_new !y_new(1:1000,1)
    
end subroutine