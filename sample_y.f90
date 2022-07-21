subroutine sample_y(gamma,y,sample_k,H,weights,type_pr)
    use nrtype; use global_var
    implicit none
    integer,dimension(indv,1),intent(inout)::y
    integer,dimension(indv,1)::y_new
    real(DP),dimension(covariates_habits,habits,types),intent(in)::gamma
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
    real(DP),dimension(covariates_habits,1)::x
    integer::h_l,c_l,g_l,e_d,age,ge_d,it,i_l,e_l,ind,changes
    real(dp)::health_d
    real(dp)::d,p,u,log_likeli
    real(DP),dimension(habits,generations,types,clusters)::alphas
    real(DP),dimension(types)::pr,filtered_pr,selection
    real(DP),dimension(generations,clusters,L_gender,L_educ,types),intent(in)::weights
    real(DP),dimension(indv,types),intent(out)::type_pr

    
    do e_l=1,types; do h_l=1,habits;do c_l=1,clusters; do g_l=1,generations
            age=initial_age+(g_l-1)*2-70
            health_d=dble(c_l-1)
            x(:,1)=(/1.0_dp,dble(age),dble(age**2.0_dp-1.0_dp),health_d/)
            alphas(h_l,g_l,e_l,c_l)=1.0_dp-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*gamma(:,h_l,e_l))/sqrt(2.0_dp)))
    end do;end do; end do;end do

    changes=0
    log_likeli=0
    
    do i_l=1,indv;
        if (race(i_l)==1) then
            pr=1.0d0
            if (sample_k(i_l,first_age(i_l))/=-1) then
                filtered_pr=weights(first_age(i_l),sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l),:) 
            else
                !Change this in the future (just one observation either way)
                filtered_pr=weights(first_age(i_l),1,gender(i_l),educ(i_l),:)
            end if
                
            do g_l=first_age(i_l),last_age(i_l)-1
                do h_l=1,habits; do e_l=1,types
                    if (data_habits(i_l,h_l,g_l)==1 .and. sample_k(i_l,g_l)/=-1 ) then 
                        pr(e_l)=pr(e_l)*alphas(h_l,g_l,e_l,sample_k(i_l,g_l))
                    elseif (data_habits(i_l,h_l,g_l)==0 .and. sample_k(i_l,g_l)/=-1) then 
                        pr(e_l)=pr(e_l)*(1.0d0-alphas(h_l,g_l,e_l,sample_k(i_l,g_l)))
                    end if
                end do; end do
                do e_l=1,types
                    if (sample_k(i_l,g_l)>=1 .and. sample_k(i_l,g_l+1)>=1) then 
                        !if (i_l<=indv_HRS) then
                            filtered_pr(e_l)=filtered_pr(e_l)*H(sample_k(i_l,g_l),sample_k(i_l,g_l+1),g_l,e_l,gender(i_l),educ(i_l))
                        !else
                        !    filtered_pr(e_l)=filtered_pr(e_l)*H(sample_k(i_l,g_l),sample_k(i_l,g_l+1),g_l,e_l,gender(i_l),educ(i_l))/(1.0d0-H(sample_k(i_l,g_l),clusters+1,g_l,e_l,gender(i_l),educ(i_l))) 
                        !end if
                    end if
                if (isnan(sum(filtered_pr))) then
                    print*,'pb sample_y'
                end if
                end do
            end do
            

            pr=pr*filtered_pr/sum(pr*filtered_pr)
            type_pr(i_l,:)=pr
            log_likeli=log_likeli+log(pr(y(i_l,1)))+log(filtered_pr(y(i_l,1)))
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
                if (ind>3) then
                    print*,'here'
                end if
            end do
        end if
    end do
    y=y_new
    
    print*,'likelihood',log_likeli
    
    !print*,'change cluster', real(changes)/real(indv), 'share',real(share)/real(indv)
    
end subroutine