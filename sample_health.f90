subroutine sample_health(H,joint_yh,y,sample_k)
    use global_var;use nrtype
    implicit none
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts),intent(in)::joint_yh
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
    integer,dimension(indv,generations),intent(out)::sample_k
    integer,dimension(indv,1),intent(in)::y
    integer::i_l, g_l,h_l,h_l2
    real(DP),dimension(clusters+1,generations)::filter_h
    real(DP),dimension(clusters+1)::smooth_pr
    integer,dimension(indv,generations)::data_shlt_aux
    
    data_shlt_aux=data_shlt
!H(1,3,:,1,1,3)
    do i_l=1,indv
        !Filter 
        filter_h=0.0d0
        do g_l=1,last_age(i_l) !first_age(i_l)
            if (obs_shlt(i_l,g_l)==1) then !obs_shlt(1,:)
                filter_h(data_shlt_aux(i_l,g_l),g_l)=1.0d0
            else
                if (g_l==1) then
                    do h_l=1,clusters
                        filter_h(h_l,g_l)=joint_yh(g_l,h_l,gender(i_l),educ(i_l),y(i_l,1),birth_cohort(i_l))/sum(joint_yh(g_l,:,gender(i_l),educ(i_l),y(i_l,1),birth_cohort(i_l)))
                    end do
                else
                    do h_l=1,clusters
                        do h_l2=1,clusters
                            filter_h(h_l,g_l)=filter_h(h_l,g_l)+H(h_l2,h_l,g_l,y(i_l,1),gender(i_l),educ(i_l))/sum(H(h_l2,1:clusters,g_l,y(i_l,1),gender(i_l),educ(i_l)))*filter_h(h_l2,g_l-1) !filter_h(:,g_l-1)
                        end do
                    end do
                end if
                !if (isnan(sum(filter_h))) then !
                !    print*,''
                !end if
            end if
        end do  
        !Sample
        do g_l=last_age(i_l),1,-1
            if (g_l==last_age(i_l)) then
                filter_h(:,g_l)=filter_h(:,g_l)/sum(filter_h(:,g_l))
                if (isnan(sum(filter_h(:,g_l)))) then
                    filter_h(:,g_l)=0.5d0
                end if
                call sample_h_i(filter_h(:,g_l),data_shlt_aux(i_l,g_l))
            else
                smooth_pr=0.0d0
                do h_l=1,clusters
                    smooth_pr(h_l)=H(h_l,data_shlt_aux(i_l,g_l+1),g_l,y(i_l,1),gender(i_l),educ(i_l))*filter_h(h_l,g_l)
                end do
                smooth_pr=smooth_pr/sum(smooth_pr)
                if (isnan(sum(smooth_pr))) then
                    smooth_pr=0.5d0
                end if
                call sample_h_i(smooth_pr,data_shlt_aux(i_l,g_l))
            end if
            if (obs_shlt(i_l,g_l)==1 .and. data_shlt_aux(i_l,g_l)/=data_shlt(i_l,g_l)) then
                print*,''
            end if
        end do
    end do       

    sample_k=data_shlt_aux
    end subroutine
    
    subroutine sample_h_i(pr_h,h_out)
    use nrtype; use global_var
    implicit none
    real(DP),dimension(clusters+1)::pr_h
    integer,intent(out):: h_out
    real(DP)::u_h
    integer::ind_h
    
    
    h_out=-9
    ind_h=1
    call RANDOM_NUMBER(u_h)
    do while (h_out==-9)
        if (u_h<sum(pr_h(1:ind_h))) then
            h_out=ind_h
        else
            ind_h=ind_h+1
        end if
    end do
                
    end subroutine
    