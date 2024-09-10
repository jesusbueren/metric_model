subroutine sim_lives(joint_yh,H)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts),intent(in)::joint_yh
    real(DP),dimension(types)::pr_y
    real(DP),dimension(clusters)::pr_h
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
    real(DP),dimension(L_educ)::LE,SD_L
    integer,parameter::indv_sim=10000
    integer,dimension(indv_sim)::lives_i
    integer::i_l,e_l,c_l,ge_l,y_l,h_l,h_l2,ind,g_l,ini_g
    real(DP)::u
    
    ge_l=1 !Males/Females
    c_l=4 !cohort
    ini_g=12 !initial age for simulation
    
    do e_l=1,L_educ
        do i_l=1,indv_sim
            y_l=-9
            h_l=-9
            h_l2=-9
    
            !Sample initial type
            pr_y=sum(joint_yh(ini_g,:,ge_l,e_l,:,c_l),1)
            call RANDOM_NUMBER(u)
            ind=1
            do while (y_l==-9)
                if (u<sum(pr_y(1:ind))) then
                    y_l=ind
                else
                    ind=ind+1
                end if
            end do
            if (ind>=types+1) then
                print*,'error sampling type: sim_lives'
                pause
            end if
    
            !Sample initial health
            pr_h=joint_yh(ini_g,:,ge_l,e_l,y_l,c_l)/sum(joint_yh(ini_g,:,ge_l,e_l,y_l,c_l))
            call RANDOM_NUMBER(u)
            ind=1
            do while (h_l==-9)
                if (u<sum(pr_h(1:ind))) then
                    h_l=ind
                else
                    ind=ind+1
                end if
            end do
            if (ind>=clusters+1) then
                print*,'error sampling initial health: sim_lives'
                pause
            end if
    
            h_l2=-9
    
            do g_l=ini_g,generations-4
                !Sample future health
                ind=1
                call RANDOM_NUMBER(u)
                do while (h_l2==-9)
                    if (u<sum(H(h_l,1:ind,g_l,y_l,ge_l,e_l))) then
                        h_l2=ind
                    else
                        ind=ind+1
                    end if
                end do
                if (h_l2==clusters+1) EXIT
                h_l=h_l2
                h_l2=-9
            end do
            lives_i(i_l)=(g_l-ini_g)*2
        end do
        LE(e_l)=dble(sum(lives_i))/dble(indv_sim)
        SD_L(e_l)=sum((dble(lives_i)-LE(e_l))**2.0d0)/dble(indv_sim)
        print*,'Av/sd L:,',e_l,LE(e_l),sqrt(SD_L(e_l))
    end do
        
            
        
    
    
    
    
    
        
    
    
    
    
end subroutine