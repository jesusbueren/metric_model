subroutine simulate_data()
    use global_var; use nrtype
    implicit none
    real(DP),dimension(clusters,L_gender,L_educ)::share_h
    real(DP),dimension(clusters,L_gender,L_educ,types)::delta
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H
    real(DP),dimension(generations,clusters,L_gender,L_educ,types)::weights,joint_yh
    real(DP),dimension(habits,generations,types,clusters)::alphas
    integer::i_l,h_l,t_l,h_l2
    integer,dimension(indv,1)::y
    real(DP)::u
    integer,dimension(indv,generations)::sample_k
    
    !Code for simulating data with two types and one education and one gender
    print*,'same gender & educ for all the sample'
    educ=1
    gender=1
    
    
    !Set share of people in gh and bh (could be a parameter to estimate too)
    share_h(1,1,1)=0.8d0
    share_h(2,1,1)=1.d0-share_h(1,1,1)
    
    !Set parameter of weights of behavior types in the initial age given intitial health
    delta(1,1,1,1)=0.8d0
    delta(1,1,1,2)=1.0d0-delta(1,1,1,1)
    
    delta(2,1,1,1)=0.2d0
    delta(2,1,1,2)=1.0d0-delta(2,1,1,1)
    
    !Set transition probabilities
    H(clusters+1,clusters+1,:,:,1,1)=1.0d0
    !for the good type in gh
    H(1,1,:,1,1,1)=0.8d0
    H(1,2,:,1,1,1)=0.15d0
    H(1,3,:,1,1,1)=0.05d0
    !for the good type in bh
    H(2,1,:,1,1,1)=0.4d0
    H(2,2,:,1,1,1)=0.5d0
    H(2,3,:,1,1,1)=0.1d0
    !for the bad type in gh
    H(1,1,:,2,1,1)=0.7d0
    H(1,2,:,2,1,1)=0.2d0
    H(1,3,:,2,1,1)=0.1d0
    !for the bad type in bh
    H(2,1,:,2,1,1)=0.35d0
    H(2,2,:,2,1,1)=0.5d0
    H(2,3,:,2,1,1)=0.15d0
    
    !Set parameters of health behaviors
    alphas(:,:,1,:)=0.8d0
    alphas(:,:,2,:)=0.4d0
    
    call compute_weights(delta,H,share_h,weights,joint_yh) 
    
    print*,'initial share of type I:',sum(joint_yh(1,:,:,:,:),2)

    
    data_habits=-1
    data_shlt=-1
    do i_l=1,indv
        !sample initial health state
        call RANDOM_NUMBER(u)
        if (u<sum(joint_yh(first_age(i_l),1,:,:,:))) then
            data_shlt(i_l,first_age(i_l))=1
        else
            data_shlt(i_l,first_age(i_l))=2
        end if
        !sample health behavior type
        call RANDOM_NUMBER(u)
        if (u<weights(first_age(i_l),data_shlt(i_l,first_age(i_l)),1,1,1)) then 
            y(i_l,1)=1
        else
            y(i_l,1)=2
        end if
        
        t_l=first_age(i_l)
        !sample healh habits given behavior type
1        do h_l=1,habits
            call RANDOM_NUMBER(u)
            if (u<alphas(h_l,t_l,y(i_l,1),data_shlt(i_l,t_l))) then
                data_habits(i_l,h_l,t_l)=1
            else
                data_habits(i_l,h_l,t_l)=0
            end if
        end do
        !Sample health transitions
        if (t_l<generations) then
            call RANDOM_NUMBER(u)
            h_l2=1
            do while (data_shlt(i_l,t_l+1)==-1 .and. h_l2<clusters+2)
                if (u<sum(H(data_shlt(i_l,t_l),1:h_l2,t_l,y(i_l,1),1,1))) then
                    if (i_l<=indv_HRS) then
                        data_shlt(i_l,t_l+1)=h_l2
                        if (h_l2==clusters+1) then
                            last_age(i_l)=t_l+1
                        end if
                    else
                        if (h_l2==clusters+1) then
                            last_age(i_l)=t_l
                        else
                            data_shlt(i_l,t_l+1)=h_l2
                        end if
                    end if                       
                end if
                h_l2=h_l2+1
            end do
            if (data_shlt(i_l,t_l+1)<clusters+1 .and. t_l+1<=last_age(i_l)) then
                t_l=t_l+1
                go to 1
            end if
        end if
    end do
        
            
   y_true=y    
   share_h_true=share_h
   delta_true=delta
        
    
end subroutine
    
