!subroutine smoothing(filtered_states,H,y,sample_k,smoothed_states,init_cond)
!    use global_var
!    implicit none
!    integer,dimension(indv,generations),intent(out)::sample_k
!    double precision,dimension(indv,clusters,generations),intent(in)::filtered_states
!    double precision,dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
!    integer,dimension(indv,1),intent(in)::y
!    real(DP),dimension(clusters,types,L_gender,L_educ),intent(out)::init_cond
!    double precision,dimension(clusters)::smoothed_states_cond
!    double precision,dimension(clusters+1,clusters)::smoothed_states_uncond
!    real(DP),dimension(indv,clusters+1,generations)::smoothed_states
!    integer,dimension(types,L_gender,L_educ)::counter
!    integer::g_l,i_l,ind,c_l2,e_l
!    double precision::u
!    
!    !Initialize smoother
!    smoothed_states=0.0d0
!    sample_k=0
!    counter=0
!    init_cond=0.0d0
!    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(g_l,i_l,u,ind,smoothed_states_cond,smoothed_states_uncond)
!    !$OMP DO 
!    do i_l=1,indv;do g_l=last_age(i_l),1,-1;
!        if (data_adls(i_l,1,g_l)==-1 .and. g_l>=first_age(i_l) .and. g_l<=last_age(i_l)) then  !Dead
!            sample_k(i_l,g_l)=clusters+1
!            smoothed_states(i_l,:,g_l)=0.0d0
!            smoothed_states(i_l,clusters+1,g_l)=1.0d0
!        elseif (g_l>last_age(i_l)) then !Not interviewed
!            sample_k(i_l,g_l)=-1
!        elseif (g_l==last_age(i_l))  then !Last interview when alive if always alive 
!            smoothed_states(i_l,:,g_l)=0.0d0
!            smoothed_states(i_l,1:clusters,g_l)=filtered_states(i_l,1:clusters,g_l)
!            call RANDOM_NUMBER(u)
!            ind=1
!            do while (sample_k(i_l,g_l)==0)
!                if (ind==1) then
!                    if (u<filtered_states(i_l,1,g_l)) then
!                        sample_k(i_l,g_l)=ind
!                    end if                
!                else
!                    if (u<sum(filtered_states(i_l,1:ind,g_l))) then 
!                        sample_k(i_l,g_l)=ind
!                    end if
!                end if
!                ind=ind+1
!            end do 
!        elseif (g_l<last_age(i_l) .or. (data_adls(i_l,1,g_l)/=-1 .and. sample_k(i_l,g_l+1)==clusters+1)) then !not last interview and alive or last alive interview
!            if (sample_k(i_l,g_l+1)==0) then
!                print*,'Error'
!            end if
!            !Kim smoother    
!            smoothed_states_cond(:)=(H(1:clusters,sample_k(i_l,g_l+1),g_l,y(i_l,1),gender(i_l),educ(i_l))*filtered_states(i_l,1:clusters,g_l))/sum((H(1:clusters,sample_k(i_l,g_l+1),g_l,y(i_l,1),gender(i_l),educ(i_l))*filtered_states(i_l,1:clusters,g_l))) !revisa
!            call RANDOM_NUMBER(u)
!            ind=1 
!            do while (sample_k(i_l,g_l)==0)
!               if (u<sum(smoothed_states_cond(1:ind))) then
!                    sample_k(i_l,g_l)=ind
!                end if
!                ind=ind+1
!            end do
!            !Hamilton Smoother
!            do c_l2=1,clusters+1
!                if (smoothed_states(i_l,c_l2,g_l+1)==0.0d0)then
!                    smoothed_states_uncond(c_l2,:)=0.0d0
!                else
!                    smoothed_states_uncond(c_l2,:)=smoothed_states(i_l,c_l2,g_l+1)* &
!                                                (H(1:clusters,c_l2,g_l,y(i_l,1),gender(i_l),educ(i_l))*filtered_states(i_l,1:clusters,g_l))/&
!                    sum((H(1:clusters,c_l2,g_l,y(i_l,1),gender(i_l),educ(i_l))*filtered_states(i_l,1:clusters,g_l)))
!                end if
!            end do
!            smoothed_states(i_l,1:clusters,g_l)=sum(smoothed_states_uncond,1)
!        end if
!    end do 
!    counter(y(i_l,1),gender(i_l),educ(i_l))=counter(y(i_l,1),gender(i_l),educ(i_l))+1
!    init_cond(:,y(i_l,1),gender(i_l),educ(i_l))=init_cond(:,y(i_l,1),gender(i_l),educ(i_l))+smoothed_states(i_l,1:clusters,1) !smoothed_states(1,1,:) 
!    end do
!    !$OMP END DO
!    !$OMP END PARALLEL
!    
!    do i_l=1,types; do g_l=1,L_gender;do e_l=1,L_educ
!        init_cond(:,i_l,g_l,e_l)=init_cond(:,i_l,g_l,e_l)/dble(counter(i_l,g_l,e_l))
!    end do;end do;end do
!    
!end subroutine
    
