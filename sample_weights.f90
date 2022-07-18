!subroutine sample_weights(weights_all,y,survival_pr) 
!    use global_var; use nrtype
!    implicit none
!    real(DP),dimension(generations,types,L_gender,L_educ),intent(inout)::weights_all
!    integer,dimension(indv,1),intent(in)::y
!    real(DP),dimension(types,L_gender,L_educ,generations),intent(in)::survival_pr
!    real(DP),dimension(L_gender,L_educ)::pdf_old,pdf_new
!    real(DP),dimension(generations,types,L_gender,L_educ)::weights_all_new
!    integer:: e_l,ge_l
!    real(DP):: u
!    
!    call compute_weights(weights_all,survival_pr) 
!    
!    call compute_weight_density(weights_all,y,pdf_old)
!    
!    call proposal_weights(weights_all,survival_pr,weights_all_new)
!    
!    call compute_weight_density(weights_all_new,y,pdf_new)
!    
!    do ge_l=1,L_gender;do e_l=1,L_educ
!        if (pdf_new(ge_l,e_l)>pdf_old(ge_l,e_l)) then
!            weights_all(:,:,ge_l,e_l)=weights_all_new(:,:,ge_l,e_l)
!        else
!            call RANDOM_NUMBER(u)
!            if (u<exp(pdf_new(ge_l,e_l)-pdf_old(ge_l,e_l))) then
!                weights_all(:,:,ge_l,e_l)=weights_all_new(:,:,ge_l,e_l)
!            end if
!        end if
!    end do; end do
!    
!end subroutine
!
!subroutine compute_weights(weights_all,survival_pr)
!    use global_var; use nrtype
!    implicit none
!    real(DP),dimension(generations,types,L_gender,L_educ),intent(inout)::weights_all
!    real(DP),dimension(types,L_gender,L_educ,generations),intent(in)::survival_pr
!    integer:: e_l,g_l,ge_l
!    
!    do ge_l=1,L_gender;do e_l=1,L_educ; do g_l=2,generations  
!        weights_all(g_l,:,ge_l,e_l)=weights_all(g_l-1,:,ge_l,e_l)*survival_pr(:,ge_l,e_l,g_l)/sum(weights_all(g_l-1,:,ge_l,e_l)*survival_pr(:,ge_l,e_l,g_l))
!    end do; end do; end do
!    
!end subroutine
!    
!    
!subroutine compute_weight_density(weights_all,y,pdf)   
!    use global_var; use nrtype
!    implicit none
!    real(DP),dimension(generations,types,L_gender,L_educ),intent(in)::weights_all
!    integer,dimension(indv,1),intent(in)::y
!    real(DP),dimension(L_gender,L_educ),intent(out)::pdf
!    integer::i_l
!    
!    pdf=0.0d0
!    do i_l=1,indv
!        pdf(gender(i_l),educ(i_l))=pdf(gender(i_l),educ(i_l))+log(weights_all(first_age(i_l),y(i_l,1),gender(i_l),educ(i_l)))
!    end do
!end subroutine
!    
!subroutine proposal_weights(weights_all,survival_pr,weights_all_new)
!use global_var; use nrtype
!implicit none
!real(DP),dimension(generations,types,L_gender,L_educ),intent(in)::weights_all
!real(DP),dimension(types,L_gender,L_educ,generations),intent(in)::survival_pr
!real(DP),dimension(generations,types,L_gender,L_educ),intent(out)::weights_all_new
!real(DP),dimension(types)::u
!integer:: ge_l,e_l
!
!do ge_l=1,L_gender;do e_l=1,L_educ
!    call RANDOM_NUMBER(u(1:types-1))
!    u(1:types-1)=(u(1:types-1)-0.5d0)*0.001d0
!    u(types)=-sum(u(1:types-1))
!    weights_all_new(1,:,ge_l,e_l)=weights_all(1,:,ge_l,e_l)+u
!end do; end do
!
!call compute_weights(weights_all_new,survival_pr)
!
!
!end subroutine