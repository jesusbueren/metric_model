subroutine sample_delta(delta,H,share_h,y,sample_k,weights,joint_yh) 
    use global_var; use nrtype
    implicit none
    real(DP),dimension(clusters,L_gender,L_educ,types),intent(inout)::delta
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
    real(DP),dimension(clusters,L_gender,L_educ),intent(in)::share_h
    integer,dimension(indv,1),intent(in)::y
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(generations,clusters,L_gender,L_educ,types),intent(out)::weights,joint_yh
    real(DP),dimension(covariates_mixture,L_gender,L_educ,types)::delta_g
    integer::i_l,e_l,ge_l,y_l,h_l
    real(DP),dimension(L_gender,L_educ)::log_likeli,log_likeli_g
    real(DP),dimension(covariates_mixture,types)::u
    real(DP),dimension(generations,clusters,L_gender,L_educ,types)::weights_g,joint_yh_g
    real(DP)::eps
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    call compute_weights(delta,H,share_h,weights,joint_yh) 
    
    !Compute likelihood of weights
    log_likeli=0.0d0
    do i_l=1,indv
        if (sample_k(i_l,first_age(i_l))>=1 .and. sample_k(i_l,first_age(i_l))<=2) then
            log_likeli(gender(i_l),educ(i_l))=log_likeli(gender(i_l),educ(i_l))+log(weights(first_age(i_l),sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l),y(i_l,1)))
        else
            log_likeli(gender(i_l),educ(i_l))=log_likeli(gender(i_l),educ(i_l))+log(weights(first_age(i_l),1,gender(i_l),educ(i_l),y(i_l,1)))
        end if
    end do
            
    !Generate proposal
    do e_l=1,L_educ;do ge_l=1,L_gender
        u=0.0d0
        do h_l=1,clusters;
            do y_l=1,types
                if (y_l<types) then
                    u(h_l,y_l)=c4_normal_01()/100.0d0        
                else
                    u(h_l,y_l)=-sum(u(h_l,:))
                end if
                delta_g(h_l,ge_l,e_l,y_l)=max(delta(h_l,ge_l,e_l,y_l)+u(h_l,y_l),0.01d0)
            end do
            delta_g(h_l,ge_l,e_l,:)=delta_g(h_l,ge_l,e_l,:)/sum(delta_g(h_l,ge_l,e_l,:))
        end do
        
        
    end do; end do
    
    call compute_weights(delta_g,H,share_h,weights_g,joint_yh_g) 
    
    !Compute likelihood of proposal
    log_likeli_g=0.0d0
    do i_l=1,indv
        if (sample_k(i_l,first_age(i_l))>=1 .and. sample_k(i_l,first_age(i_l))<=2) then
            log_likeli_g(gender(i_l),educ(i_l))=log_likeli_g(gender(i_l),educ(i_l))+log(weights_g(first_age(i_l),sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l),y(i_l,1)))
        else
            log_likeli_g(gender(i_l),educ(i_l))=log_likeli_g(gender(i_l),educ(i_l))+log(weights_g(first_age(i_l),1,gender(i_l),educ(i_l),y(i_l,1)))
        end if
    end do
    
    !Accept/reject proposal using Metropolis algorithm
    do e_l=1,L_educ;do ge_l=1,L_gender
        if(log_likeli_g(ge_l,e_l)>log_likeli(ge_l,e_l)) then
            delta(:,ge_l,e_l,:)=delta_g(:,ge_l,e_l,:)
            weights(:,:,ge_l,e_l,:)=weights_g(:,:,ge_l,e_l,:)
            joint_yh(:,:,ge_l,e_l,:)=joint_yh_g(:,:,ge_l,e_l,:)
        else
           call RANDOM_NUMBER(eps)
           if (eps<exp(log_likeli_g(ge_l,e_l)-log_likeli(ge_l,e_l))) then
               delta(:,ge_l,e_l,:)=delta_g(:,ge_l,e_l,:)
                weights(:,:,ge_l,e_l,:)=weights_g(:,:,ge_l,e_l,:)
                joint_yh(:,:,ge_l,e_l,:)=joint_yh_g(:,:,ge_l,e_l,:)
           end if
        end if            
    end do; end do
    
    
    

end subroutine