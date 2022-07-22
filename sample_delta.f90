subroutine sample_delta(delta,H,share_h,y,sample_k,weights,joint_yh) 
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates_mixture,L_gender,L_educ,types),intent(inout)::delta
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
    real(DP),dimension(clusters,L_gender,L_educ,types)::fraction
    real(DP)::eps
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    call delta_2_fraction(delta,fraction)
    call compute_weights(fraction,H,share_h,weights,joint_yh) 
    
    !Compute likelihood of weights
    log_likeli=0.0d0
    do i_l=1,indv
        if (sample_k(i_l,first_age(i_l))>=1 .and. sample_k(i_l,first_age(i_l))<=2) then
            log_likeli(gender(i_l),educ(i_l))=log_likeli(gender(i_l),educ(i_l))+log(weights(first_age(i_l),sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l),y(i_l,1)))
        else
            log_likeli(gender(i_l),educ(i_l))=log_likeli(gender(i_l),educ(i_l))+log(weights(first_age(i_l),1,gender(i_l),educ(i_l),y(i_l,1)))
        end if
    end do
            
    !Generate new guess from normal proposal
    do e_l=1,L_educ;do ge_l=1,L_gender
        u=0.0d0
        do h_l=1,covariates_mixture;
            do y_l=1,types
                if (y_l<types) then
                    u(h_l,y_l)=c4_normal_01()/100.0d0        
                end if
                delta_g(h_l,ge_l,e_l,y_l)=delta(h_l,ge_l,e_l,y_l)+u(h_l,y_l)
            end do
        end do
        
        
    end do; end do
    
    call delta_2_fraction(delta_g,fraction)
    call compute_weights(fraction,H,share_h,weights_g,joint_yh_g) 
    
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
    
    
subroutine delta_2_fraction(delta,fraction)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates_mixture,L_gender,L_educ,types),intent(in)::delta
    real(DP),dimension(clusters,L_gender,L_educ,types),intent(out)::fraction
    integer::ge_l,e_l,t_l,h_l,y_l,y_l2
    real(DP),dimension(types)::y_star
    real(DP),dimension(covariates_mixture,1)::x
    integer,parameter::nodes=5
    real(DP),dimension(nodes):: xs=(/0.117581320211778,	1.0745620124369,	3.08593744371755,	6.41472973366203,	11.8071894899717/), &
                                weight=(/1.22172526747065,	0.480277222164629,	0.0677487889109621,	0.00268729149356246,	1.52808657104652E-05/),prod1,prod2
    
    do ge_l=1,L_gender; do e_l=1,L_educ; do h_l=1,clusters
        !By numerical integration
        x(:,1)=(/1.0d0, dble(h_l)-1.0d0/)
        do t_l=1,types
            y_star(t_l)=sum(x(:,1)*delta(:,ge_l,e_l,t_l))
        end do
        do y_l=1,types
            prod1=1
            prod2=1
            do y_l2=1,types
                if (y_l/=y_l2) then
                    prod1=prod1*0.5d0*(1.0d0+erf((-sqrt(2.0d0*xs)-(y_star(y_l2)-y_star(y_l)))/sqrt(2.0d0)))
                    prod2=prod2*0.5d0*(1.0d0+erf(( sqrt(2.0d0*xs)-(y_star(y_l2)-y_star(y_l)))/sqrt(2.0d0)))
                end if
            end do
            fraction(h_l,ge_l,e_l,y_l)=0.5d0/sqrt(pi)*sum(weight*(prod1+prod2))
        end do 
    end do;end do;end do
    
    
end subroutine