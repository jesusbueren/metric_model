subroutine sample_delta(delta,mean_delta,cov_delta,H,share_h,y,sample_k,weights,joint_yh,it,acc_delta,shrinkage) 
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates_mixture*(types-1),L_gender,L_educ),intent(inout)::delta
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
    real(DP),dimension(clusters,L_gender,L_educ),intent(in)::share_h
    integer,dimension(indv,1),intent(in)::y
    integer,dimension(indv,generations),intent(in)::sample_k
    integer,intent(in)::it
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts),intent(out)::weights,joint_yh
    real(DP),dimension(covariates_mixture*(types-1),L_gender,L_educ),intent(inout)::mean_delta
    real(DP),dimension(covariates_mixture*(types-1),covariates_mixture*(types-1),L_gender,L_educ),intent(inout)::cov_delta
    integer,dimension(L_gender,L_educ),intent(inout)::acc_delta
    real(DP),dimension(L_gender,L_educ),intent(inout)::shrinkage
    real(DP),dimension(covariates_mixture*(types-1),L_gender,L_educ)::old_mean,dif
    real(DP),dimension(covariates_mixture*(types-1),covariates_mixture*(types-1))::sigma,eps_var
    real(DP),dimension(covariates_mixture*(types-1),L_gender,L_educ)::delta_g
    integer::i_l,e_l,ge_l,y_l,h_l,cov_l
    real(DP),dimension(L_gender,L_educ)::log_likeli,log_likeli_g
    real(DP),dimension(covariates_mixture*(types-1))::u
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts)::weights_g,joint_yh_g
    real(DP),dimension(clusters,L_gender,L_educ,types,cohorts)::fraction
    real(DP)::eps
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
        double precision function gamma_rng(alpha)
            use nrtype
            implicit none
            real(DP), intent(in) :: alpha
        end function
    end interface
    
    call delta_2_fraction(delta,fraction)
    call compute_weights(fraction,H,share_h,weights,joint_yh) 
    
    !Compute likelihood of weights
    log_likeli=0.0d0
    do i_l=1,indv
        if (sample_selection(i_l)) then
            if (sample_k(i_l,first_age(i_l))>=1 .and. sample_k(i_l,first_age(i_l))<=2 ) then
                log_likeli(gender(i_l),educ(i_l))=log_likeli(gender(i_l),educ(i_l))+log(weights(first_age(i_l),sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l),y(i_l,1),birth_cohort(i_l))) !weights_g(first_age(97),sample_k(97,first_age(97)),gender(97),educ(97),y(97,1),birth_cohort(97))
            end if
        end if
    end do
            
    !Generate new guess from normal proposal
    delta_g=0.0d0
    do e_l=1,L_educ;do ge_l=1,L_gender
            u=0.0d0            
            sigma=0.0d0
            eps_var=0.0d0
            do cov_l=1,covariates_mixture*(types-1)
                eps_var(cov_l,cov_l)=1.0d-8
                u(cov_l)=c4_normal_01( )
                sigma(cov_l,cov_l)=1.0d-3
            end do
            
            if (it>=500) then
                if (it==500) then
                    shrinkage(ge_l,e_l)=1.0d0
                end if
                sigma=cov_delta(:,:,ge_l,e_l)/dble(it-1)*2.4**2.0d0/dble(covariates_mixture*(types-1))
            end if
            sigma=sigma+eps_var
            call choldc(sigma,covariates_mixture*(types-1))
            delta_g(:,ge_l,e_l)=delta(:,ge_l,e_l)+matmul(sigma*shrinkage(ge_l,e_l),u)
        end do
    end do
    
    call delta_2_fraction(delta_g,fraction)
    call compute_weights(fraction,H,share_h,weights_g,joint_yh_g) 
    
    !Compute likelihood of proposal
    log_likeli_g=0.0d0
    do i_l=1,indv
        if (sample_selection(i_l)) then
            if (sample_k(i_l,first_age(i_l))>=1 .and. sample_k(i_l,first_age(i_l))<=2) then
                log_likeli_g(gender(i_l),educ(i_l))=log_likeli_g(gender(i_l),educ(i_l))+log(weights_g(first_age(i_l),sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l),y(i_l,1),birth_cohort(i_l))) 

            end if
        end if
    end do


    !Accept/reject proposal using Metropolis algorithm
    do e_l=1,L_educ;do ge_l=1,L_gender
        log_likeli(ge_l,e_l)=log_likeli(ge_l,e_l)-1.0d0/(2.0d0*10.0d0)*sum(delta(:,ge_l,e_l)**2.0d0)
        log_likeli_g(ge_l,e_l)=log_likeli_g(ge_l,e_l)-1.0d0/(2.0d0*10.0d0)*sum(delta_g(:,ge_l,e_l)**2.0d0)
        call RANDOM_NUMBER(eps)
        if (log(eps)<log_likeli_g(ge_l,e_l)-log_likeli(ge_l,e_l) .and. log_likeli(ge_l,e_l)/=0.0d0) then
            delta(:,ge_l,e_l)=delta_g(:,ge_l,e_l)
            weights(:,:,ge_l,e_l,:,:)=weights_g(:,:,ge_l,e_l,:,:)
            joint_yh(:,:,ge_l,e_l,:,:)=joint_yh_g(:,:,ge_l,e_l,:,:)
            acc_delta(ge_l,e_l)=acc_delta(ge_l,e_l)+1
        end if          
    end do; end do   
    

    old_mean=mean_delta
    mean_delta=mean_delta+(delta-mean_delta)/dble(it)
    dif=delta-mean_delta
    if (it>1) then
        do e_l=1,L_educ;do ge_l=1,L_gender
            cov_delta(:,:,ge_l,e_l)=cov_delta(:,:,ge_l,e_l)+ &
                                        matmul(reshape(delta(:,ge_l,e_l)-old_mean(:,ge_l,e_l),[covariates_mixture*(types-1),1]), reshape(delta(:,ge_l,e_l)-mean_delta(:,ge_l,e_l),[1,covariates_mixture*(types-1)]))  
        end do;end do
    end if

    
    if (mod(it,100) == 0) then
        print*, "Iteration ", it
        print*, ": acceptance rate = ", dble(acc_delta(1,1))/100
        print*, ": cov = ", cov_delta(1,1,1,1)/dble(it), cov_delta(2,2,1,1)/dble(it), cov_delta(1,2,1,1)/dble(it)
        print*, ": shrinkage = ", shrinkage(1,1)
        print*,": delta = ",delta(1,1,1),delta_g(1,1,1)
        print*,":  L= ",log_likeli(1,1),log_likeli_g(1,1)
        do e_l=1,L_educ;do ge_l=1,L_gender
            if (dble(acc_delta(ge_l,e_l))/100.0d0>0.6d0) then
                shrinkage(ge_l,e_l)=shrinkage(ge_l,e_l)*2.0d0
            elseif (dble(acc_delta(ge_l,e_l))/100.0d0<0.15d0 ) then
                shrinkage(ge_l,e_l)=shrinkage(ge_l,e_l)/1.5d0
            end if 
        end do; end do
        acc_delta=0
    end if
    
    !print*,'delta',delta(1,1,1,1)
    !print*,'cov',cov_delta(1,1,1,1,1)
    
    end subroutine

subroutine delta_2_fraction(delta_in,fraction)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates_mixture*(types-1),L_gender,L_educ),intent(in)::delta_in
    real(DP),dimension(clusters,L_gender,L_educ,types,cohorts),intent(out)::fraction
    real(DP),dimension(covariates_mixture,types)::delta
    integer::ge_l,e_l,t_l,h_l,y_l,y_l2,co_l
    real(DP),dimension(types)::y_star
    real(DP),dimension(cohorts)::cohort_d
    real(DP),dimension(covariates_mixture,1)::x
    integer,parameter::nodes=5
    real(DP),dimension(nodes):: xs=(/0.117581320211778,	1.0745620124369,	3.08593744371755,	6.41472973366203,	11.8071894899717/), &
                                weight=(/1.22172526747065,	0.480277222164629,	0.0677487889109621,	0.00268729149356246,	1.52808657104652E-05/),prod1,prod2
    
    do ge_l=1,L_gender; do e_l=1,L_educ; do h_l=1,clusters; do co_l=1,cohorts
        !reshape
        delta(:,1:types-1)=reshape(delta_in(:,ge_l,e_l),(/covariates_mixture,types-1/))
        delta(:,types)=0.0d0
        !By numerical integration
        cohort_d=0.0d0
        cohort_d(co_l)=1.0d0
        x(:,1)=[dble(h_l)-1.0d0,cohort_d]
        do t_l=1,types
            y_star(t_l)=sum(x(:,1)*delta(:,t_l))
        end do
        do y_l=1,types
            !prod1=1.0d0
            !prod2=1.0d0
            !do y_l2=1,types
            !    if (y_l/=y_l2) then
            !        prod1=prod1*0.5d0*(1.0d0+erf((-sqrt(2.0d0*xs)-(y_star(y_l2)-y_star(y_l)))/sqrt(2.0d0)))
            !        prod2=prod2*0.5d0*(1.0d0+erf(( sqrt(2.0d0*xs)-(y_star(y_l2)-y_star(y_l)))/sqrt(2.0d0)))
            !    end if
            !end do
            fraction(h_l,ge_l,e_l,y_l,co_l)=1.0d0/sum(exp(y_star(:)-y_star(y_l))) !max(0.5d0/sqrt(pi)*sum(weight*(prod1+prod2)),1.0d-8)
        end do 
        fraction(h_l,ge_l,e_l,:,co_l)=fraction(h_l,ge_l,e_l,:,co_l)/sum(fraction(h_l,ge_l,e_l,:,co_l))
        if (isnan(sum(fraction(h_l,ge_l,e_l,:,co_l)))) then
            print*,'problem delta2 fraction'
        end if
    end do;end do;end do;end do
    
    
end subroutine

    
    

    

