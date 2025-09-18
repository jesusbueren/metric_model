subroutine sample_delta(delta,mean_delta,cov_delta,H,share_h,y,sample_k,weights,joint_yh,it,acc_delta,shrinkage) 
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates_mixture,L_gender,L_educ,types),intent(inout)::delta
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
    real(DP),dimension(clusters,L_gender,L_educ),intent(in)::share_h
    integer,dimension(indv,1),intent(in)::y
    integer,dimension(indv,generations),intent(in)::sample_k
    integer,intent(in)::it
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts),intent(out)::weights,joint_yh
    real(DP),dimension(covariates_mixture,L_gender,L_educ,types),intent(inout)::mean_delta
    real(DP),dimension(covariates_mixture,covariates_mixture,L_gender,L_educ,types),intent(inout)::cov_delta
    integer,dimension(L_gender,L_educ,types),intent(inout)::acc_delta
    real(DP),dimension(L_gender,L_educ),intent(inout)::shrinkage
    real(DP),dimension(covariates_mixture,covariates_mixture)::sigma,eps_var
    real(DP),dimension(covariates_mixture,L_gender,L_educ,types)::delta_g,dif
    integer::i_l,e_l,ge_l,y_l,h_l,cov_l
    real(DP),dimension(L_gender,L_educ)::log_likeli,log_likeli_g
    real(DP),dimension(covariates_mixture)::u
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
        if (sample_k(i_l,first_age(i_l))>=1 .and. sample_k(i_l,first_age(i_l))<=2 .and. race(i_l)==1) then
            log_likeli(gender(i_l),educ(i_l))=log_likeli(gender(i_l),educ(i_l))+log(weights(first_age(i_l),sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l),y(i_l,1),birth_cohort(i_l)))
        elseif ( race(i_l)==1) then
            log_likeli(gender(i_l),educ(i_l))=log_likeli(gender(i_l),educ(i_l))+log(weights(first_age(i_l),1,gender(i_l),educ(i_l),y(i_l,1),birth_cohort(i_l)))
        end if
    end do

    log_likeli=log_likeli-sum(reshape(delta,[covariates_mixture*L_gender*L_educ*types,1])/(2.0d0*100.0d0))
    
            
    !Generate new guess from normal proposal
    delta_g=0.0d0
    do e_l=1,L_educ;do ge_l=1,L_gender;do y_l=1,types-1
            u=0.0d0            
            sigma=0.0d0
            eps_var=0.0d0
            do cov_l=1,covariates_mixture
                eps_var(cov_l,cov_l)=1.0d-8
                !if (it<200) then
                    sigma(cov_l,cov_l)=1d-4
                !end if
                u(cov_l)=c4_normal_01( )
            end do
            !if (it>=200) then
            !    sigma=cov_delta(:,:,ge_l,e_l,y_l)/dble(covariates_mixture)+eps_var
            !end if
            call choldc(sigma,covariates_mixture)
            delta_g(:,ge_l,e_l,y_l)=delta(:,ge_l,e_l,y_l)+matmul(sigma*shrinkage(ge_l,e_l),u)
        end do
    end do; end do
    
    call delta_2_fraction(delta_g,fraction)
    call compute_weights(fraction,H,share_h,weights_g,joint_yh_g) 
    
    !Compute likelihood of proposal
    log_likeli_g=0.0d0
    do i_l=1,indv
        if (sample_k(i_l,first_age(i_l))>=1 .and. sample_k(i_l,first_age(i_l))<=2 .and. race(i_l)==1) then
            log_likeli_g(gender(i_l),educ(i_l))=log_likeli_g(gender(i_l),educ(i_l))+log(weights_g(first_age(i_l),sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l),y(i_l,1),birth_cohort(i_l))) 
        elseif ( race(i_l)==1) then
            log_likeli_g(gender(i_l),educ(i_l))=log_likeli_g(gender(i_l),educ(i_l))+log(weights_g(first_age(i_l),1,gender(i_l),educ(i_l),y(i_l,1),birth_cohort(i_l)))
        end if
    end do
    log_likeli_g=log_likeli_g-sum(reshape(delta_g**2.0d0,[covariates_mixture*L_gender*L_educ*types,1])/(2.0d0*100.0d0))
    
    !Accept/reject proposal using Metropolis algorithm
    do e_l=1,L_educ;do ge_l=1,L_gender
        if(log_likeli_g(ge_l,e_l)>log_likeli(ge_l,e_l)) then
            delta(:,ge_l,e_l,:)=delta_g(:,ge_l,e_l,:)
            weights(:,:,ge_l,e_l,:,:)=weights_g(:,:,ge_l,e_l,:,:)
            joint_yh(:,:,ge_l,e_l,:,:)=joint_yh_g(:,:,ge_l,e_l,:,:)
            acc_delta(ge_l,e_l,1)=acc_delta(ge_l,e_l,1)+1
        else
           call RANDOM_NUMBER(eps)
           if (eps<exp(log_likeli_g(ge_l,e_l)-log_likeli(ge_l,e_l))) then
               delta(:,ge_l,e_l,:)=delta_g(:,ge_l,e_l,:)
                weights(:,:,ge_l,e_l,:,:)=weights_g(:,:,ge_l,e_l,:,:)
                joint_yh(:,:,ge_l,e_l,:,:)=joint_yh_g(:,:,ge_l,e_l,:,:)
                acc_delta(ge_l,e_l,1)=acc_delta(ge_l,e_l,1)+1
           end if
        end if            
    end do; end do   
    
    !if (it>100) then
    !    mean_delta=mean_delta+(delta-mean_delta)/dble(it-100)
    !    dif=delta-mean_delta
    !    if (it>101) then
    !        do e_l=1,L_educ;do ge_l=1,L_gender;do y_l=1,types-1
    !            cov_delta(:,:,ge_l,e_l,y_l)=dble(it-2-100)/dble(it-1-100)*cov_delta(:,:,ge_l,e_l,y_l)+ &
    !                                        dble(1)/dble(it-100)*matmul(reshape(dif(:,ge_l,e_l,y_l),[covariates_mixture,1]), reshape(dif(:,ge_l,e_l,y_l),[1,covariates_mixture]) )
    !        end do;end do;end do
    !    end if
    !end if
    
    if (mod(it,100) == 0) then
        print*,shrinkage(1,2)
        print*, "Iteration ", it, ": acceptance rate = ", dble(acc_delta(1,2,1))/100, ": cov = ", cov_delta(1,1,1,2,1),':shrinkage=',shrinkage(1,1),": dif = ",dif(1,1,2,1),":  L= ",log_likeli(1,1)
        do e_l=1,L_educ;do ge_l=1,L_gender
            if (dble(acc_delta(ge_l,e_l,1))/100.0d0>0.4d0 .and. it<5000) then
                shrinkage(ge_l,e_l)=shrinkage(ge_l,e_l)*1.3d0
            elseif (dble(acc_delta(ge_l,e_l,1))/100.0d0<0.2d0 .and. it<5000) then
                shrinkage(ge_l,e_l)=shrinkage(ge_l,e_l)*0.8d0
            end if 
        end do; end do
        acc_delta=0
    end if
    
    !print*,'delta',delta(1,1,1,1)
    !print*,'cov',cov_delta(1,1,1,1,1)
    
end subroutine

subroutine delta_2_fraction(delta,fraction)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates_mixture,L_gender,L_educ,types),intent(in)::delta
    real(DP),dimension(clusters,L_gender,L_educ,types,cohorts),intent(out)::fraction
    integer::ge_l,e_l,t_l,h_l,y_l,y_l2,co_l
    real(DP),dimension(types)::y_star
    real(DP),dimension(cohorts)::cohort_d
    real(DP),dimension(covariates_mixture,1)::x
    integer,parameter::nodes=5
    real(DP),dimension(nodes):: xs=(/0.117581320211778,	1.0745620124369,	3.08593744371755,	6.41472973366203,	11.8071894899717/), &
                                weight=(/1.22172526747065,	0.480277222164629,	0.0677487889109621,	0.00268729149356246,	1.52808657104652E-05/),prod1,prod2
    
    do ge_l=1,L_gender; do e_l=1,L_educ; do h_l=1,clusters; do co_l=1,cohorts
        !By numerical integration
        cohort_d=0.0d0
        cohort_d(co_l)=1.0d0
        x(:,1)=[dble(h_l)-1.0d0,cohort_d]
        do t_l=1,types
            y_star(t_l)=sum(x(:,1)*delta(:,ge_l,e_l,t_l))
        end do
        do y_l=1,types
            prod1=1.0d0
            prod2=1.0d0
            do y_l2=1,types
                if (y_l/=y_l2) then
                    prod1=prod1*0.5d0*(1.0d0+erf((-sqrt(2.0d0*xs)-(y_star(y_l2)-y_star(y_l)))/sqrt(2.0d0)))
                    prod2=prod2*0.5d0*(1.0d0+erf(( sqrt(2.0d0*xs)-(y_star(y_l2)-y_star(y_l)))/sqrt(2.0d0)))
                end if
            end do
            fraction(h_l,ge_l,e_l,y_l,co_l)=max(0.5d0/sqrt(pi)*sum(weight*(prod1+prod2)),1.0d-8)
        end do 
        fraction(h_l,ge_l,e_l,:,co_l)=fraction(h_l,ge_l,e_l,:,co_l)/sum(fraction(h_l,ge_l,e_l,:,co_l))
    end do;end do;end do;end do
    
    
end subroutine
    
    
double precision function gamma_rng(alpha)
    use nrtype
    implicit none
    real(DP), intent(in) :: alpha
    real(DP) ::  u
    interface
        double precision function gamma_rng_core(alpha) 
        use nrtype
        implicit none
        real(DP), intent(in) :: alpha
        end function
    end interface

    if (alpha <= 0.d0) then
        print*, "Error: alpha must be > 0"
        stop
    end if

    if (alpha < 1.d0) then
        ! For alpha < 1: boost to alpha+1, then adjust
        call random_number(u)
        gamma_rng = gamma_rng_core(alpha + 1.d0) * u**(1.d0/alpha)
    else
        gamma_rng = gamma_rng_core(alpha)
    end if
end function gamma_rng

    
double precision function  gamma_rng_core(alpha) 
    use nrtype
    implicit none
    real(DP), intent(in) :: alpha
    real(DP) :: g, d, c, x, v, u
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface

    d = alpha - 1.d0/3.d0
    c = 1.d0 / sqrt(9.d0*d)

    do
        x = c4_normal_01( )
        v = (1.d0 + c*x)**3
        if (v > 0.d0) then
            call random_number(u)
            if (u < 1.d0 - 0.0331d0*(x**4)) then
                gamma_rng_core = d*v
                return
            end if
            if (log(u) < 0.5d0*x**2 + d*(1.d0 - v + log(v))) then
                gamma_rng_core = d*v
                return
            end if
        end if
    end do
end function gamma_rng_core
    

