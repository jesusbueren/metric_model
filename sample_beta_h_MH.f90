subroutine sample_beta_h_MH(beta_h,beta_d,share_h,H,y,sample_k,weights,joint_yh,beta_h_mean,sigma_h,it,shrinkage_h)
    use nrtype; use global_var
    implicit none
    real(DP),dimension(covariates,clusters,L_gender,L_educ),intent(inout)::beta_h,beta_h_mean
    real(DP),dimension(covariates,clusters,L_gender,L_educ),intent(in)::beta_d
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(inout)::H
    integer,dimension(indv,1),intent(in)::y
    integer,dimension(indv,generations),intent(in)::sample_k
    integer,intent(in)::it
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts),intent(inout)::weights,joint_yh
     real(DP),dimension(clusters,L_gender,L_educ),intent(in)::share_h
     real(DP),dimension(clusters,L_gender,L_educ),intent(inout)::shrinkage_h
    real(DP),dimension(clusters,L_gender,L_educ)::log_L,log_L_new
    integer::i_l,g_l,h_l,c_l,e_l,ge_l,t_l
    real(DP),dimension(covariates,covariates,clusters,L_gender,L_educ),intent(inout)::sigma_h
    real(DP),dimension(covariates,covariates,clusters,L_gender,L_educ)::sigma
    real(DP),dimension(covariates)::u,var_proposal
    real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_g
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H_g
    real(DP),dimension(types,L_gender,L_educ,clusters+1)::LE
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts)::joint_yh_g,weights_g
    real(DP)::u_MH
    real(DP),dimension(covariates,clusters,L_gender,L_educ)::old_mean
    
    
    call compute_Likelihood_h(H,y,sample_k,weights,log_L)
    
    if (it==51) then
        beta_h_mean=0.0d0
        sigma_h=0.0d0
    end if
    old_mean=beta_h_mean
    beta_h_mean=beta_h_mean+(beta_h-beta_h_mean)/dble(it-50)
    
    if (it>51) then
        do h_l=1,clusters; do ge_l=1,L_gender;do e_l=1,L_educ
            sigma_h(:,:,h_l,ge_l,e_l)=sigma_h(:,:,h_l,ge_l,e_l)+ &
                                      matmul(reshape(beta_h(:,h_l,ge_l,e_l)-old_mean(:,h_l,ge_l,e_l),[covariates,1]), reshape(beta_h(:,h_l,ge_l,e_l)-beta_h_mean(:,h_l,ge_l,e_l),[1,covariates])) 
        end do;end do;end do
    end if
    
    var_proposal=1.0d-8
    sigma=0.0d0
    
    if (it>=1000) then
        if (it==1000) then
            shrinkage_h=1.0d0
        end if
        sigma=sigma_h/dble(it-1-50)*2.4**2.0d0/dble(covariates)
    else
        do c_l=1,covariates
            sigma(c_l,c_l,:,:,:)=var_proposal(c_l)
        end do
    end if

    do h_l=1,clusters; do ge_l=1,L_gender;do e_l=1,L_educ
        do c_l=1,covariates
            call normal_01_sample (u(c_l))
        end do
        sigma(:,:,h_l,ge_l,e_l)=sigma(:,:,h_l,ge_l,e_l)*shrinkage_h(h_l,ge_l,e_l)
        call choldc(sigma(:,:,h_l,ge_l,e_l),covariates)
        beta_g(:,h_l,ge_l,e_l)= matmul(sigma(:,:,h_l,ge_l,e_l),u) +beta_h(:,h_l,ge_l,e_l) 
    end do;end do;end do
    
    compute_LE=0
    call transitions(beta_g,beta_d,H_g,LE,joint_yh_g) 
    call compute_weights(weights(1,:,:,:,:,:),H_g,share_h,weights_g,joint_yh_g) 
    call compute_Likelihood_h(H_g,y,sample_k,weights_g,log_L_new)
    
    do h_l=1,clusters; do ge_l=1,L_gender;do e_l=1,L_educ
        log_L(h_l,ge_l,e_l)=log_L(h_l,ge_l,e_l)-1.0d0/(2.0d0*10.0d0)*sum(beta_h(:,h_l,ge_l,e_l)**2.0d0)
        log_L_new(h_l,ge_l,e_l)=log_L_new(h_l,ge_l,e_l)-1.0d0/(2.0d0*10.0d0)*sum(beta_g(:,h_l,ge_l,e_l)**2.0d0)
        call RANDOM_NUMBER(u_mh)
        if (log(u_mh)<log_L_new(h_l,ge_l,e_l)-log_L(h_l,ge_l,e_l) .and. log_L(h_l,ge_l,e_l)/=0.0d0 ) then
            beta_h(:,h_l,ge_l,e_l)=beta_g(:,h_l,ge_l,e_l)
            acc_h(h_l,ge_l,e_l)=acc_h(h_l,ge_l,e_l)+1
        end if  
    end do;end do;end do
    
    compute_LE=0
    call transitions(beta_h,beta_d,H,LE,joint_yh_g)
    call compute_weights(weights(1,:,:,:,:,:),H,share_h,weights,joint_yh) 
    if (mod(it,100) == 0) then
        print*,'acc rate h %',acc_h(1,1,1)
        do h_l=1,clusters; do ge_l=1,L_gender;do e_l=1,L_educ
            if (acc_h(h_l,ge_l,e_l)<10) then
                shrinkage_h(h_l,ge_l,e_l)=shrinkage_h(h_l,ge_l,e_l)/1.5d0
            elseif (acc_h(h_l,ge_l,e_l)>60) then
                shrinkage_h(h_l,ge_l,e_l)=shrinkage_h(h_l,ge_l,e_l)*2.0
            end if
        end do;end do;end do
        acc_h=0
    end if
    
    end subroutine
    
    subroutine compute_Likelihood_h(H,y,sample_k,weights,log_L)
    use nrtype;use global_var
    implicit none
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
    integer,dimension(indv,1),intent(in)::y
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts),intent(in)::weights
    real(DP),dimension(clusters,L_gender,L_educ),intent(out)::log_L
    integer::i_l,g_l
    
    log_L=0.0d0
    do i_l=1,indv
        if (sample_selection(i_l)) then
            if (sample_k(i_l,first_age(i_l))>=1) then
                log_L(sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l))=log_L(sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l))+log(weights(first_age(i_l),sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l),y(i_l,1),birth_cohort(i_l)))
            end if
            do g_l=first_age(i_l),last_age(i_l)-1
                if (sample_k(i_l,g_l)>=1 .and. sample_k(i_l,g_l+1)>=1 .and. sample_k(i_l,g_l+1)/=clusters+1) then 
                    log_L(sample_k(i_l,g_l),gender(i_l),educ(i_l))=log_L(sample_k(i_l,g_l),gender(i_l),educ(i_l))+log(H(sample_k(i_l,g_l),sample_k(i_l,g_l+1),g_l,y(i_l,1),gender(i_l),educ(i_l))/(1.0d0-H(sample_k(i_l,g_l),clusters+1,g_l,y(i_l,1),gender(i_l),educ(i_l))))
                end if
            end do
        end if
    end do
    
    
    
    end subroutine