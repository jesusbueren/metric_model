subroutine sample_beta_d_MH(beta_d,beta_h,share_h,H,y,sample_k,weights,joint_yh,beta_d_mean,sigma_d,it,shrinkage_d)
    use nrtype; use global_var
    implicit none
    real(DP),dimension(covariates,clusters,L_gender,L_educ),intent(inout)::beta_d,beta_d_mean
    real(DP),dimension(covariates,clusters,L_gender,L_educ),intent(in)::beta_h
    integer,intent(in)::it
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(inout)::H
    integer,dimension(indv,1),intent(in)::y
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts),intent(inout)::weights,joint_yh
     real(DP),dimension(clusters,L_gender,L_educ),intent(in)::share_h
     real(DP),dimension(clusters,L_gender,L_educ),intent(inout)::shrinkage_d
    real(DP),dimension(clusters,L_gender,L_educ)::log_L,log_L_new
    integer::i_l,g_l,h_l,c_l,e_l,ge_l,t_l
    real(DP),dimension(covariates,covariates,clusters,L_gender,L_educ)::sigma
    real(DP),dimension(covariates,covariates,clusters,L_gender,L_educ),intent(inout)::sigma_d
    real(DP),dimension(covariates)::u,var_pro,mu_0,var_0
    real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_g
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H_g
    real(DP),dimension(types,L_gender,L_educ,clusters+1)::LE
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts)::joint_yh_g,weights_g
    real(DP)::u_MH
    real(DP),dimension(covariates,clusters,L_gender,L_educ)::old_mean
    
    ! build prior mean and variance vectors for all stacked types
    mu_0  = 0.0d0
    var_0 = 0.0d0
    do t_l = 1, types
        mu_0(2*(t_l-1)+1)  = -4.5d0
        mu_0(2*(t_l-1)+2)  = 0.05d0
        var_0(2*(t_l-1)+1) = 1.0d0**2
        var_0(2*(t_l-1)+2) = 0.02d0**2
    end do
    
    call compute_Likelihood_tr(H,y,sample_k,weights,log_L)
    
    if (it==1) then
        beta_d_mean=0.0d0
        sigma_d=0.0d0
    end if
    old_mean=beta_d_mean
    beta_d_mean=beta_d_mean+(beta_d-beta_d_mean)/dble(it)
    
    if (it>1) then
        do h_l=1,clusters; do ge_l=1,L_gender;do e_l=1,L_educ
            sigma_d(:,:,h_l,ge_l,e_l)=sigma_d(:,:,h_l,ge_l,e_l)+ &
                                      matmul(reshape(beta_d(:,h_l,ge_l,e_l)-old_mean(:,h_l,ge_l,e_l),[covariates,1]), reshape(beta_d(:,h_l,ge_l,e_l)-beta_d_mean(:,h_l,ge_l,e_l),[1,covariates])) 
        end do;end do;end do
    end if
    
    var_pro=reshape(spread((/1.0d-3,1.0d-7/), 2, types), (/ covariates /))

    sigma=0.0d0
    
    if (it>=5000) then
        if (it==5000) then
            shrinkage_d=1.0d0
        end if
        sigma=sigma_d/dble(it-1)*2.4**2.0d0/dble(covariates)
    else
        do c_l=1,covariates
            sigma(c_l,c_l,:,:,:)=var_pro(c_l)
        end do
    end if
    
    do h_l=1,clusters; do ge_l=1,L_gender;do e_l=1,L_educ
        do c_l=1,covariates
            call normal_01_sample ( u(c_l))
        end do
        sigma(:,:,h_l,ge_l,e_l)=sigma(:,:,h_l,ge_l,e_l)*shrinkage_d(h_l,ge_l,e_l)
        call choldc(sigma(:,:,h_l,ge_l,e_l),covariates)
        beta_g(:,h_l,ge_l,e_l)= matmul(sigma(:,:,h_l,ge_l,e_l),u) +beta_d(:,h_l,ge_l,e_l) 
    end do;end do;end do
    
    compute_LE=0
    call transitions(beta_h,beta_g,H_g,LE,joint_yh_g) 
    call compute_weights(weights(1,:,:,:,:,:),H_g,share_h,weights_g,joint_yh_g) 
    call compute_Likelihood_tr(H_g,y,sample_k,weights_g,log_L_new)
    
    do h_l=1,clusters; do ge_l=1,L_gender;do e_l=1,L_educ
        log_L(h_l,ge_l,e_l)=log_L(h_l,ge_l,e_l)-1.0d0/(2.0d0)*sum((beta_d(:,h_l,ge_l,e_l)-mu_0)**2.0d0/var_0)
        log_L_new(h_l,ge_l,e_l)=log_L_new(h_l,ge_l,e_l)-1.0d0/(2.0d0)*sum((beta_g(:,h_l,ge_l,e_l)-mu_0)**2.0d0/var_0)
        call RANDOM_NUMBER(u_mh)
        if (log(u_mh)<log_L_new(h_l,ge_l,e_l)-log_L(h_l,ge_l,e_l) .and. log_L(h_l,ge_l,e_l)/=0.0d0 ) then
            beta_d(:,h_l,ge_l,e_l)=beta_g(:,h_l,ge_l,e_l)
            acc_d(h_l,ge_l,e_l)=acc_d(h_l,ge_l,e_l)+1
        end if  
    end do;end do;end do
    
    compute_LE=0
    call transitions(beta_h,beta_d,H,LE,joint_yh_g)
    call compute_weights(weights(1,:,:,:,:,:),H,share_h,weights,joint_yh) 
  

    if (mod(it,100) == 0) then
        print*,'acc rate d %',acc_d(1,1,1)
        do h_l=1,clusters; do ge_l=1,L_gender;do e_l=1,L_educ
            if (acc_d(h_l,ge_l,e_l)<10) then
                shrinkage_d(h_l,ge_l,e_l)=shrinkage_d(h_l,ge_l,e_l)/1.5d0
            elseif (acc_d(h_l,ge_l,e_l)>60) then
                shrinkage_d(h_l,ge_l,e_l)=shrinkage_d(h_l,ge_l,e_l)*2.0
            end if
        end do;end do;end do
        acc_d=0
    end if
    
    end subroutine
    
    subroutine compute_Likelihood_tr(H,y,sample_k,weights,log_L)
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
                do g_l=first_age(i_l),last_age(i_l)-1
                    if (sample_k(i_l,g_l)>=1 .and. sample_k(i_l,g_l+1)>=1) then 
                        if (i_l<=indv_HRS) then
                            log_L(sample_k(i_l,g_l),gender(i_l),educ(i_l))=log_L(sample_k(i_l,g_l),gender(i_l),educ(i_l))+log(H(sample_k(i_l,g_l),sample_k(i_l,g_l+1),g_l,y(i_l,1),gender(i_l),educ(i_l)))
                        else
                            log_L(sample_k(i_l,g_l),gender(i_l),educ(i_l))=log_L(sample_k(i_l,g_l),gender(i_l),educ(i_l))+ &
                                                        log(H(sample_k(i_l,g_l),sample_k(i_l,g_l+1),g_l,y(i_l,1),gender(i_l),educ(i_l))/(1.0d0-H(sample_k(i_l,g_l),clusters+1,g_l,y(i_l,1),gender(i_l),educ(i_l))))
                        end if
                    end if
                end do
            end if
        end if
    end do
    
    
    
    end subroutine