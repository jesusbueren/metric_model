subroutine sample_beta_h_MH(beta_h,beta_d,share_h,H,y,sample_k,weights,joint_yh,beta_h_mean,sigma_h)
    use nrtype; use global_var
    implicit none
    real(DP),dimension(covariates,clusters,L_gender,L_educ),intent(inout)::beta_h,beta_h_mean
    real(DP),dimension(covariates,clusters,L_gender,L_educ),intent(in)::beta_d
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(inout)::H
    integer,dimension(indv,1),intent(in)::y
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts),intent(inout)::weights,joint_yh
     real(DP),dimension(clusters,L_gender,L_educ),intent(in)::share_h
    real(DP),dimension(clusters,L_gender,L_educ)::log_L,log_L_new
    integer::i_l,g_l,h_l,c_l,e_l,ge_l
    real(DP),dimension(covariates,covariates,clusters,L_gender,L_educ),intent(inout)::sigma_h
    real(DP),dimension(covariates,covariates,clusters,L_gender,L_educ)::sigma
    real(DP),dimension(types-1)::dummy_type,dummy_type_x_age
    real(DP),dimension(covariates)::u,var_proposal
    real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_g
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H_g
    real(DP),dimension(types,L_gender,L_educ,clusters+1)::LE
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts)::joint_yh_g,weights_g
    real(DP)::u_MH
    
    
    call compute_Likelihood_h(H,y,sample_k,weights,log_L)
    
    dummy_type=5.0d-6
    dummy_type_x_age=5.0d-6
    var_proposal=[(/1.0d-3,5.0d-6/),dummy_type,dummy_type_x_age]   
    

    beta_h_mean=dble(it_h2-1)/dble(it_h2)*beta_h_mean+dble(1)/dble(it_h2)*beta_h


    do h_l=1,clusters; do ge_l=1,L_gender;do e_l=1,L_educ
        sigma_h(:,:,h_l,ge_l,e_l)=dble(it_h2-1)/dble(it_h2)*sigma_h(:,:,h_l,ge_l,e_l)+dble(1)/dble(it_h2)*matmul(reshape(beta_h(:,clusters,L_gender,L_educ)-beta_h_mean(:,clusters,L_gender,L_educ),(/covariates,1/)),&
                                reshape(beta_h(:,clusters,L_gender,L_educ)-beta_h_mean(:,clusters,L_gender,L_educ),(/1,covariates/)))
    end do;end do;end do
    
    sigma=0.0d0
    sigma=sigma_h*1.0d-1
    do h_l=1,clusters; do ge_l=1,L_gender;do e_l=1,L_educ
        do c_l=1,covariates
            sigma(c_l,c_l,h_l,ge_l,e_l)=sigma(c_l,c_l,h_l,ge_l,e_l)+1.0d-6 !var_proposal(c_l)
            call normal_01_sample (  u(c_l))
        end do
        call choldc(sigma(:,:,h_l,ge_l,e_l),covariates)
        beta_g(:,h_l,ge_l,e_l)= matmul(sigma(:,:,h_l,ge_l,e_l),u) +beta_h(:,h_l,ge_l,e_l) 
    end do;end do;end do
    
    compute_LE=0
    call transitions(beta_g,beta_d,H_g,LE,joint_yh_g) 
    call compute_weights(weights(1,:,:,:,:,:),H_g,share_h,weights_g,joint_yh_g) 
    call compute_Likelihood_h(H_g,y,sample_k,weights_g,log_L_new)
    
    do h_l=1,clusters; do ge_l=1,L_gender;do e_l=1,L_educ
        call RANDOM_NUMBER(u_mh)
        if (log(u_mh)<log_L_new(h_l,ge_l,e_l)-log_L(h_l,ge_l,e_l)) then
            beta_h(:,h_l,ge_l,e_l)=beta_g(:,h_l,ge_l,e_l)
            if (h_l==1 .and. ge_l==1 .and. e_l==1) then
                acc_h=acc_h+1
                !print*,'accept draw'
            end if
        end if  
    end do;end do;end do
    
    compute_LE=0
    call transitions(beta_h,beta_d,H,LE,joint_yh_g)
    call compute_weights(weights(1,:,:,:,:,:),H,share_h,weights,joint_yh) 
    
    it_h=it_h+1
    it_h2=it_h2+1
    if (it_h==100) then
        print*,'acc rate h %',acc_h
        it_h=0
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
        if (sample_k(i_l,first_age(i_l))>=1) then
            log_L(sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l))=log_L(sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l))+log(weights(first_age(i_l),sample_k(i_l,first_age(i_l)),gender(i_l),educ(i_l),y(i_l,1),birth_cohort(i_l)))
        end if
        do g_l=first_age(i_l),last_age(i_l)-1
            if (sample_k(i_l,g_l)>=1 .and. sample_k(i_l,g_l+1)>=1 .and. sample_k(i_l,g_l+1)/=clusters+1) then 
                log_L(sample_k(i_l,g_l),gender(i_l),educ(i_l))=log_L(sample_k(i_l,g_l),gender(i_l),educ(i_l))+log(H(sample_k(i_l,g_l),sample_k(i_l,g_l+1),g_l,y(i_l,1),gender(i_l),educ(i_l))/(1.0d0-H(sample_k(i_l,g_l),clusters+1,g_l,y(i_l,1),gender(i_l),educ(i_l))))
            end if
        end do
    end do
    
    
    
    end subroutine