!subroutine logit(c_tr,posterior_gamma)
!    use global_ini;use global_var; use nrtype
!    implicit none
!    real(DP),dimension(covariates_habits*habits*types,1),intent(in)::c_tr
!    real(DP),intent(out)::posterior_gamma
!    integer::h_l,c_l,g_l,e_d,ge_l,age,ge_d,it,i_l,e_l,health_d
!    real(DP),dimension(habits,clusters,generations,types,L_gender)::alphas
!    real(DP),dimension(covariates_habits,habits,types)::gamma
!    real(DP),dimension(indv)::posterior_gamma_i
!    real(DP),dimension(covariates_habits,1)::x
!    
!    gamma=reshape(c_tr,(/covariates_habits,habits,types/)) 
!    do e_l=1,types
!        do h_l=1,habits;do c_l=1,clusters; do g_l=1,generations; do ge_l=1,2
!            age=initial_age+(g_l-1)*2
!            ge_d=ge_l-1
!            health_d=c_l-1
!            x(:,1)=(/1.0_dp,dble(age),dble(age**2.0_dp),dble(health_d),dble(ge_d),dble(age*ge_d),dble(age*health_d)/)
!            alphas(h_l,c_l,g_l,e_l,ge_l)=&
!        1.0_dp/(1.0_dp+exp(sum(x(:,1)*gamma(:,h_l,e_l))))
!        end do;end do; end do; end do
!    end do
!
!    do i_l=1,indv;
!        posterior_gamma_i(i_l)=1.0_dp
!        do g_l=first_age(i_l),last_age(i_l);do h_l=1,habits
!            if (data_habits(i_l,h_l,g_l)==1 ) then
!                posterior_gamma_i(i_l)=posterior_gamma_i(i_l)* &
!                alphas(h_l,sample_k_ini(i_l,g_l),g_l,y_ini(i_l,1)+1,gender(i_l))
!            elseif (data_habits(i_l,h_l,g_l)==0 ) then
!                posterior_gamma_i(i_l)=posterior_gamma_i(i_l)* &
!                (1.0_dp-alphas(h_l,sample_k_ini(i_l,g_l),g_l,y_ini(i_l,1)+1,gender(i_l)))
!            end if
!        end do; end do
!    end do
!    posterior_gamma=sum(log(posterior_gamma_i))
!end subroutine
!    
!subroutine prior_habits(c_tr,prior)
!    use global_var; use nrtype
!    real(DP),dimension(covariates_habits*habits*types,1),intent(in)::c_tr
!    real(DP),intent(out)::prior
!    integer::i
!    prior=0.0_dp
!    do i=1,covariates_habits*habits*types
!        prior=prior+log(1.0_dp/sqrt(2.0_dp*100.0_dp*4.0_dp*atan(1.0_dp)))-(c_tr(i,1)**2.0_dp/2.0_dp/100.0_dp)
!    end do
!end 
!    
!    