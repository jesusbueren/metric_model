!subroutine E_part(p,pi_k,gamma_ik,E_L_new)
!    use global_var; use nrtype
!    implicit none
!    real(DP),dimension(adls,clusters),intent(in)::p
!    real(DP),dimension(clusters),intent(in)::pi_k
!    real(DP),dimension(indv,generations,clusters),intent(out)::gamma_ik
!    real(DP),intent(out)::E_L_new
!    real(DP),dimension(indv,clusters,generations)::L_i
!    real(DP),dimension(indv,generations,clusters)::E_L_i
!    integer::g_l,i_l,c_l
!    real(DP)::A
!    
!    call likelihood_all(p,L_i)
!    E_L_i=0.0d0
!    gamma_ik=0.0d0
!    !$OMP PARALLEL NUM_THREADS(4) DEFAULT(SHARED) PRIVATE(g_l,i_l,A,c_l)
!    !$OMP DO
!    do i_l=1,indv; do g_l=1,generations
!      if (data_adls(i_l,1,g_l) /= -1) then
!        A=sum(pi_k*L_i(i_l,:,g_l))
!        do c_l=1,clusters
!            gamma_ik(i_l,g_l,c_l)=(pi_k(c_l)*L_i(i_l,c_l,g_l))/A
!            E_L_i(i_l,g_l,c_l)=gamma_ik(i_l,g_l,c_l)*(log(pi_k(c_l))+log(L_i(i_l,c_l,g_l)))
!            if (isnan(E_L_i(i_l,g_l,c_l))) then
!                print*,'error'
!            end if
!        end do
!      end if
!    end do; end do
!    !$OMP END DO NOWAIT
!    !$OMP END PARALLEL
!
!    E_L_new=sum(E_L_i)
!    
!end subroutine
!    
!subroutine M_part(indv_all,gamma_ik,p,pi_k)
!    use global_var; use nrtype
!    implicit none
!    real(DP),dimension(adls,clusters),intent(out)::p
!    real(DP),dimension(clusters),intent(out)::pi_k
!    real(DP),dimension(indv,generations,clusters),intent(in)::gamma_ik
!    real(DP),intent(in)::indv_all
!    real(DP),dimension(clusters)::N_k
!    real(DP),dimension(indv,generations)::p_i
!    integer::c_l,v_l,g_l,i_l
!    
!    N_k=0
!    do c_l=1,clusters
!      N_k(c_l)=sum(gamma_ik(:,:,c_l))
!      pi_k(c_l)=N_k(c_l)/indv_all
!    end do
!
!    do c_l=1,clusters; do v_l=1,adls
!      p_i=0
!      do i_l=1,indv;do g_l=1,generations
!        if (data_adls(i_l,1,g_l)/=-1 .and. data_adls(i_l,v_l,g_l)/=-9) then
!          p_i(i_l,g_l)=gamma_ik(i_l,g_l,c_l)*data_adls(i_l,v_l,g_l)
!        elseif (data_adls(i_l,v_l,g_l)==-9) then
!          p_i(i_l,g_l)=gamma_ik(i_l,g_l,c_l)*1/dble(clusters)
!        end if
!      end do; end do
!      p(v_l,c_l)=1/N_k(c_l)*sum(p_i(:,:))
!    end do; end do
!    
!end subroutine
    