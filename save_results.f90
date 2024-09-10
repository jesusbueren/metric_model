subroutine save_results(beta_h,beta_d,gamma,delta,LE,fraction_t,fraction_h,H,it)
    use global_var; use nrtype
    implicit none
    real(dp),dimension(covariates,clusters,L_gender,L_educ),intent(in)::beta_h
    real(DP),dimension(covariates,clusters,L_gender,L_educ),intent(in)::beta_d
    real(dp),dimension(covariates_habits,habits,types),intent(in)::gamma
    real(DP),dimension(covariates_mixture,L_gender,L_educ,types),intent(in)::delta
    real(DP),dimension(types,L_gender,L_educ,clusters+1),intent(in)::LE
    integer,intent(in)::it
    real(DP),dimension(generations,L_gender,L_educ,types,cohorts),intent(in)::fraction_t
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts),intent(in)::fraction_h 
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
    

    if (it==1) then
        open(unit=9,file=path_s//'c_tr_'//types_s//'.txt')
            write(9,'(F20.8)') beta_h
        close(9)
        open(unit=10,file=path_s//'c_tr_d_'//types_s//'.txt')
            write(10,'(F20.8)') beta_d
        close(10)
        open(unit=11,file=path_s//'c_habits_'//types_s//'.txt')
            write(11,'(F20.8)') gamma
        close(11)
        open(unit=12,file=path_s//'LE_'//types_s//'.txt')
            write(12,'(F6.3)') LE
        close(12)
        open(unit=13,file=path_s//'fraction_t_'//types_s//'.txt')
            write(13,'(F7.4)') fraction_t
        close(13)
        open(unit=14,file=path_s//'H_'//types_s//'.txt')
            write(14,'(F7.4)') H
        close(14)
        open(unit=15,file=path_s//'fraction_h_'//types_s//'.txt')
            write(15,'(F7.4)') fraction_h(:,:,:,:,:,3)
        close(15)
        open(unit=16,file=path_s//'delta_'//types_s//'.txt')
            write(16,'(F20.8)') delta
        close(16)
    else
        open(unit=9,file=path_s//'c_tr_'//types_s//'.txt',access='append')
            write(9,'(F20.8)') beta_h
        close(9)
        open(unit=10,file=path_s//'c_tr_d_'//types_s//'.txt',access='append')
            write(10,'(F20.8)') beta_d
        close(10)
        open(unit=11,file=path_s//'c_habits_'//types_s//'.txt',access='append')
            write(11,'(F20.8)') gamma
        close(11)
        open(unit=12,file=path_s//'LE_'//types_s//'.txt',access='append')
            write(12,'(F6.3)') LE
        close(12)
        open(unit=13,file=path_s//'fraction_t_'//types_s//'.txt',access='append')
            write(13,'(F7.4)') fraction_t
        close(13)
        open(unit=14,file=path_s//'H_'//types_s//'.txt',access='append')
            write(14,'(F7.4)') H
        close(14)
        open(unit=15,file=path_s//'fraction_h_'//types_s//'.txt',access='append')
            write(15,'(F7.4)') fraction_h(:,:,:,:,:,3)
        close(15)
        open(unit=16,file=path_s//'delta_'//types_s//'.txt',access='append')
            write(16,'(F20.8)') delta
        close(16)
    end if
        
end subroutine