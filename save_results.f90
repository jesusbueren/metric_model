subroutine save_results(beta_h,beta_d,gamma,gamma_med,delta,LE,fraction_t,fraction_h,H,it)
    use global_var; use nrtype
    implicit none
    real(dp),dimension(covariates,clusters,L_gender,L_educ),intent(in)::beta_h,beta_d
    real(dp),dimension(covariates_habits,habits_nomed,types),intent(in)::gamma
    real(dp),dimension(covariates_habits_med,habits_med,types),intent(in)::gamma_med
    real(DP),dimension(covariates_mixture*(types-1),L_gender,L_educ),intent(in)::delta
    real(DP),dimension(types,L_gender,L_educ,clusters+1),intent(in)::LE
    integer,intent(in)::it
    real(DP),dimension(generations,L_gender,L_educ,types,cohorts),intent(in)::fraction_t
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts),intent(in)::fraction_h 
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
    integer:: file_num
    
    file_num=9+(types-1)
    
    if (it==0) then
        open(unit=file_num,file=path_s//'c_tr_'//types_s//'.txt')
            write(file_num,'(F20.8)') beta_h
        close(file_num)
        open(unit=file_num,file=path_s//'c_tr_d_'//types_s//'.txt')
            write(file_num,'(F20.8)') beta_d
        close(file_num)
        open(unit=file_num,file=path_s//'c_habits_'//types_s//'.txt')
            write(file_num,'(F20.8)') gamma
        close(file_num)
        open(unit=file_num,file=path_s//'c_habits_med_'//types_s//'.txt')
            write(file_num,'(F20.8)') gamma_med
        close(file_num)
        open(unit=file_num,file=path_s//'LE_'//types_s//'.txt')
            write(file_num,'(F6.3)') LE
        close(12)
        open(unit=file_num,file=path_s//'fraction_t_'//types_s//'.txt')
            write(file_num,'(F7.4)') fraction_t !fraction_t(1,1,1,1,1)
        close(file_num)
        open(unit=file_num,file=path_s//'H_'//types_s//'.txt')
            write(file_num,'(F7.4)') H
        close(file_num)
        open(unit=file_num,file=path_s//'fraction_h_'//types_s//'.txt')
            write(file_num,'(F7.4)') fraction_h(:,:,:,:,:,3)
        close(file_num)
        open(unit=file_num,file=path_s//'delta_'//types_s//'.txt')
            write(file_num,'(F20.8)') delta
        close(file_num)
    else
        open(unit=file_num,file=path_s//'c_tr_'//types_s//'.txt',access='append')
            write(file_num,'(F20.8)') beta_h
        close(file_num)
        open(unit=file_num,file=path_s//'c_tr_d_'//types_s//'.txt',access='append')
            write(file_num,'(F20.8)') beta_d
        close(file_num)
        open(unit=file_num,file=path_s//'c_habits_'//types_s//'.txt',access='append')
            write(file_num,'(F20.8)') gamma
        close(file_num)
        open(unit=file_num,file=path_s//'c_habits_med_'//types_s//'.txt',access='append')
            write(file_num,'(F20.8)') gamma_med
        close(file_num)
        open(unit=file_num,file=path_s//'LE_'//types_s//'.txt',access='append')
            write(file_num,'(F6.3)') LE
        close(file_num)
        open(unit=file_num,file=path_s//'fraction_t_'//types_s//'.txt',access='append')
            write(file_num,'(F7.4)') fraction_t
        close(file_num)
        open(unit=file_num,file=path_s//'H_'//types_s//'.txt',access='append')
            write(file_num,'(F7.4)') H
        close(file_num)
        open(unit=file_num,file=path_s//'fraction_h_'//types_s//'.txt',access='append')
            write(file_num,'(F7.4)') fraction_h(:,:,:,:,:,3)
        close(file_num)
        open(unit=file_num,file=path_s//'delta_'//types_s//'.txt',access='append')
            write(file_num,'(F20.8)') delta
        close(file_num)
    end if
        
end subroutine