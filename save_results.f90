subroutine save_results(beta_h,beta_d,gamma,LE,fraction_t,H,it)
    use global_var; use nrtype
    implicit none
    real(dp),dimension(covariates,clusters,L_gender,L_educ),intent(in)::beta_h
    real(DP),dimension(covariates,clusters,L_gender,L_educ),intent(in)::beta_d
    real(dp),dimension(covariates_habits,habits,types),intent(in)::gamma
    real(DP),dimension(types,L_gender,L_educ,clusters+1),intent(in)::LE
    integer,intent(in)::it
    real(DP),dimension(generations,L_gender,L_educ,types,cohorts),intent(in)::fraction_t
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
    

    if (it==1) then
        open(unit=9,file=path_s//'c_tr.txt')
            write(9,'(F20.8)') beta_h
        close(9)
        open(unit=9,file=path_s//'c_tr_d.txt')
            write(9,'(F20.8)') beta_d
        close(9)
        open(unit=9,file=path_s//'c_habits.txt')
            write(9,'(F20.8)') gamma
        close(9)
        open(unit=9,file=path_s//'LE.txt')
            write(9,'(F6.3)') LE
        close(9)
        open(unit=9,file=path_s//'fraction_t.txt')
            write(9,'(F7.4)') fraction_t
        close(9)
        open(unit=9,file=path_s//'H.txt')
            write(9,'(F7.4)') H
        close(9)
    else
        open(unit=9,file=path_s//'c_tr.txt',access='append')
            write(9,'(F20.8)') beta_h
        close(9)
        open(unit=9,file=path_s//'c_tr_d.txt',access='append')
            write(9,'(F20.8)') beta_d
        close(9)
        open(unit=9,file=path_s//'c_habits.txt',access='append')
            write(9,'(F20.8)') gamma
        close(9)
        open(unit=9,file=path_s//'LE.txt',access='append')
            write(9,'(F6.3)') LE
        close(9)
        open(unit=9,file=path_s//'fraction_t.txt',access='append')
            write(9,'(F7.4)') fraction_t
        close(9)
        open(unit=9,file=path_s//'H.txt',access='append')
            write(9,'(F7.4)') H
        close(9)
    end if
        
end subroutine