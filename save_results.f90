subroutine save_results(beta,gamma,p,posterior_mu,LE,fraction_t,it)
    use global_var; use nrtype
    implicit none
    real(dp),dimension(covariates,clusters,clusters+1),intent(in)::beta
    real(dp),dimension(covariates_habits,habits,types),intent(in)::gamma
    real(DP),dimension(types,L_gender,L_educ,clusters+1),intent(in)::LE
    integer,intent(in)::it
    real(DP),dimension(types,L_gender,L_educ),intent(in)::fraction_t
    real(dp),dimension(clusters**2*covariates,1)::c_tr
    real(dp),dimension(covariates_habits*habits*types,1)::c_habits
    real(dp),dimension(adls,clusters)::p
    real(dp),intent(in)::posterior_mu
    real(DP),dimension(types*L_gender*L_educ*(clusters+1),1)::LE_v
    real(DP),dimension(types*L_gender*L_educ,1)::fraction_t_v
    
    c_tr=reshape(beta(:,:,1:clusters),(/covariates*clusters*clusters,1/))
    c_habits=reshape(gamma,(/covariates_habits*habits*types,1/))
    LE_v=reshape(LE,(/types*L_gender*L_educ*(clusters+1),1/))
    fraction_t_v=reshape(fraction_t,(/types*L_gender*L_educ,1/))
    
    if (it==1) then
        open(unit=9,file=path_s//'c_tr.txt')
            write(9,'(F20.8)') c_tr
        close(9)
        open(unit=9,file=path_s//'c_habits.txt')
            write(9,'(F20.8)') c_habits
        close(9)
        open(unit=9,file=path_s//'p.txt')
            write(9,'(F20.8)') p
        close(9)
        open(unit=9,file=path_s//'posterior.txt')
            write(9,'(F14.2,F14.2,F14.2)') posterior_mu
        close(9)
        open(unit=9,file=path_s//'LE.txt')
            write(9,'(F6.3)') LE_v
        close(9)
        open(unit=9,file=path_s//'fraction_t.txt')
            write(9,'(F7.4)') fraction_t_v
        close(9)
    else
        open(unit=9,file=path_s//'c_tr.txt',access='append')
            write(9,'(F20.8)') c_tr
        close(9)
        open(unit=9,file=path_s//'c_habits.txt',access='append')
            write(9,'(F20.8)') c_habits
        close(9)
        open(unit=9,file=path_s//'p.txt',access='append')
            write(9,'(F20.8)') p
        close(9)
        open(unit=9,file=path_s//'posterior.txt',access='append')
            write(9,'(F14.2,F14.2,F14.2)') posterior_mu
        close(9)
        open(unit=9,file=path_s//'LE.txt',access='append')
            write(9,'(F6.3)') LE_v
        close(9)
        open(unit=9,file=path_s//'fraction_t.txt',access='append')
            write(9,'(F7.4)') fraction_t_v
        close(9)
    end if
        
end subroutine