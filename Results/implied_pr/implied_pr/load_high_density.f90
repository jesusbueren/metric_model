subroutine load_high_density(c_tr,c_tr_d,gamma_mean)
    use global_var;use nrtype
    implicit none
    double precision,dimension(covariates,types,clusters,clusters),intent(out)::c_tr
    double precision,dimension(covariates,types,clusters),intent(out)::c_tr_d
    real(DP),dimension(covariates_habits,habits,types),intent(out)::gamma_mean
    integer,parameter::burn=200
    integer,parameter::iterations=3997
    double precision,dimension(clusters*(clusters-1)*types*covariates,iterations)::c_vec
    double precision,dimension(clusters*types*covariates,iterations)::c_vec_d
    double precision,dimension(covariates_habits*habits*types,iterations)::gamma_vec
        

    open(unit=9,file=path_s//'c_tr.txt')
        read(9,'(F20.5)') c_vec
    close(9)
    
    open(unit=9,file=path_s//'c_tr_d.txt')
        read(9,'(F20.5)') c_vec_d
    close(9)
    
    open(unit=9,file=path_s//'c_habits.txt')
        read(9,'(F20.5)') gamma_vec
    close(9)
        
    c_tr=0.0d0
    c_tr(:,:,:,1:clusters-1)=reshape(sum(c_vec(:,burn:iterations),2)/dble(iterations-burn+1),(/covariates,types,clusters,clusters-1/))
    
    c_tr_d=reshape(sum(c_vec_d(:,burn:iterations),2)/dble(iterations-burn+1),(/covariates,types,clusters/))
    
    
    gamma_mean=reshape(sum(gamma_vec(:,burn:iterations),2)/dble(iterations-burn+1),(/covariates_habits,habits,types/))

end subroutine