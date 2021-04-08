subroutine load_high_density(c_mean,gamma_mean)
    use global_var;use nrtype
    implicit none
    double precision,dimension(covariates,clusters,clusters+1),intent(out)::c_mean
    real(DP),dimension(covariates_habits,habits,types),intent(out)::gamma_mean
    integer,parameter::burn=200
    integer,parameter::iterations=58640
    double precision,dimension(clusters**2*covariates,iterations)::c_vec
    double precision,dimension(covariates_habits*habits*types,iterations)::gamma_vec
        

    open(unit=9,file=path_s//'c_tr.txt')
        read(9,'(F20.5)') c_vec
    close(9)
    
    open(unit=9,file=path_s//'c_habits.txt')
        read(9,'(F20.5)') gamma_vec
    close(9)
        
    c_mean=0.0d0
    c_mean(:,:,1:clusters)=reshape(sum(c_vec(:,burn:iterations),2)/dble(iterations-burn+1),(/covariates,clusters,clusters/))
    
    gamma_mean=reshape(sum(gamma_vec(:,burn:iterations),2)/dble(iterations-burn+1),(/covariates_habits,habits,types/))

end subroutine