subroutine load_high_density(p,c_mean,gamma_mean)
    use global_var;use nrtype
    implicit none
    double precision,dimension(adls,clusters),intent(out)::p
    double precision,dimension(covariates,clusters,clusters+1),intent(out)::c_mean
    real(DP),dimension(covariates_habits,habits,types),intent(out)::gamma_mean
    integer,parameter::burn=600000
    integer,parameter::iterations=696055
    double precision,dimension(adls,clusters,iterations)::p_all
    double precision,dimension(clusters**2*covariates,iterations)::c_vec
    double precision,dimension(covariates_habits*habits*types,iterations)::gamma_vec
    
    integer::c_l,c_l2,v_l,i,j
    
    open(unit=9,file=path_s//'p.txt')
        read(9,'(F20.5)') p_all
    close(9)
    open(unit=9,file=path_s//'c_tr.txt')
        read(9,'(F20.5)') c_vec
    close(9)
    
    open(unit=9,file=path_s//'c_habits.txt')
        read(9,'(F20.5)') gamma_vec
    close(9)
    
    !Compute high density point as mean of posterior density
    do v_l=1,adls; do c_l=1,clusters
        p(v_l,c_l)=sum(p_all(v_l,c_l,burn:iterations))/dble(iterations-burn+1)
    end do; end do
    
    c_mean=0.0d0
    c_mean(:,:,1:clusters)=reshape(sum(c_vec(:,burn:iterations),2)/dble(iterations-burn+1),(/covariates,clusters,clusters/))
    
    gamma_mean=reshape(sum(gamma_vec(:,burn:iterations),2)/dble(iterations-burn+1),(/covariates_habits,habits,types/))

end subroutine