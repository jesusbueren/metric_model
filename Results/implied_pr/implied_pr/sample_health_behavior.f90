subroutine sample_health_behavior(type_pr,y) 
    use global_var;use nrtype
    double precision,dimension(indv,types),intent(in)::type_pr
    integer,dimension(indv,1),intent(out)::y
    integer::i_l,ind
    real(DP)::u
    
    do i_l=1,indv;
        y(i_l,1)=-9
        call RANDOM_NUMBER(u)
        ind=1
        do while (y(i_l,1)==-9)
            if (u<sum(type_pr(i_l,1:ind),1) .or. ind==types) then
                y(i_l,1)=ind
            else
                ind=ind+1
            end if 
        end do
    end do
    

end subroutine