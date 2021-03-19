subroutine C_random_walk(C_rnd_W,dim,factor)
    use nrtype
    integer,intent(in)::dim
    real(dp),dimension(dim,dim),intent(out)::C_rnd_W
    real(dp),dimension(dim,1)::factor
    integer::c_l,c_l2
    
    C_rnd_W=0.0d0
    do c_l=1,dim
        do c_l2=1,dim
            if(c_l==c_l2) then
                C_rnd_W(c_l,c_l2)=factor(c_l,1)
            end if
        end do
    end do
    
end subroutine