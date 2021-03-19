subroutine compute_mean(mean,c,it,dim,new_mean)
    use nrtype
    implicit none
    integer,intent(in)::dim,it
    real(DP),dimension(dim,1),intent(in)::mean,c
    real(DP),dimension(dim,1),intent(out)::new_mean
    
    new_mean=dble(it)/dble(it+1)*mean+1.0d0/dble(it+1)*c 
    
end subroutine
    
subroutine compute_cov(C_A,it,s_d,mean,new_mean,c,dim)
use nrtype
    implicit none
    integer,intent(in)::dim,it
    real(DP),intent(in)::s_d
    real(DP),dimension(dim,dim),intent(inout)::C_A
    real(DP),dimension(dim,1),intent(in)::mean,new_mean,c
    real(DP)::it_d,eps=1.0D-8
    double precision,dimension(dim,dim)::I
    integer::c_l,c_l2
    
    I=0.0d0
    do c_l=1,dim; do c_l2=1,dim
        if (c_l==c_l2) then
            I(c_l,c_l2)=1.0d0
        end if
    end do; end do   
    
    it_d=dble(it)
    if (it==0) then
        C_A=0.0d0
    else
        C_A=(it_d-1)/it_d*C_A+s_d/it_d*(it_d*matmul(mean,transpose(mean))-(it_d+1)*matmul(new_mean,transpose(new_mean))+matmul(c,transpose(c))+eps*I)
    end if
    

end subroutine