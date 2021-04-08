module global_var
    use nrtype
    implicit none
    !Define number of health clusters and health behavior types
    integer,parameter::clusters=2,types=2
    integer,parameter::adls=12,habits=5, indv=27090, generations=25, initial_age=50,L_gender=2,covariates_habits=3,L_educ=3
    integer,parameter::covariates=13
    integer,dimension(indv,generations)::data_shlt
    integer,dimension(indv,adls,generations)::data_adls
    integer,dimension(indv,habits,generations)::data_habits
    real(DP),dimension(indv*generations,types,habits,covariates_habits)::big_X
    real(DP),dimension(indv*generations,types,habits)::big_Y
    integer,dimension(types,habits)::counter_big_X
    
    real(DP),dimension(indv*generations,clusters,covariates)::big_X_h
    real(DP),dimension(indv*generations,clusters,clusters+1)::big_Y_h
    integer,dimension(clusters)::counter_big_X_h
    
    integer,dimension(indv)::first_age,last_age,gender,high_school,college,educ
    character(LEN=42)::path="C:\Users\jbueren\Google Drive\endo_health\"  
    character(LEN=71)::path_s="C:\Users\jbueren\OneDrive - Istituto Universitario Europeo\endo_health\"
end module
    
MODULE nrtype
	INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
	INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
	INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
	INTEGER, PARAMETER :: SP = KIND(1.0)
	INTEGER, PARAMETER :: DP = KIND(1.0D0)
	INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
	INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
	INTEGER, PARAMETER :: LGT = KIND(.true.)
	REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
	REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
	REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
	REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
	REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
	REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
	REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
	REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
	TYPE sprs2_sp
		INTEGER(I4B) :: n,len
		REAL(SP), DIMENSION(:), POINTER :: val
		INTEGER(I4B), DIMENSION(:), POINTER :: irow
		INTEGER(I4B), DIMENSION(:), POINTER :: jcol
	END TYPE sprs2_sp
	TYPE sprs2_dp
		INTEGER(I4B) :: n,len
		REAL(DP), DIMENSION(:), POINTER :: val
		INTEGER(I4B), DIMENSION(:), POINTER :: irow
		INTEGER(I4B), DIMENSION(:), POINTER :: jcol
	END TYPE sprs2_dp
END MODULE nrtype