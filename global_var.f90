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
    
module global_var
    use nrtype
    implicit none
    !Define number of health clusters and health behavior types
    integer,parameter::clusters=2,types=3,cohorts=5 ! make cohort either 1 or 5
    integer,parameter::adls=12,habits=6,indv_HRS=27090,indv_PSID=7304,generations=37,initial_age=26,indv=indv_HRS+indv_PSID,L_gender=2,L_educ=3
    integer,parameter::covariates=2+(types-1)*2,covariates_habits=4,covariates_mixture=2+cohorts-1
    integer,parameter::g_max=10
    integer,dimension(indv,generations)::data_shlt
    integer,dimension(indv,habits,generations)::data_habits
    integer,dimension(indv)::first_age,last_age,gender,high_school,college,educ,race,birth_cohort
    real(DP),dimension(indv)::h_bar,a_bar
    integer,dimension(indv,1)::y_true
    real(DP),dimension(clusters,L_gender,L_educ)::share_h_true
    real(DP),dimension(clusters,L_gender,L_educ,types,cohorts)::delta_true
    character(LEN=42)::path="C:\Users\jbueren\Google Drive\endo_health\"  
    character(LEN=71)::path_s="C:\Users\jbueren\OneDrive - Istituto Universitario Europeo\endo_health\"
end module global_var
    
