module get_params

! This module retrieves model parameters. 
!
use Lib_kindSet
implicit none

!precision for ergodic distribution
real(rk),    parameter  ::  precerg = 1e-10
real(rk),    parameter  ::  precision_gss = 1.0e-4;  !precision level for golden section search
real(rk),    parameter  ::  precision = 1.0e-7;      !precision level for labor bisection

!Change these to run different cases. Parameter values change depending on cases
integer(ik), parameter  :: iergovlini = 0_ik     ! if it is 0, re-draw labor productivity from low-skill distribution again.
integer(ik), parameter  :: iarray = 0_ik         ! h5 file. decision check 1- report 0- noreport
integer(ik), parameter  :: iarraysim = 0_ik      ! h5 file for simulation  1- report 0- noreport
integer(ik), parameter  :: itext = 0_ik          ! steady state result text file report
integer(ik), parameter  :: itran = 0_ik          ! steady state data for transition code
integer(ik), parameter  :: isim = 1_ik
integer(ik), parameter  :: itransim = 1_ik
integer(ik), parameter  :: inormal = 1_ik
integer(ik), parameter  :: irouwen = 0_ik
integer(ik), parameter  :: iwelfare = 0_ik

!parameters
integer(ik), parameter  ::  jnum= 67_ik         ! Final age (individual starts its age at 18 and retire at age 65 and live until the age 85)
integer(ik), parameter  ::  jrnum = 46_ik       ! Retirement age (individuals retire at jrnum+1)
integer(ik), parameter  ::  jc = 5_ik           ! Initial working age
integer(ik), parameter  ::  jd = 5_ik           ! Drop out age
integer(ik), parameter  ::  jpaynum = jrnum     ! jrnum + 1_ik    ! age at which individuals have no choice but to pay the student loan debt
integer(ik), parameter  ::  nT = 13_ik          ! final period of repayment

!parameters for grid
integer(ik), parameter  ::  epnum = 7_ik; 
integer(ik), parameter  ::  etnum = 3_ik; 
integer(ik), parameter  ::  e_num = epnum*etnum;
integer(ik), parameter  ::  ncost = 5_ik
integer(ik), parameter  ::  chibnum = 10_ik
real(rk),    parameter  ::  al = 0.0_rk;
real(rk),    parameter  ::  ah = 30.0_rk 
real(rk),    parameter  ::  blow = -6.0; 
real(rk),    parameter  ::  bhigh = 0.0_rk; 
real(rk),    parameter  ::  bbar = -0.24
real(rk),    parameter  ::  fam = 0.7_rk ! 1- percentage of family contribution in NLSY97                    
     
!parameter values
real(rk),    parameter  ::  beta = 0.95!5!0.95_rk                ! subjective discount factor
real(rk),    parameter  ::  sigma = 2.0                          ! CRRA parameter
real(rk),    parameter  ::  ssrr = 0.40_rk*0.3_rk                ! replacement rate of social security
real(rk),    parameter  ::  alpha = 0.36_rk                      ! capital share
real(rk),    parameter  ::  delta = 0.069_rk                     ! capital depreciation rate
real(rk),    parameter  ::  theta = 1.67_rk                      ! Elasticity of Subs. between skilled and unskilled labor (i.e. Krussell et al 2000)
real(rk),    parameter  ::  r  = 0.03_rk                         ! interest rate on savings
real(rk),    parameter  ::  tau = 0.0_rk!0.26_rk                 ! tax rate on labor income
real(rk),    parameter  ::  eta = 2.0!172_rk                     ! frische elasticity of labor supply 0.48=(1/sigma)*((1-n)/n)
real(rk),    parameter  ::  psi = 70.0_rk!1.85!1.22_rk           ! utility from leisure
real(rk),    parameter  ::  psi_c = 480.0_rk                     ! utility from leisure during college            
real(rk),    parameter  ::  mean = 0.0_rk
real(rk),    parameter  ::  multiple = 2.575_rk; 
real(rk),    parameter  ::  meana0 = 0.0_rk
real(rk),    parameter  ::  stda =1.0_rk 
real(rk),    parameter  ::  nworker = 1.0_rk
real(rk),    parameter  ::  ncollege = 1.0_rk!0.2_rk 

real(rk),    parameter  ::  nskilled = 0.3_rk!0.2750_rk
real(rk),    parameter  ::  nunskilled = 0.3_rk!0.3229_rk

real(rk),    parameter  ::  rho1 = 0.9792_rk!0.9784_rk; 
real(rk),    parameter  ::  rho2 = 0.9792_rk!0.9831_rk
real(rk),    parameter  ::  varp1 =  0.004765223 !0.01147! 0.026471175 !0.016552673!
real(rk),    parameter  ::  varp2 =  varp1!*5.1!don't go over 6.0 
real(rk),    parameter  ::  vart1 =  0.06852299  !0.0448!  0.078077301 !0.061126764! 
real(rk),    parameter  ::  vart2 =  vart1!
real(rk),    parameter  ::  wgar = 0.0_rk!0.025_rk !max is 15%. 

!data moments for normalization    
real(rk),    parameter ::   yaggdata = 7339.580987            !total GDP in billion from world bank data in 1984(2004 dollars)
real(rk),    parameter ::   gdp_data = 31122.99793            !GDP per capita in 1984 (2004 dollars)
real(rk),    parameter ::   totalpopdata = 235825000          !Total population in 1984from world bank

!parameters for transtion      
!index value to fix the sources
integer(ik), parameter  :: idefault = 0_ik
integer(ik), parameter  :: ipsycost = 1_ik
integer(ik), parameter  :: ituition = 1_ik
integer(ik), parameter  :: iwp = 1_ik
integer(ik), parameter  :: ishock = 1_ik
integer(ik), parameter  :: iability = 1_ik
integer(ik), parameter  :: ilaborfix = 0_ik
integer(ik), parameter  :: tnum = 5!100_ik
integer(ik), parameter  :: tmid = 5!2015-1984+1 ! = 32_ik
integer(ik), parameter  :: tsim = tmid!100_ik
end module get_params