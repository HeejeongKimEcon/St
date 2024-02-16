subroutine tauchen(meanz, stdinnov, rho, multiple, znum, zvec, piz)
use Lib_kindset
implicit none
  
! G. Tauchen (1986) 'Finite State Markov-Chain Approximations to 
! Univariate and Vector Autoregressions' Economics Leters 20: 177-181 
! 
! initially programmed by Aubhik Khan in MATLAB and translated to Fortran 90 by Carol Cui
! The program contains numerical integration of normal density to yield distribution.


!               meanz        mean of stochastic process z
!               stdinnov     standard deviation of innovations to z
!               rho          persistence of stochastic process z
!               multiple     numerical parameter that determines bounds as multiples of the standard deviation of z
!               znum         the number of grid points in the discretised space
!               zvec         the discretised space
!               piz          Markov transition matrix on zvec, piz(i,j) = Pr{z' = z(j) | z = z(i)}

! Bounds for tauchen grid are multiples of standard deviation of the stochastic process.

! parameters & variables
real(rk), intent(IN) :: meanz, stdinnov, rho, multiple
integer(ik), intent(IN) :: znum
real(rk), dimension(znum), intent(OUT) :: zvec
real(rk), dimension(znum,znum), intent(OUT) :: piz
real(rk) stdz,zlow,zhigh,meaninnov,lowerbound,w,z,F1,F0
integer(ik) i,j,k,gridsize

stdz = stdinnov**2_ik
stdz = stdz/(1.0_rk - rho**2_ik)
stdz = sqrt(stdz)
zlow = meanz - stdz*multiple
zhigh = meanz + stdz*multiple
meaninnov = meanz*(1.0_rk - rho)

! lowerbound and grid size for numerical integration
lowerbound = meaninnov - stdinnov*max(10.0_rk, 2.0_rk*multiple)
gridsize = 10000_ik

w = (zhigh - zlow)/(znum-1) ! increment of z
zvec(1) = zlow
zvec(znum) = zhigh
do i = 2,znum-1
zvec(i) = zvec(i-1)+w
end do

piz = 0.0

! This is equations (3a) and (3b) from Tauchen (1986)
do j = 1,znum
	
	z = zvec(1) - rho*zvec(j)
	F1 = normal(z + w/2.0_rk, meaninnov, stdinnov, gridsize, lowerbound)
	piz(j,1) = F1
	
	do k = 2,znum-1
		z = zvec(k) - rho*zvec(j)
		F1 = normal(z + w/2.0_rk, meaninnov, stdinnov, gridsize, lowerbound)
		F0 = normal(z - w/2.0_rk, meaninnov, stdinnov, gridsize, lowerbound)
		piz(j,k) = F1 - F0
    end do

	z = zvec(znum) - rho*zvec(j)
	F0 = normal(z - w/2.0_rk, meaninnov, stdinnov, gridsize, lowerbound)
	piz(j,znum) = 1.0_rk - F0

end do



contains
function normal(upperbound, mean, sd, gridsize, lowerbound)
! declaration
real(rk), intent(IN) :: upperbound, mean, sd, lowerbound
integer(ik), intent(IN) :: gridsize
real(rk), dimension(gridsize):: X0DL,X1DU,F0DL,F1DU,F
real(rk) normal,increment
integer(ik) i

! width of Darboux Rectangles
increment = (upperbound - lowerbound)/gridsize

! create left and right vertices for all Darboux rectangles
X0DL(1) = lowerbound
X0DL(gridsize) = upperbound-increment
do i = 2,gridsize-1
X0DL(i) = X0DL(i-1)+(X0DL(gridsize)-X0DL(1))/(gridsize-1_ik)
end do

X1DU(1) = lowerbound+increment
X1DU(gridsize) = upperbound
do i = 2,gridsize-1
X1DU(i) = X1DU(i-1)+(X1DU(gridsize)-X1DU(1))/(gridsize-1_ik)
end do

! obtain normal density at vertices
F0DL = normaldensity(X0DL, mean, sd)
F1DU = normaldensity(X1DU, mean, sd)

! Average the function values.
F = (F0DL + F1DU)/2.0_rk

! weight each one by the width of the rectangle of integration.
normal = sum(F*increment)

end function normal



function normaldensity(vectorofpoints, mean, sd)

! From the Matlab version of this program, the reference is apparently 
! page 244, chapter 4 of Bopwerman and O'Connell (1997) Applied Statistics

! declaration
real(rk), intent(IN) :: mean, sd
real(rk), dimension(:), intent(IN) :: vectorofpoints
real(rk), dimension(size(vectorofpoints)) :: transcend,normaldensity
real(rk) variance,pi,coefficient

variance = sd**2_ik
pi = 2.0_rk*asin(1.0_rk)
coefficient = variance*2.0_rk*pi
coefficient = 1.0_rk/sqrt(coefficient)

! note the vector-valued operations as a scalar, mean, is subtracted from each 
! element of a vector, then this resultant vector is divided by two scalars and
! finally used as a vector of exponents for the exponential function.

transcend = (vectorofpoints - mean)**2_ik;
transcend = transcend/variance;
transcend = transcend/2.0_rk;
transcend = -1.0_rk*transcend;
normaldensity = coefficient*exp(transcend);

end function normaldensity


end subroutine tauchen