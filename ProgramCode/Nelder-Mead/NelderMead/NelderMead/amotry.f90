    FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
    use kindset
    INTEGER(ik):: ihi,mp,ndim,np,NMAX
    REAL(rk):: amotry,fac,p(mp,np),psum(np),y(mp),funk
    PARAMETER (NMAX=20)
    EXTERNAL funk
    
    !Extrapolates by a factor fac through the face of the simplex across from the high point,
    !tries it, and replaces the high point if the new point is better.
    
    INTEGER(ik):: j
    REAL(rk):: fac1,fac2,ytry,ptry(NMAX)
    fac1=(1.0_rk-fac)/ndim
    fac2=fac1-fac
    do  j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
    end do 
    ytry=funk(ptry) !Evaluate the function at the trial point.
    if (ytry.lt.y(ihi)) then !If it’s better than the highest, then replace the highest.
        y(ihi)=ytry
        do  j=1,ndim
            psum(j)=psum(j)-p(ihi,j)+ptry(j)
            p(ihi,j)=ptry(j)
    end do 
    endif
    amotry=ytry
    return
    END