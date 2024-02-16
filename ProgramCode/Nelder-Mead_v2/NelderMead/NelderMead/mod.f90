

    SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter)
        use kindset
    implicit none
    !use amotry
    !use funk
    INTEGER(ik):: iter,mp,ndim,np,NMAX,ITMAX
    REAL(rk):: ftol,p(mp,np),y(mp),funk,TINY
    PARAMETER (NMAX=20,ITMAX=5000,TINY=1.e-10) !Maximum allowed dimensions and function evaluations, and a small number.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Multidimensional minimization of the function funk(x) where x(1:ndim) is a vector
    !in ndim dimensions, by the downhill simplex method of Nelder and Mead. The matrix
    !p(1:ndim+1,1:ndim) is input. Its ndim+1 rows are ndim-dimensional vectors which are
    !the vertices of the starting simplex. Also input is the vector y(1:ndim+1),
    !whose components must be pre-initialized to the values of funk evaluated at the ndim+1 vertices (rows)
    !of p; and ftol the fractional convergence tolerance to be achieved in the function value
    !(n.b.!). On output, p and y will have been reset to ndim+1 new points all within ftol of
    !a minimum function value, and iter gives the number of function evaluations taken.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    INTEGER(ik):: i,ihi,ilo,inhi,j,m,n
    REAL(rk):: rtol,sum,swap,ysave,ytry,psum(NMAX),amotry

    iter=0
1    do n=1,ndim    !Enter here when starting or have just overall contracted.
        sum=0.0_rk  !Recompute psum.
        do m=1,ndim+1
            sum=sum+p(m,n)
        end do
        psum(n)=sum
    end do

2   ilo=1_ik ! Enter here when have just changed a single point.
    if (y(1).gt.y(2)) then !Determine which point is the highest (worst), next-highest, and lowest (best).
        ihi=1
        inhi=2
    else
        ihi=2
        inhi=1
    end if

    do i=1,ndim+1 !by looping over the points in the simplex.
        if(y(i).le.y(ilo))  then
            ilo=i
        else
            if(y(i).gt.y(ihi)) then
                inhi=ihi
                ihi=i
            elseif(y(i).gt.y(inhi)) then
                if(i.ne.ihi) then
                    inhi=i
                end if
            end if
        end if
    end do
    rtol=2.0_rk*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
    !Compute the fractional range from highest to lowest and return if satisfactory.
    if (rtol.lt.ftol) then !If returning, put best point and value in slot 1.
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do n=1,ndim
            swap=p(1,n)
            p(1,n)=p(ilo,n)
            p(ilo,n)=swap
        end do
    return
    end if
    if (iter.ge.ITMAX) then 
        pause
        write(*,*) 'ITMAX exceeded in amoeba'
    else
        iter=iter+2 !Begin a new iteration. First extrapolate by a factor −1 through the face of the simplex across
                    !from the high point, i.e., reflect the simplex from the high point.
   ! ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0)
    
    if (ytry.le.y(ilo)) then   !Gives a result better than the best point, so try an additional extrapolation by a factor 2.
    !    ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0)
    elseif (ytry.ge.y(inhi)) then !The reflected point is worse than the second-highest, so look for an intermediate lower point,
                                   !i.e., do a one-dimensional contraction.
        ysave=y(ihi)
       ! ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5)
        if (ytry.ge.ysave) then !Can’t seem to get rid of that high point. Better contract around the lowest (best) point.
            do i=1,ndim+1     
                if(i.ne.ilo)then
                    do  j=1,ndim
                        psum(j)=0.5*(p(i,j)+p(ilo,j))
                        p(i,j)=psum(j)
                end do 
          !      y(i)=funk(psum)
        end if
    end do 
    iter=iter+ndim ! Keep track of function evaluations.
    goto 1 !Go back for the test of doneness and the next iteration.
    endif
    else
        iter=iter-1 !Correct the evaluation count.
    endif
    goto 2
    end if
    
    end subroutine amoeba
