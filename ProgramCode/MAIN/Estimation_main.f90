    program estimation_main

    use estimation_mod
    use Lib_kindset
    use Lib_conshock
    use Lib_inequality
    use Lib_sort
    use Lib_grid
    use omp_lib
    use get_params

    !To optimize programs to run faster, set the following
    !Configuation Properties\Fortran\Run-time\set Check Stack Frame to no
    !Fortran\General\Set Debug information format to none
    !Fortran\ Optimization to Maximize Speed
    !Use Intel Math Kernel Library in Properties\Fortran\Libraries
    !In properties\linker\manifest\generate manifest file 'NO'

    implicit none

    integer(ik) :: nabil, anum, bnum, adnum, bdnum, olgo, nparam, nsimp, ntarg, maxiters,shrink, thisreport, thisstep
    integer(ik) :: j, i, t, ib, ief,  nsim, chinum, ia, ie, iteration,  ic,   iread, ichib
    integer(ik), allocatable:: indices(:)


    real(rk) ::  r_b,  wp, wl, wh,   min, max, lb, ub,  chid, tbar,  pi, chir, chihigh, ebar, &
        eval, psi1, psi2, psi3, rhonm, chinm, gammanm, sigmanm, fval, distance, &
        fr, fe, fc, fsol, temp, temp2, temp3, phic0, phic1, phic2, phic3, phic4, phi0, phi1, phi2, scale10, scale100, scale1000, &
        chid1, chid2, chid3, chid4, chid5, chib

    real(rk), allocatable ::  piea0(:,:),  pop(:), lambda(:), agrid(:), adgrid(:), bgrid(:), bdgrid(:),   &
        egridc(:),  egridl(:), piec(:,:), tempa(:), piel(:,:),  a0grid(:,:), piel0(:), piec0(:),  pia0(:),  &
        expel(:), expeh(:),  netphi(:,:),  psycost(:,:,:), chi_random(:), pia(:,:), x0(:), x(:,:), targ(:), weights(:),  &
        xold(:,:), fold(:), xr(:), xe(:), xc(:), v(:,:), f(:), fv(:), xbar(:), xsol(:), probcost(:,:,:), ctemp(:), &
        scalev(:), chidage(:), eagridj(:,:,:,:), abilgrid(:), pabil(:), chibprob(:), chibshock(:)

    character(30):: resultfile, datafile, datestring, timestring

    call date_and_time(datestring, timestring)

    resultfile = 'stddebt_est.txt'
    datafile = 'KK_stddebt_'
    olgo = len_trim(datafile)
    datafile(olgo + 1:olgo+4) = timestring(1:4)



    !!!!!Read parameters from steady state text file!!!!!
    OPEN(33, FILE='ss_par.txt', STATUS='OLD', ACTION='READ')
    read(33,*) r_b
    read(33,*) wp
    read(33,*) wl
    read(33,*) wh
    read(33,*) tbar
    read(33,*) pi
    read(33,*) chir
    read(33,*) psi1
    read(33,*) psi2
    read(33,*) psi3
    read(33,*) phi0
    read(33,*) phi1
    read(33,*) phi2
    read(33,*) phic0
    read(33,*) phic1
    read(33,*) phic2
    read(33,*) phic3
    read(33,*) phic4
    read(33,*) chid1
    read(33,*) chid2
    read(33,*) chid3
    read(33,*) chid4
    read(33,*) chid5
    read(33,*) anum
    read(33,*) bnum
    read(33,*) adnum
    read(33,*) bdnum
    read(33,*) ebar
    read(33,*) chib
    read(33,*) nabil

    allocate(pop(jnum), lambda(nT), expel(jrnum), expeh(jrnum), chidage(jnum))

    do j = 1,jnum
        read(33,*) chidage(j)
    end do
    do j = 1,jnum
        read(33,*) pop(j)
    end do
    do i = 1,nT
        read(33,*) lambda(i)
    end do
    do j = 1,jrnum
        read(33,*) expel(j)
    end do
    do j = 1,jrnum
        read(33,*) expeh(j)
    end do
    allocate(agrid(anum), bgrid(bnum), adgrid(adnum), bdgrid(bdnum))
    do ia = 1,anum
        read(33,*) agrid(ia)
    end do
    do ib = 1,bnum
        read(33,*) bgrid(ib)
    end do
    do ia = 1,adnum
        read(33,*) adgrid(ia)
    end do
    do ib = 1,bdnum
        read(33,*) bdgrid(ib)
    end do
    allocate(egridc(e_num), egridl(e_num), piec(e_num,e_num), piel(e_num,e_num), piel0(e_num), piec0(e_num), &
        tempa(anum), pia0(anum),  eagridj(jnum,nabil,e_num,2), abilgrid(nabil), pabil(nabil))
    do ia= 1, anum
        read(33,*) tempa(ia)
    end do
    do ia= 1,anum
        read(33,*) pia0(ia)
    end do

    do ie = 1, e_num
        read(33,*) egridc(ie)
    end do
    do ie = 1, e_num
        read(33,*) egridl(ie)
    end do
    do ie = 1, e_num
        do ief = 1, e_num
            read(33,*) piec(ie, ief)
        end do
    end do
    do ie = 1, e_num
        do ief = 1, e_num
            read(33,*) piel(ie, ief)
        end do
    end do
    do ie = 1, e_num
        read(33,*) piel0(ie)
    end do
    do ie = 1, e_num
        read(33,*) piec0(ie)
    end do
    do j = 1, jnum
        do ia = 1, nabil
            do ie = 1,e_num
                do ief = 1, 2
                    read(33,*) eagridj(j,ia,ie,ief)
                end do
            end do
        end do
    end do
    do ie = 1, nabil
        read(33,*) abilgrid(ie)
    end do
    do ie = 1, nabil
        read(33,*) pabil(ie)
    end do

    allocate(chibprob(chibnum), chibshock(chibnum))
    do ichib = 1, chibnum
        read(33,*) chibshock(ichib)
    end do

    do ichib = 1, chibnum
        read(33,*) chibprob(ichib)
    end do
    close(33)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !!     start estimation here
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nsim = 10000   ! number of simulation
    t =0_ik
    iread = 0_ik    ! read from the previous estimation

    !Nelder-Mead parameter
    rhonm = 1.0_rk      !parameter for relfection
    chinm = 2.0_rk      !parameter for expansion
    gammanm = 0.5_rk    !parameter for contraction
    sigmanm = 0.5_rk    !parameter for shrink

    nparam = 12_ik       !number of estimating parameters
    nsimp = nparam+1_ik !dimension for simplex
    ntarg = 13_ik       !number of moments to target

    !create simplex
    allocate(x0(nparam), x(nsimp,nparam), targ(ntarg), weights(ntarg), f(nsimp), xold(nsimp,nparam), fold(nsimp), indices(nsimp), &
        xr(nparam), xe(nparam), xc(nparam), v(nsimp,nparam),  fv(nsimp), xbar(nparam), xsol(nparam), scalev(nparam))

    ! target moments from data - college completion rates
    targ(1) = 0.0103_rk
    targ(2) = 0.0801_rk
    targ(3) = 0.3315_rk
    targ(4) = 0.0903_rk
    targ(5) = 0.0783_rk
    targ(6) = 0.1165_rk
    targ(7) = 0.2466_rk

    !write(*,*) 'college completion rate in data (1T, 2T, 3T) = (0.0619 , 0.1655, 0.3906)'
    !write(*,*) 'college completion rate in data (1Q, 2Q, 3Q, 4Q) = (0.1191 , 0.1722, 0.1826, 0.3501)'

    targ(8) = 0.1563_rk                       ! college-graduation rate
    targ(9) = 0.1554_rk                       ! average net college tuition
    targ(10) = 0.1553_rk                       ! average parental transfer
    targ(11) = 0.1800_rk                       ! cumulative student debt at the graduation CHECK!
    targ(12) = 0.40_rk                         ! graduating seniors with student debt
    targ(13) = 0.0075_rk                         ! aggregate student debt ratio

    ! weights for sum of weighted squared residuals between target and model
    weights = 0.0_rk
    do i = 1, ntarg
        weights(i) = 1.0_rk!/(targ(i))!**2.0_rk)
    end do
    !weights = weights/sum(weights)

    scale10 = 10.0_rk
    scale100 = 100.0_rk
    scale1000 = 1000.0_rk
    !initial guess for paramters
    x0 = (/                   &
        psi1,                     &     !1. psi1
        psi2,                     &     !2. psi2
        psi3,                     &     !3. psi3
        phi0,                     &     !4. phi0
        phi1,                     &     !5. phi1
        phi2,                     &     !6. phi2,
        phic0,                    &     !7. phic0
        phic1,                    &     !8. phic1
        phic2,                    &     !9. phic2
        phic3,                    &     !10. phic3
        pi,                       &     !11. adjusting parameter for earnings for college students
        chib                      &     !12
        /)

    if ( iread .eq. 1_ik) then
        OPEN(40, FILE=resultfile, STATUS='OLD', ACTION='READ')
        do i = 1, nparam
            read(40, *) x0(i)
        end do
        close(40)
    end if


    scalev(:) = 1.0_rk
    scalev(1:3) = scalev(1:3)*scale10
    scalev(4:6) =scalev(4:6)*scale100
    !scalev(2:6) =scalev(2:6)*scale10
    scalev(11) = scalev(11)*scale10
    scalev(12) = scalev(12)*scale100

    x(1,1:nparam) =x0;
    do i= 2, nsimp
        x(i, 1:nparam) = x0;
        x(i,i-1) = x0(i-1)+1.0_rk/scalev(i-1)
    end do

    !education cost
    allocate(netphi(nabil,anum), psycost(nabil, anum, ncost), probcost(nabil, anum, ncost), a0grid(nabil, anum), piea0(nabil,anum))

    do i= 1, nsimp
        write(*,'(1x, a, i2.0)')    'initial estimation i  :  ', i
        pi = x(i,11)
        chib = x(i,12)

        call educost(nabil, abilgrid, pabil, anum, nparam,   x(i,1:nparam), egridc,  pia0, piec0, tempa,  a0grid, psycost, probcost, netphi, piea0 )

        call ss_solve_est(nabil, abilgrid, pabil, eagridj, chibshock, anum, bnum, adnum, bdnum,ebar,  pop,   weights, targ, ntarg,   &
            chir,psycost, probcost, netphi,  pi,  tbar, chidage,   t, nsim,  r_b, &
            wp, wl, wh, lambda, agrid, adgrid, bgrid, bdgrid, egridc,  egridl,piec, piel, a0grid, &
            piel0, piec0, pia0, piea0, fval,chibprob)

        f(i) = fval
        write(*,'(1x, a, f8.4)')    'initial estimation f(i)  :  ', f(i)
        write(*,*)
    end do


    distance = dmax1(1.0_rk, precision_gss)

    iteration = 1_ik
    maxiters = 500_ik

    do while (distance.gt.precision_gss.and.iteration.lt.maxiters)

        !if (distance.lt.precision.or.iteration.gt.maxiters) then
        !    exit
        !end if

        xold = x; fold = f


        ! step 1, sort the simplex
        call sortme(nsimp, fold, indices)

        do i = 1, nsimp
            x(i,1:nparam) = xold(indices(i), 1:nparam); f(i) = fold(indices(i))
        end do

        ! print the existing results
        if (modulo(iteration, 1_ik).eq.0_ik ) then
            thisreport = 1_ik
            !write(*,*)
            write(*,'(1x,a,i4,a)', advance = 'no') ' continuous iteration ', iteration, ' f simplex = '
            do i = 1, nsimp
                write(*,'(1x,e12.4)', advance = 'no') f(i)
            end do
            write(*,*)
        end if

        ! Generate the centroid by averaging the nparam best points, excluding the worst
        xbar = sum(x(1:nparam,1:nparam),dim = 1); xbar = xbar/dble(nparam)

        ! reflection of worst point through centre
        xr = xbar + rhonm*(xbar - x(nsimp, 1:nparam))
        call educost(nabil, abilgrid, pabil,  anum, nparam, xr, egridc, pia0, piec0, tempa,  a0grid, psycost, probcost, netphi, piea0 )
        pi = xr(11)
        chib = xr(12)
        call ss_solve_est(nabil, abilgrid, pabil, eagridj, chibshock, anum, bnum, adnum, bdnum,ebar,  pop,   weights, targ, ntarg,   &
            chir,psycost, probcost, netphi,  pi,  tbar, chidage,   t, nsim,  r_b, &
            wp, wl, wh, lambda, agrid, adgrid, bgrid, bdgrid, egridc,  egridl,piec, piel, a0grid, &
            piel0, piec0, pia0, piea0, fval,chibprob)

        fr = fval



        ! reflection point is neither new worst nor best
        if (f(1).le.fr.and.fr.lt.f(nparam)) then
            thisstep = 0_ik
            x(nsimp,1:nparam) = xr;
            f(nsimp) = fr

            ! reflection point is new best, so attempt expansion in its direction
        elseif(fr.lt.f(1)) then
            xe = xr + chinm*(xr - xbar)
            call educost(nabil, abilgrid, pabil,  anum, nparam,  xe, egridc,  pia0, piec0, tempa,  a0grid, psycost, probcost, netphi, piea0  )
            pi = xe(11)
            chib = xe(12)
            call ss_solve_est(nabil, abilgrid, pabil, eagridj, chibshock, anum, bnum, adnum, bdnum,ebar,  pop,   weights, targ, ntarg,   &
                chir,psycost, probcost, netphi,  pi,  tbar, chidage,   t, nsim,  r_b, &
                wp, wl, wh, lambda, agrid, adgrid, bgrid, bdgrid, egridc,  egridl,piec, piel, a0grid, &
                piel0, piec0, pia0, piea0, fval,chibprob)
            fe = fval

            ! expansion point is better than reflection point, use it to replace worst point
            if(fe.lt.fr) then
                thisstep = 1_ik
                x(nsimp,1:nparam) = xe; f(nsimp) = fe

                ! expansion failed, use reflection point to replace worst point
            else
                thisstep = 2_ik
                x(nsimp,1:nparam) = xr; f(nsimp) = fr
            end if

            ! fr.ge.f(nparam): reflection point would become new worst point
        else

            shrink = 0
            ! perform an outside contraction
            if(fr.lt.f(nsimp)) then
                xc = xbar + gammanm*(xr - xbar)

                call educost(nabil, abilgrid, pabil,  anum, nparam,  xc, egridc, pia0, piec0, tempa,  a0grid, psycost, probcost, netphi, piea0 )
                pi = xc(11)
                chib = xc(12)
                call ss_solve_est(nabil, abilgrid, pabil, eagridj, chibshock, anum, bnum, adnum, bdnum,ebar,  pop,   weights, targ, ntarg,   &
                    chir,psycost, probcost, netphi,  pi,  tbar, chidage,   t, nsim,  r_b, &
                    wp, wl, wh, lambda, agrid, adgrid, bgrid, bdgrid, egridc,  egridl,piec, piel, a0grid, &
                    piel0, piec0, pia0, piea0, fval,chibprob)
                fc = fval

                if (fc.le.fr) then
                    thisstep = 3_ik
                    x(nsimp,1:nparam) = xc; f(nsimp) = fc
                else
                    shrink = 1
                end if

                ! perform an inside contraction
            else
                xc = xbar - gammanm*(xbar - x(nsimp, 1:nparam))
                call educost(nabil, abilgrid, pabil,  anum, nparam,  xc, egridc,  pia0, piec0, tempa,  a0grid, psycost, probcost, netphi, piea0  )
                pi = xc(11)
                chib = xc(12)
                call ss_solve_est(nabil, abilgrid, pabil, eagridj, chibshock, anum, bnum, adnum, bdnum,ebar,  pop,   weights, targ, ntarg,   &
                    chir,psycost, probcost, netphi,  pi,  tbar, chidage,   t, nsim,  r_b, &
                    wp, wl, wh, lambda, agrid, adgrid, bgrid, bdgrid, egridc,  egridl,piec, piel, a0grid, &
                    piel0, piec0, pia0, piea0, fval,chibprob)
                fc = fval
                if (fc.lt.f(nsimp)) then
                    thisstep = 4_ik
                    x(nsimp,1:nparam) = xc; f(nsimp) = fc
                else
                    shrink = 1
                end if

            end if

            if (shrink.eq.1) then
                thisstep = 5_ik
                v(1, 1:nparam) = x(1, 1:nparam)
                fv(1) = f(1)
                do i = 2, nsimp
                    v(i,1:nparam) = x(1,1:nparam) + sigmanm*(x(i,1:nparam) - x(1,1:nparam))
                    call educost(nabil, abilgrid, pabil, anum,  nparam,  v(i,1:nparam), egridc,  pia0, piec0, tempa,  a0grid, psycost, probcost, netphi, piea0 )
                    pi = v(i,11)
                    chib = v(i, 12)
                    call ss_solve_est(nabil, abilgrid, pabil, eagridj, chibshock, anum, bnum, adnum, bdnum,ebar,  pop,   weights, targ, ntarg,   &
                        chir,psycost, probcost, netphi,  pi,  tbar, chidage,   t, nsim,  r_b, &
                        wp, wl, wh, lambda, agrid, adgrid, bgrid, bdgrid, egridc,  egridl,piec, piel, a0grid, &
                        piel0, piec0, pia0, piea0, fval,chibprob)
                    fv(i) = fval
                end do
                x = v; f  = fv
            end if

        end if

        distance = 0.0; iteration = iteration + 1
        do i = 1, nsimp
            distance = sum((x(i,1:nparam) - xold(i,1:nparam))**2.0_rk) + distance
        end do
        distance = dsqrt(distance/(dble(nparam) + 1.0_rk))

        if (thisreport.eq.1_ik) then
            if (thisstep.eq.0_ik) then
                write(*,*) ' reflection '
            elseif (thisstep.eq.1_ik) then
                write(*,*) ' reflect with expansion '
            elseif (thisstep.eq.2_ik) then
                write(*,*) ' reflect without expansion '
            elseif (thisstep.eq.3_ik) then
                write(*,*) ' outside contraction'
            elseif (thisstep.eq.4_ik) then
                write(*,*) ' inside contraction'
            else
                write(*,*) ' shrink'
            end if
        end if



        ! thisreport = 0_ik
        open(unit = 50, file = resultfile, status = 'unknown', action = 'write')

        do i = 1,nparam
            write(50,*) (x(1,i))
        end do

        write(50,*) f(1)
        close(50)


    end do

    xsol = x(1, 1:nparam); fsol = f(1)



    end program estimation_main