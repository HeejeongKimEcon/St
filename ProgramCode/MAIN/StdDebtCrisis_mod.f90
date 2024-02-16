    module StdDebtCrisis_mod



    ! OpenMP is used in 3 places, retire, workerinitial, worker.
    ! Enable Process OpenMP Directives in Properties\ Fortran\ Language.


    ! When maximize for speed, Disable check stack frame and check array and string bounds in Properties\ Fortan\ Run-Time

    use Lib_kindset
    use Lib_conshock
    use Lib_inequality
    use Lib_sort
    use Lib_grid
    use omp_lib
    use get_params
    use lib_rwhdf5
    implicit none

    public:: workretire, expectedvalue, decisionrule, goldensection, distini, disth, distl,entiredist, defaultdecision, simulate, &
        indexcal, goldensectionb, sim_moments, bisectionlabor, simul_sample

    contains


    subroutine steadyshock(pi, nabil, grah, gral, abilgrid, expeh, expel, egridc,epgridc, etgridc, egridl, &
        piec, piel, piel0, piec0, eagridj)

    integer(ik)::  ip, ipf, it, itf, place, placef, ie, j, ia, nabil
    real(rk)::  grah, gral, sigmap(2), sigmat(2), ability, pi
    real(rk):: egridc(e_num), epgridc(epnum), etgridc(etnum), egridl(e_num), epgridl(epnum), etgridl(etnum),  &
        piepc(epnum,epnum), pietc(etnum,etnum), piec(e_num,e_num), piepl(epnum,epnum), pietl(etnum,etnum), piel(e_num,e_num), &
        piel0(e_num), piec0(e_num), expeh(jnum), expel(jnum), eagridj(jnum, nabil, e_num, 2), abilgrid(nabil)

    intent(in):: grah, gral, abilgrid, expeh, expel, nabil, pi
    intent(out):: egridc,epgridc, etgridc, egridl, piec, piel, piel0, piec0, eagridj

    !Labor productivity shock process for college-educated & non-college-educated households
    !1 is for college educated, 2 is for non-college educated
    sigmap(1) = dsqrt(varp1); sigmap(2) = dsqrt(varp2)
    sigmat(1) = dsqrt(vart1); sigmat(2) = dsqrt(vart2)

    !first, fix the support to possible widest range using the largest variances of shocks in sample periods.

    if (irouwen.eq.1_ik) then
        call rouwenhorst(rho1, sigmap(1), epnum, epgridc, piepc)
        call rouwenhorst(rho2, sigmap(2), epnum, epgridl, piepl)
        call rouwenhorst(0.0_rk, sigmat(1), etnum, etgridc, pietc)
        call rouwenhorst(0.0_rk, sigmat(2), etnum, etgridl, pietl)
    else
        call tauchen(mean, dsqrt(varp1) , rho1,  multiple, epnum, epgridc, piepc)
        call tauchen(mean, dsqrt(vart1), 0.0_rk, multiple, etnum, etgridc, pietc)
       
        call tauchen(mean, dsqrt(varp2) , rho2,  multiple, epnum, epgridl, piepl)
        call tauchen(mean, dsqrt(vart2), 0.0_rk, multiple, etnum, etgridl, pietl)        
        
        !piepl = piepc
        !pietl = pietc
        
        !call tauchenonsupport(mean, 2.0*sigmap(1), rho1, epnum, epgridc, piepc)
        !call tauchenonsupport(mean, 2.0*sigmap(2), rho2, epnum, epgridc, piepl)
        !call tauchenonsupport(mean, 2.0*sigmat(1), 0.0_rk , etnum, etgridc, pietc)
        !call tauchenonsupport(mean, 2.0*sigmat(2), 0.0_rk , etnum, etgridc, pietl)
    end if


    !Combining persistent and transitional shock
    do ip = 1, epnum
        do it = 1, etnum
            place = (ip - 1_ik)*etnum + it
            egridc(place) = epgridc(ip) + etgridc(it)
            egridl(place) = epgridl(ip) + etgridl(it) !!added
            if (irouwen.eq.1_ik) then
                egridl(place) = epgridl(ip) + etgridl(it)
            end if
            do ipf = 1, epnum
                do itf = 1, etnum
                    placef = (ipf - 1_ik)*etnum + itf
                    piec(place,placef) = piepc(ip,ipf)*pietc(it,itf)
                    piel(place,placef) = piepl(ip,ipf)*pietl(it,itf)
                end do
            end do
        end do
    end do
    do ie = 1, e_num
        piec(ie, 1:e_num) = piec(ie,1:e_num)/sum(piec(ie,1:e_num))
        piel(ie, 1:e_num) = piel(ie,1:e_num)/sum(piel(ie,1:e_num))
    end do

    egridc = exp(egridc)
    if (irouwen.eq.1_ik) then
        egridl = exp(egridl)
    else !tauchen
        egridl = exp(egridl)!egridc!
        !epgridl = epgridc
        !etgridl = etgridc
    end if
    !Ergodic distribution for income shock and parental transfer
    call  ergodicdist(e_num, piel, precerg, piel0)
    call  ergodicdist(e_num, piec, precerg, piec0)

    ! whern eagridj(:,:,1) for skilled/ eagidj(:,:,2) for unskilled
    do j = 1, jnum
        do ia = 1, nabil
            ability = abilgrid(ia)
            do ie = 1,e_num
                if ( j .lt. jc) then
                    eagridj(j,ia, ie,1) =pi*expeh(1)*exp(grah*ability)
                elseif ( j .le. jrnum) then
                    eagridj(j,ia,ie,1) = expeh(j)*egridc(ie)*exp(grah*ability)
                else
                    eagridj(j,ia,ie, 1) = egridc(ie)*(sum(expeh(jc:jrnum))/(jrnum-jc+1))*exp(grah*ability)
                end if

                if ( j .le. jrnum) then
                    eagridj(j,ia,ie,2) = expel(j)*egridl(ie)*exp(gral*ability)
                else
                    eagridj(j,ia,ie,2) = egridl(ie)*(sum(expel(jc:jrnum))/(jrnum-jc+1))*exp(gral*ability)
                end if
            end do
        end do
    end do




    end subroutine steadyshock





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!              BISECTION FOR LABOR           !!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine bisectionlabor(anext, term0, term1, term2, term3, nupper, nval, consumption)
    integer(ik):: iter
    real(rk):: anext, term0, term1, term2, term3, nval, consumption, nupper
    real(rk):: nlow, nhigh, dist2, fval, flow, fhigh

    intent(in):: anext, term0, term1, term2, term3, nupper
    intent(out)::nval, consumption

    dist2 = 2.0_rk*precision
    nlow = precision; nhigh = nupper -precision;
    nval = nlow;
    !fval =term3+term0*nval -term1*((nupper-nval)**term2)-anext;
    fval =term3+term0*nval -term1*(nval**term2)-anext;
    flow = fval;

    nval =nhigh;
    !fval =term3+term0*nval -term1*((nupper-nval)**term2)-anext;
    fval =term3+term0*nval -term1*(nval**term2)-anext;
    fhigh = fval;

    if (flow*fhigh .gt. 0.0_rk) then
        !write(*,'(1x,a)', advance = 'no') ' Could not bisect '
        !write(*,*) flow, fhigh

        consumption = 0.0!term1*((nupper-nval)**term2)
        nval = 0.0_rk!nupper-(consumption/term1)**(1.0_rk/term2)

    else
        iter =1_ik
        do while (dist2 .gt. precision)

            nval = (nlow + nhigh)/2.0_rk;
            !fval =term3+term0*nval -term1*((nupper-nval)**term2)-anext;
            fval =term3+term0*nval -term1*(nval**term2)-anext;
            if (fval*flow.gt.0) then
                nlow = nval;
                flow = fval;
            else
                nhigh = nval;
                fhigh = fval;
            end if

            dist2 = nhigh - nlow;
            iter = iter+1_ik
            if ( iter .gt. 100_ik) then
                write(*,*) 'dist2', dist2
            end if

        end do

        consumption =term3+term0*nval-anext ;
        if (consumption .le. 0.0_rk) then
            consumption =0.0!term1*((nupper-nval)**term2)
            nval = 0.0_rk!nupper-(consumption/term1)**(1.0_rk/term2)
        end if


    end if

    end subroutine bisectionlabor


    subroutine decisionrule(tbar, nabil, eagridj,chibshock, anum, bnum, chir, psycost,  piel0,piec0,  netphi, bbar,  &
        a0grid,   t, vhf_t, vlf_t, vhcollegef_t, chidage, lambda, agrid, bgrid,  wh, wl,  piec, piel, phid_abil, &
        payindex, payweight, nopayindex, nopayweight, ghworker, glworker, vhworker, vlworker,&
        pchoicec, pchoicenc, laborh, laborl, edudecision,bfchoice,bhini, vhini,  vhcollege,  laborhcollege, ghcollege)

    integer(ik):: anum, bnum, e,  t, ib,  ia, ie, j,   bindex,   ic, iab, nabil, ichib
    integer(ik):: pchoicec(jpaynum,nT,nabil, e_num,anum,bnum), pchoicenc(jpaynum,nT,nabil,e_num,anum,bnum), &
        payindex(bnum, nT), nopayindex(bnum), pchoice(nT, bnum)

    real(rk)::  wh,wl, aval, bval,  chid,  evalc, evalnc, phi,  &
        term0, term1, term2,   chir,  cval, yval, vval,  bbar, bweight,  vltemp, estat, tbar

    real(rk)::  eagridj(jnum, nabil, e_num, 2), agrid(anum), a0grid(nabil, anum), lambda(nT),bgrid(bnum), chibshock(chibnum), &
        ghworker(jnum,nT,nabil,e_num,anum,bnum), glworker(jnum,nT,nabil,e_num,anum,bnum), &
        piec(e_num,e_num), piel(e_num,e_num),  payweight(bnum, nT), nopayweight(bnum), &
        vhf_t(jnum, nT, nabil,e_num, anum, bnum), vlf_t(jnum, nT, nabil,e_num, anum, bnum), &
        laborh(jrnum,nT,nabil,e_num,anum,bnum), laborl(jrnum,nT,nabil,e_num,anum,bnum), piel0(e_num),  piec0(e_num),&
        vhini(nabil,anum,chibnum), bhini(nabil,anum, chibnum),ghini(nabil,anum, chibnum), laborhini(nabil,anum, chibnum), netphi(nabil,anum),evlini(nabil, anum),&
        edudecision( nabil,anum, ncost, chibnum),bfchoice(nabil,anum, ncost,chibnum), ev(anum, bnum, nabil), psycost(nabil,anum, ncost), chidage(jnum), labor(nT, bnum), &
        vhcollege(jc-1_ik,nabil,anum,bnum,chibnum), ghcollege(jc-1_ik,nabil,anum,bnum,chibnum), laborhcollege(jc-1_ik,nabil,anum,bnum,chibnum), &
        vhcollegef_t(jc-1_ik,nabil,anum,bnum,chibnum), vhcollegein(anum, bnum), phid_abil(nabil)

    real(rk), allocatable:: vhworker(:,:,:,:,:,:), vlworker(:,:,:,:,:,:), evh(:,:,:,:,:), evl(:,:,:,:,:)

    intent(in):: tbar, nabil, chibshock, anum, bnum, chir, a0grid,  t,  vhf_t, vlf_t, chidage, lambda, agrid, bgrid, &
        wh, wl, piec, piel, phid_abil,  payindex, payweight, nopayindex, nopayweight, &
        bbar,  netphi,  piel0, piec0,  psycost, eagridj, vhcollegef_t
    intent(out)::   ghworker, glworker, vhworker, vlworker, pchoicec, pchoicenc, laborh, laborl, edudecision, &
        bfchoice, bhini, vhini,  vhcollege, laborhcollege, ghcollege

    pchoicec= 1_ik; pchoicenc = 1_ik;
    laborh = 0.0_rk; laborl = 0.0_rk; edudecision = 0.0_rk;  bfchoice = 0.0_rk; vhcollege = 0.0_rk; ev = 0.0_rk


    allocate(vhworker(jnum,nT,nabil, e_num,anum,bnum), vlworker(jnum,nT,nabil, e_num,anum,bnum), &
        evh(nT,e_num,anum,bnum, nabil), evl(nT,e_num,anum,bnum, nabil))

    vhworker =0.0_rk;
    vlworker = 0.0_rk;


    !Solve the last age's goldensection first - we need to make sure that everyone pays back before the last age
    do ie = 1, e_num
        do iab  =1, nabil
            evalnc = eagridj(jnum,iab, ie, 2); evalc =  eagridj(jnum,iab, ie, 1);
            do ia = 1, anum
                aval = agrid(ia)
                do ib = 1, bnum
                    bval = bgrid(ib)
                    bindex = payindex(ib,nT); bweight = payweight(ib,nT)

                    !unskilled
                    call lastagesub(chir,  wl,  aval, bval,  evalnc, vval, cval)
                    vlworker(jnum, 1:nT, iab, ie, ia, ib) = vval

                    !skilled
                    call lastagesub(chir, wh,  aval, bval,  evalc, vval, cval)
                    vhworker(jnum, 1:nT, iab, ie, ia, ib) = vval

                end do
            end do
        end do
    end do


    !Solve goldensection rules backward by age
    j = jnum - 1_ik

    do while (j .ge. 1)
        chid = chidage(j)
        do iab = 1, nabil
            if (t .eq. 0_ik .or. t .eq. tnum) then
                ! for steady-state
                e = 1_ik                ! skilled workers
                if (j.ge.jc-1_ik) then
                    call expectedvalue(e, j, piec0, anum, bnum,  vhworker(j+1_ik,1:nT,iab,1:e_num,1:anum,1:bnum), piec, evh(1:nT,1:e_num,1:anum,1:bnum,iab))
                end if

                e = 0_ik                ! unskilled workers
                call expectedvalue(e, j, piec0,anum, bnum,   vlworker(j+1_ik,1:nT,iab,1:e_num,1:anum,1:bnum), piel, evl(:,:,:,:,iab))
            else
                ! for transition
                e = 1_ik
                if (j.ge.jc-1_ik) then
                    call expectedvalue(e, j, piec0, anum, bnum,   vhf_t(j+1_ik,1:nT,iab,1:e_num,1:anum,1:bnum), piec, evh(:,:,:,:,iab))
                end if

                e = 0_ik
                call expectedvalue(e, j, piec0, anum, bnum,   vlf_t(j+1_ik,1:nT,iab,1:e_num,1:anum,1:bnum), piel, evl(:,:,:,:,iab))
            end if
        end do

        ! When students drop out, they redraw their productivity from unskilled-workers' shock processes
        if ( j .eq. jd-1_ik) then! jc-1_ik) then
            ev = 0.0_rk
            do ia = 1,anum
                do ib = 1,bnum
                    do iab= 1, nabil
                        ev(ia, ib, iab) = dot_product(piel0(1:e_num),evl(1,1:e_num, ia,ib, iab))
                    end do
                end do
            end do
        end if


        do iab = 1, nabil
            !$omp parallel do private(ia, ie, evalc, evalnc, aval, phi, yval, term0, term1, labor, pchoice,estat)
            do ie = 1, e_num
                if (j.ge.jc) then
                    evalnc = eagridj(j, iab, ie, 2); evalc =  eagridj(j,iab, ie, 1);
                else !j is between 1 and 4
                    evalc = eagridj(j,iab, ie, 1)
                    evalnc = eagridj(j,iab, ie, 2)
                end if

                do ia = 1, anum
                    if (j.ge.jc) then ! j >= 5. solve for both skilled and unskilled workers
                        aval = agrid(ia)
                        estat = 0.0_rk !work status (i.e. 1.0_rk = during college/ 0.0_rk = after college or unskilled)
                        !skilled
                        call workretire(0.0_rk, anum, bnum, chir, estat,  chid,  j,   aval, wh,  evalc,  agrid, &
                            evh(1:nT,ie,1:anum,1:bnum, iab), lambda, bgrid,  payindex, payweight, nopayindex, nopayweight, pchoice, &
                            vhworker(j,1:nT,iab,ie,ia,1:bnum), ghworker(j,1:nT,iab,ie,ia,1:bnum), labor)

                        if ( j .le. jrnum) then
                            laborh(j,1:nT,iab,ie,ia,1:bnum) = labor(1:nT, 1:bnum)
                        end if

                        if( j .le. jpaynum) then
                            pchoicec(j,1:nT,iab,ie,ia,1:bnum) = pchoice(1:nT, 1:bnum)
                        end if
                        estat = 0.5_rk
                        !unskilled
                        call workretire(0.0_rk,  anum, bnum, chir, estat,  chid,  j,   aval,  wl,  evalnc,  agrid, &
                            evl(1:nT,ie,1:anum,1:bnum, iab), lambda, bgrid, payindex, payweight, nopayindex, nopayweight, pchoice, &
                            vlworker(j,1:nT,iab,ie,ia,1:bnum), glworker(j,1:nT,iab,ie,ia,1:bnum), labor)
                        if ( j .le. jrnum) then
                            laborl(j,1:nT,iab,ie,ia,1:bnum) = labor(1:nT, 1:bnum)
                        end if
                        if( j .le. jpaynum) then
                            pchoicenc(j,1:nT,iab,ie,ia,1:bnum) = pchoice(1:nT, 1:bnum)
                        end if

                    elseif (j.lt.jc .and. j.gt.1_ik) then ! j = 2, 3, 4
                        aval = agrid(ia)
                        !unskilled
                        estat = 0.5_rk
                        call workretire(0.0_rk, anum, bnum, chir, estat,  chid, j, aval,  wl,  evalnc,  agrid, &
                            evl(1:nT,ie,1:anum,1:bnum,iab), lambda, bgrid,  payindex, payweight, nopayindex, nopayweight, pchoice, &
                            vlworker(j,1:nT,iab,ie,ia,1:bnum), glworker(j,1:nT,iab,ie,ia,1:bnum), labor)

                        laborl(j,1:nT,iab,ie,ia,1:bnum) = labor(1:nT, 1:bnum)
                        if( j .le. jpaynum) then
                            pchoicenc(j,1:nT,iab,ie,ia,1:bnum) = pchoice(1:nT, 1:bnum)
                        end if

                        !skilled
                        estat = 1.0_rk
                        aval = a0grid(iab, ia);
                        phi = netphi(iab,ia)

                        if ( ie .eq. 1_ik) then
                            do ichib = 1, chibnum
                                if (j.eq.jc-1_ik) then
                                    vhcollegein = evh(1,ie,1:anum, 1:bnum, iab)
                                else
                                    if (t .eq. 0_ik .or. t .eq. tnum) then
                                        vhcollegein = vhcollege(j+1_ik,iab, 1:anum,1:bnum,ichib)
                                    else
                                        vhcollegein = vhcollegef_t(j+1_ik,iab, 1:anum,1:bnum,ichib)
                                    end if
                                end if

                                call college_student(tbar, estat, chibshock(ichib), anum, bnum, evalc, phi, phid_abil(iab),  j,  ia, aval, bgrid, agrid, wl,  &
                                    vhcollegein, ev(:,:, iab),ghcollege(j,iab,ia,1:bnum,ichib), vhcollege(j,iab,ia,1:bnum,ichib) , &
                                    laborhcollege(j,iab,ia,1:bnum,ichib))!
                            end do
                        end if


                    else !j = 1.
                        aval = a0grid(iab, ia);

                        !non-college
                        estat = 0.5_rk
                        yval = aval + (1-tau)*wl*evalnc*nworker
                        term0 = (1-tau)*wl*evalnc; term1 = (psi/term0)**(-1/sigma); term2=-eta/sigma;
                        call initialage_nc(estat, anum, bnum,  j,  yval, aval, agrid, term0, term1, term2, evl(1:nT,ie, 1:anum,1:bnum, iab), &
                            glworker(j,1:nT,iab,ie,ia,1:bnum), vlworker(j,1:nT,iab,ie,ia,1:bnum),  labor)
                        laborl(j,1:nT,iab,ie,ia,1:bnum) = labor(1:nT, 1:bnum)

                        !college
                        if ( ie .eq. 1_ik) then
                            aval = a0grid(iab,ia)
                            term0 = (1-tau)*wl*evalc;  term1 = (psi_c/term0)**(-1/sigma); term2=-eta/sigma; phi = netphi(iab,ia)

                            do ichib = 1, chibnum

                                if (t .eq. 0_ik .or. t .eq. tnum) then
                                    vhcollegein = vhcollege(j+1_ik,iab, 1:anum,1:bnum,ichib)
                                else
                                    vhcollegein = vhcollegef_t(j+1_ik,iab, 1:anum,1:bnum,ichib)
                                end if

                                call initialage_c(tbar, chibshock(ichib), anum, bnum, ia,  bbar,  bgrid, term0, term1, term2, aval, phi,   &
                                    vhcollegein, bhini(iab,ia,ichib), ghini(iab,ia,ichib), laborhini(iab,ia,ichib), &
                                    vhini(iab,ia,ichib))
                                laborhcollege(1,iab,ia,1:bnum,ichib) = laborhini(iab,ia, ichib)
                                ghcollege(1, iab, ia, 1:bnum, ichib) = ghini(iab,ia, ichib)
                            end do


                        end if
                    end if

                end do

            end do

            !$omp end parallel do
        end do

        j = j - 1_ik

    end do

    edudecision =0.0_rk; bfchoice = 0.0_rk;

    evlini = 0
    do iab = 1, nabil
        do ia = 1,anum
            evlini(iab, ia) =  dot_product(vlworker(1,1,iab, 1:e_num,ia,bnum),piel0(1:e_num))
        end do
    end do


    do ic = 1,ncost
        do iab = 1, nabil
            do ia = 1, anum
                if (iergovlini.eq.1_ik) then
                    vltemp = vlworker(1,1,iab,1,ia,bnum)
                else
                    vltemp = evlini(iab, ia)
                end if
                do ichib = 1, chibnum
                    if (vhini(iab,ia,ichib) - psycost(iab,ia,ic).ge. vltemp) then !If college is better
                        edudecision(iab,ia,ic,ichib) = 1.0_rk
                        bfchoice(iab,ia,ic,ichib) = bhini(iab,ia,ichib)
                    else !If work is better
                        edudecision(iab,ia,ic,ichib) = 0.0_rk
                        bfchoice(iab,ia,ic,ichib) = 0.0_rk
                    end if
                end do
            end do
        end do
    end do


    end subroutine decisionrule


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     decision rule during the college     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine college_student(tbar, estat, chib, anum, bnum,   evalc, phi, phid, j, ia,  aval, bgrid,  agrid,  wl,   &
        evh, evl,  ghworker, vhworker, laborh)

    integer(ik):: anum, bnum, j,ia,bindex,ib

    real(rk)::  evalc, phi, phid, aval, bval,   wl,    yval, term0, term1, term2, term3, alow, ahigh, &
        bweight,objective,afval,cval,nval,chib,estat, tbar

    real(rk)::   bgrid(bnum), agrid(anum), evh(anum,bnum),evl(anum,bnum), ghworker(bnum), vhworker(bnum), &
        laborh(bnum),ev(anum,bnum), evhtemp(anum,2)

    intent(in)::  chib, anum, bnum,  evalc, phi, j, ia,aval, bgrid,  agrid, wl, evh,evl,estat, tbar, phid
    intent(out):: ghworker,vhworker,laborh

    !Now solve college from j = 4 to j = 2 backwards
    ghworker=0.0_rk; vhworker=0.0_rk;  laborh =0.0_rk

    if (j.eq. jd-1_ik) then!(jc-1)) then
        ev(1:anum,1:bnum) = phid*evl(1:anum,1:bnum) +(1.0_rk-phid)*evh(1:anum,1:bnum)
    else
        ev(1:anum,1:bnum) = evh(1:anum,1:bnum)
    end if

    ! !$omp parallel do private(ib,bval, yval, term3, ahigh, evhtemp, bindex, bweight, afval, objective, cval, nval )
    do ib = 1,bnum
        bval = bgrid(ib);
        term0 = (1-tau)*wl*evalc; term1 = (psi_c/term0)**(-1/sigma); term2=-eta/sigma;

        yval = aval + (1-tau)*wl*evalc*(ncollege-tbar) - dmax1(phi + bval/(dble(jc)-1.0_rk), 0.0_rk)+chib*bval !wage is low during education period.
        term3 = aval -dmax1(phi + bval/(dble(jc)-1.0_rk), 0.0_rk)+chib*bval

        if (j .eq. jc-1_ik) then

            if (yval .gt. 0.0_rk .and. bval .ge. dmax1(bbar, -fam*phi*(dble(jc)-1.0_rk))) then
                alow = agrid(1); ahigh = dmin1(agrid(anum), yval - precision)
                if (ahigh .gt. alow) then
                    if (ib .lt. bnum) then
                        evhtemp(1:anum,1:2) = ev(1:anum,ib:ib+1)
                        bindex = 1; bweight = 1.0_rk
                    else
                        evhtemp(1:anum,1:2) = ev(1:anum,ib-1:ib)
                        bindex = 1; bweight = 0.0_rk
                    end if
                    call goldensection(tbar, estat , term0, term1, term2, term3, j, bindex, bweight,  alow, ahigh, anum, 2_ik,  agrid, &
                        evhtemp, yval, objective, afval, cval,nval)
                    ghworker(ib) = afval
                    if (bval .lt. 0.0_rk) then
                        vhworker(ib) = objective!-chib
                    else
                        vhworker(ib) = objective
                    end if
                    laborh(ib)  = nval
                else
                    vhworker(ib)  = -9999999_rk
                end if
            else
                vhworker(ib)  = -9999999_rk
            end if

        else !j = 2,3
            !
            !nval = ncollege
            !yval=term3+term0*nval;
            call bisectionlabor(0.0_rk, term0, term1, term2, term3, 1.0_rk-tbar, nval, yval)

            if (yval .gt. 0.0_rk  .and. bval .ge. dmax1(bbar, -fam*phi*(dble(jc)-1.0_rk))) then
                alow = agrid(1); ahigh = dmin1(agrid(anum), yval - precision)
                ghworker(ib) =0.0_rk
                laborh(ib) = nval
                vhworker(ib) = (yval**(1.0_rk-sigma))/(1.0_rk-sigma) - psi_c*(((nval)**(1.0_rk+eta))/(1+eta)) + beta*ev(ia,ib)
                if (bval .lt. 0.0_rk) then
                    vhworker(ib) = vhworker(ib)!-chib
                else
                    vhworker(ib) = vhworker(ib)
                end if
            else
                vhworker(ib)  = -9999999_rk
            end if

        end if
    end do
    ! !$omp end parallel do

    end subroutine college_student

    subroutine lastagesub( chir,  w, aval, bval, eval, vval, cval)
    real(rk)::  chir,  w
    real(rk):: aval, bval, zeta1, yval, vval,  cval,  eval
    intent(in):: chir, w, aval, bval,   eval
    intent(out):: vval, cval

    if (bval.lt.0.0_rk) then
        zeta1 = chir
    else
        zeta1 = 0.0_rk
    end if

    !repays the outstanding student debt balance at last period.
    yval  = (1.0_rk+r)*aval + (1.0_rk-tau)*w*ssrr*eval + bval
    call lastage(yval, cval, vval)
    vval = vval-zeta1


    end subroutine lastagesub

    subroutine lastage(yval, cval, vval)
    real(rk)::yval, cval, vval
    intent(in):: yval
    intent(out):: cval, vval

    if (yval .gt. 0.0_rk) then
        vval = (yval**(1.0_rk-sigma))/(1.0_rk-sigma)
        !vval = (yval**(1.0_rk-sigma))/(1.0_rk-sigma)+psi*(((1.0_rk)**(1.0_rk-eta))/(1-eta)) !optimal choice at the last period is to spend everything.
        cval = yval
    else
        vval =  -9999999.0_rk
    end if
    end subroutine lastage

    subroutine initialage_nc(estat,  anum, bnum,  j,  yval, aval, agrid, term0, term1, term2,  &
        evl,  gl, vl, laborl)

    integer(ik)::anum, bnum,  ip, ib,  j

    real(rk):: alow, ahigh, yval,  agrid(anum), term0, term1, term2, term3,  evl(nT, anum, bnum),aval,  &
        objective, afval, cval, nval, gl(nT,bnum), vl(nT,bnum),laborl(nT,bnum), estat

    intent(in):: anum, bnum, j,  yval, agrid, term0, term1, term2,  evl,  aval, estat
    intent(out)::gl, vl, laborl

    gl=0.0_rk; vl=0.0_rk; laborl=0.0_rk;

    term3 = aval
    do ip = 1,nT
        do ib = 1,bnum
            if (yval .gt. 0.0_rk) then
                alow = agrid(1); ahigh = dmin1(agrid(anum), yval - precision)

                if (ahigh.gt.alow) then
                    call goldensection(0.0_rk, estat , term0, term1, term2, term3, j,  bnum-1_ik, 0.0_rk,  alow, ahigh, anum, bnum, agrid, evl(ip, 1:anum, 1:bnum), yval,  objective, afval, cval, nval)
                    gl(ip,ib) = afval
                    vl(ip,ib) = objective
                    laborl(ip,ib) = nval
                else
                    vl(ip,ib) = -9999999_rk
                end if

            else
                vl(ip,ib) = -9999999_rk
            end if
        end do
    end do
    end subroutine initialage_nc

    subroutine initialage_c(tbar, chib, anum, bnum, ia,  bbar, bgrid, term0, term1, term2, aval, phi,    evhjc, &
        bhini, ghini, laborhini, vhini)

    integer(ik)::   anum, bnum,   ia
    real(rk)::     bgrid(bnum), bval, term0, term1, term2, term3, aval, phi,   nval, yval,  &
        evhjc(anum,bnum), objective, bhini, ghini, laborhini, vhini, bbar, blow, bhigh, chib, tbar, cval

    intent(in)::    chib, anum, bnum, ia, bgrid, term0, term1, term2, aval, phi, tbar, evhjc, bbar
    intent(out)::   bhini, ghini, laborhini, vhini

    term3 =aval
    blow = dmax1(bbar, -fam*phi*(dble(jc)-1.0_rk))
    bhigh = 0.0_rk
    yval = term0*(ncollege-tbar) + term3
    if ( blow .le. bhigh .and. yval .gt. 0.0_rk) then
        call goldensectionb(tbar,chib, phi, term0, term1, term2, term3,  blow, bhigh, bnum, bgrid, &
            evhjc(ia,1:bnum), objective, bval, cval,nval)
    else
        objective = -9999999_rk
    end if

    if (objective .gt. -9999999_rk) then
        bhini = bval
        ghini  = 0.0_rk
        laborhini = nval
        vhini = objective;
    else
        laborhini = 0.0_rk
        bhini = 0.0_rk
        ghini = 0.0_rk
        vhini = -9999999_rk
    end if
    end subroutine initialage_c

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Subroutine: goldensection rules for workers and retirees !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine workretire(tbar, anum, bnum, chir, estat, chid,  j,  aval, wage,  eval, agrid,  &
        ev, lambda, bgrid, payindex, payweight, nopayindex,  nopayweight, pchoice, vworker, gworker, labor)

    integer(ik):: anum, bnum, bindex, ib, ip, bindexnp, j
    integer(ik):: pchoice(nT,bnum),  payindex(bnum, nT), nopayindex(bnum)

    real(rk):: aval, wage, eval, agrid(anum),  bweight,  bval, yvalp, yvalnp,  chid, &
        bweightnp, term0, term1, term2, term3, term3np,estat, chir, zeta1,tbar, term0np, term1np
    real(rk):: ev(nT,anum,bnum), vworker(nT,bnum), gworker(nT,bnum),  labor(nT,bnum), lambda(nT),bgrid(bnum),  &
        payweight(bnum, nT), nopayweight(bnum)

    intent(in):: anum, bnum, chir, estat,   j,   aval, wage,  eval,  agrid, &
        ev, lambda, bgrid, payindex, nopayindex, payweight, nopayweight, chid,tbar
    intent(out):: pchoice, vworker, gworker, labor

    pchoice = 0_ik; vworker = -9999999_rk; gworker =0.0_rk;

    term0 = (1-tau)*wage*eval; term0np = (1-tau)*(1-wgar)*wage*eval
    term1 = (psi/term0)**(-1/sigma); term1np = (psi/term0np)**(-1/sigma); 
    term2=-eta/sigma; term3np = (1+r)*aval
    !This module is used for both college-graduates (from j=5 to 64) and non-college graduates(from j = 1 to 64)
    do ip = 1,nT
        do ib = 1,bnum

            bval = bgrid(ib)
            if (bval.lt.0.0_rk) then
                zeta1 = chir
            else
                zeta1 = 0.0_rk
            end if

            !dont subsidize it after college education periods
            bindex = payindex(ib,ip); bweight = payweight(ib,ip)

            !all retirees pay according to the schedules.
            if (j .gt. jpaynum) then
                if ( j .le. jrnum) then
                    if (ilaborfix.eq.1_ik) then
                        if (estat .eq. 0.0_rk) then
                            yvalp  = (1.0_rk+r)*aval + (1.0_rk-tau)*wage*eval*nskilled + lambda(ip)*bval !debt payment
                        elseif (estat .eq. 0.5_rk) then
                            yvalp  = (1.0_rk+r)*aval + (1.0_rk-tau)*wage*eval*nunskilled + lambda(ip)*bval !debt payment
                        end if
                    else
                            yvalp  = (1.0_rk+r)*aval + (1.0_rk-tau)*wage*eval*nworker + lambda(ip)*bval !debt payment    
                    end if
                else
                    yvalp = (1.0_rk+r)*aval + (1.0_rk-tau)*wage*ssrr*eval + lambda(ip)*bval
                end if
                yvalnp = yvalp
                bindexnp = bindex; bweightnp = bweight
            else
                if (ilaborfix.eq.1_ik) then
                    if (estat .eq. 0.0_rk) then
                        yvalp  = (1.0_rk+r)*aval + (1.0_rk-tau)*wage*eval*nskilled + lambda(ip)*bval !debt payment
                        yvalnp = (1.0_rk+r)*aval + (1.0_rk-tau)*wage*eval*nskilled !if delinquent
                    elseif (estat.eq. 0.5_rk) then
                        yvalp  = (1.0_rk+r)*aval + (1.0_rk-tau)*wage*eval*nunskilled + lambda(ip)*bval !debt payment
                        yvalnp = (1.0_rk+r)*aval + (1-tau)*(1-wgar)*wage*eval*nunskilled !if delinquent
                    end if
                else
                    yvalp  = (1.0_rk+r)*aval + (1.0_rk-tau)*wage*eval*nworker + lambda(ip)*bval !debt payment
                    yvalnp = (1.0_rk+r)*aval + (1.0_rk-tau)*(1-wgar)*wage*eval*nworker !if delinquent                    
                end if
                bindexnp = nopayindex(ib); bweightnp = nopayweight(ib)
            end if

            term3 = (1.0_rk+r)*aval + lambda(ip)*bval

            call defaultdecision(tbar, anum, bnum, zeta1, term3, term3np, yvalp, yvalnp, term0, term0np, term1np, term1, term2, ip, estat,  &
                j,  agrid,   ev(1:nT,1:anum,1:bnum),  chid,bweight, bweightnp,bindex, bindexnp, &
                pchoice(ip,ib), vworker(ip,ib), gworker(ip,ib), labor(ip,ib) )

        end do
    end do

    end subroutine workretire


    subroutine defaultdecision(tbar, anum, bnum,  chir, term3, term3np,  yvalp, yvalnp, term0, term0np, term1np, term1, term2,ip, estat, j,   &
        agrid,   ev,  chid,bweight, bweightnp,bindex, bindexnp, pchoice, vworker, gworker, labor )

    integer(ik)::  anum, bnum, j,   pchoice,  bindex, ip, ipf, bindexnp

    real(rk):: agrid(anum),   ev(nT, anum,bnum), vworker, gworker,  labor,tbar
    real(rk)::   bweight, alow, yvalp, yvalnp, objectivep, objectivenp, chid, afvalp, afvalnp, cvalp, cvalnp, &
        ahighp, ahighnp, bweightnp, term0,term0np, term1, term1np, term2, term3, term3np, estat, nvalp, nvalnp, chir

    intent(in)::  anum, bnum,  chir, estat,   j,    agrid,   ev, tbar,&
        chid,bweight, bweightnp,bindex, bindexnp, ip, yvalp, yvalnp, term0,term0np, term1,term1np, term2,term3, term3np
    intent(out):: pchoice, vworker, gworker, labor

    pchoice=0_ik; gworker=0.0_rk; labor=0.0_rk; afvalp = 0.0_rk; cvalp =0.0_rk;



    if (idefault.eq.1_ik) then
        if (yvalp .gt. 0.0_rk) then

            !solve when it pays!
            alow = agrid(1); ahighp = dmin1(agrid(anum), yvalp - precision)

            ipf = min(ip+1,nT)
            if (ahighp .gt. alow) then
                call goldensection(tbar, estat, term0, term1, term2, term3, j, bindex, bweight, &
                    alow, ahighp, anum, bnum,  agrid, ev(ipf, 1:anum, 1:bnum), yvalp, objectivep, afvalp, cvalp, nvalp)
            else
                objectivep = -9999999
            end if
            !retirement - always pays
            if (j .gt. jrnum ) then
                gworker = afvalp
                vworker = objectivep - chir
                pchoice = 1_ik
                labor = 0.0_rk
            else !working age
                gworker = afvalp
                vworker = objectivep
                pchoice = 1_ik
                labor = nvalp
            end if

        else
            vworker = -9999999_rk
        end if
    else
        if (yvalp .gt. 0.0_rk) then

            !solve when it pays!
            alow = agrid(1); ahighp = dmin1(agrid(anum), yvalp - precision)

            ipf = min(ip+1,nT)
            if (ahighp .gt. alow) then
                call goldensection(tbar, estat, term0, term1, term2, term3, j, bindex, bweight, &
                    alow, ahighp, anum, bnum,  agrid, ev(ipf, 1:anum, 1:bnum), yvalp, objectivep, afvalp, cvalp, nvalp)
            else
                objectivep = -9999999
                pchoice = 0_ik
            end if

            !retirement - always pays
            if (j .gt. jrnum ) then
                gworker = afvalp
                vworker = objectivep - chir
                pchoice = 1_ik
                labor = 0.0_rk
            else !working age
                if (j .gt. jpaynum) then !pay if age has reached the force payment age - always pays
                    gworker = afvalp
                    vworker = objectivep
                    pchoice = 1_ik
                    labor = nvalp
                else
                    !have option not to pay - solve when it doesnt pay
                    if (yvalnp .gt. 0.0_rk) then
                        alow = agrid(1); ahighnp = dmin1(agrid(anum), yvalnp - precision) !this one should change
                        if (ahighnp .gt. alow) then
                            call goldensection( tbar, estat, term0np, term1np, term2, term3np, j, bindexnp, bweightnp, &
                                alow, ahighnp, anum, bnum, agrid, ev(ip,1:anum, 1:bnum), yvalnp, objectivenp, afvalnp, cvalnp, nvalnp)
                        else
                            objectivenp = -9999999
                        end if
                    else
                        objectivenp = -9999999
                    end if

                    if (objectivenp .ne. -9999999_rk) then
                        if (objectivep .gt. objectivenp-chid) then
                            gworker = afvalp
                            vworker = objectivep
                            pchoice = 1_ik
                            labor = nvalp
                        else
                            gworker = afvalnp
                            vworker= objectivenp - chid
                            labor = nvalnp
                            pchoice = 0_ik
                        end if
                    else
                        vworker = -9999999
                    end if
                end if
            end if

        elseif (yvalnp .gt. 0.0_rk) then  !here, remember that yvalp is zero, so only option available is not to pay
            alow = agrid(1); ahighnp = dmin1(agrid(anum), yvalnp - precision)
            pchoice = 0_ik
            if (ahighnp .gt. alow) then

                call goldensection( tbar, estat, term0np, term1np, term2, term3np, j, bindexnp, bweightnp, &
                    alow, ahighnp, anum, bnum, agrid, ev(ip,1:anum, 1:bnum), yvalnp, objectivenp, afvalnp, cvalnp,nvalnp)
                gworker = afvalnp

                if (j .gt. jrnum) then
                    vworker = objectivenp - chir !since we made yvalnp as paying the debt after retirement.
                    labor = 0.0_rk
                elseif (j.le.jrnum .and. j.gt.jpaynum) then
                    vworker = objectivenp -chid
                    labor = nvalnp
                else
                    vworker = objectivenp - chid
                    labor = nvalnp
                end if
            else
                vworker = -9999999

            end if

        else
            vworker = -9999999_rk
        end if
    end if

    end subroutine defaultdecision
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Subroutine: calculate expected value function   !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine expectedvalue(e, j, piec0, anum, bnum,   vworker, pie, ev)
    integer(ik)::  anum, bnum, j
    integer(ik):: ip, ie, ia, ib, e
    real(rk):: vbal
    real(rk):: vworker(nT,e_num,anum,bnum), vidx(nT,e_num,anum,bnum), piec0(e_num),  &
        pietemp(e_num), vtemp(e_num), pie(e_num,e_num), pietemp2(e_num), ev(nT,e_num,anum,bnum)

    intent(in)::  anum, bnum, e,  vworker, pie, j, piec0
    intent(out):: ev

    ev = 0.0_rk
    vidx = 0.0_rk; vtemp =0.0_rk;
    do ie =1,e_num
        do ia = 1,anum
            do ib =1,bnum
                do ip = 1,nT
                    vbal = vworker(ip, ie, ia, ib)
                    if (vbal .le. -9999999_rk) then
                        vidx(ip,ie,ia,ib) = 0.0_rk
                    else
                        vidx(ip,ie,ia,ib) = 1.0_rk
                    end if
                end do
            end do
        end do
    end do


    do ie = 1, e_num
        if (j.gt.jrnum) then  !Retirees do not experience income shock.
            pietemp(1:e_num) = 0.0_rk
            pietemp(ie) = 1.0_rk
        else !Income shock applies to workers only! Also to students (i.e.dropout)
            if (e.eq.0_ik) then
                pietemp(1:e_num) = 0.0_rk
                pietemp(1:e_num) = pie(ie,1:e_num)
            elseif (e .eq. 1_ik .and. j .lt. jc-1_ik) then
                pietemp(1:e_num) = 0.0_rk
                pietemp(ie) = 1.0_rk
            elseif (e .eq. 1_ik .and. j .eq. jc-1_ik) then
                pietemp = piec0
            else
                pietemp(1:e_num) = 0.0_rk
                pietemp(1:e_num) = pie(ie,1:e_num)
            end if
        end if

        do ia = 1, anum
            do ib = 1, bnum
                do ip = 1,nT
                    vtemp(1:e_num) = vworker(ip, 1:e_num, ia, ib)  !Choosing last period
                    pietemp2(1:e_num) = vidx(ip,1:e_num,ia,ib)*pietemp(1:e_num) !Giving zero weight to those with -9999999 values
                    ev(ip,ie,ia,ib) = dot_product(pietemp2,vtemp)/sum(pietemp2) !expected conditional value function when income shock transition is pie
                    if (isnan(ev(ip,ie,ia,ib))) then
                        ev(ip,ie,ia,ib) =-9999999_rk
                    end if
                end do
            end do
        end do
    end do

    end subroutine expectedvalue




    subroutine goldensection(tbar, estat, term0, term1, term2, term3, j,  bindex, bweight,  &
        alow, ahigh, anum, bnum, agrid, ev, yval,  objective, afval, cval,nval)

    integer(ik)::    anum, aindex, bindex, bnum, j
    real(rk)::  a, alow, b, yval,  ahigh, rg, c, agrid(anum), aweight, evalf, ev(anum, bnum), fval, tbar, &
        fc, d, s1, fd, z,  objective, aval, afval, cval, bweight, term0, term1, term2, term3, estat,consumption,nval, nupper

    intent(in):: estat,term0, term1, term2, term3, j,   alow, ahigh, anum, bnum, agrid, &
        ev, yval,  bindex, bweight,tbar
    intent(out):: objective, afval, cval, nval


    !First get optimal labor and corresponding income(i.e. yval)

    a = alow; b = (yval - precision*10.0_rk); b = dmin1(b,ahigh)
    rg = (3.0_rk - dsqrt(5.0_rk))/2.0_rk

    !First c evaluation
    c = a + rg*(b-a)

    aval = c

    if ( j.le. jrnum .and. j.ge.jc) then !working
        nupper = 1.0_rk
        if (ilaborfix.eq.1_ik) then
            if (estat.eq.0.0_rk) then !skilled
                nval = nskilled
            elseif (estat.eq. 0.5_rk) then !unskilled
                nval = nunskilled
            end if
            consumption =term3+term0*nval-aval
        else
            call bisectionlabor( aval, term0, term1, term2, term3, nupper, nval, consumption)
        end if
    elseif (j.lt.jc .and. j.ge.1_ik) then !studying
        nupper = 1.0_rk -tbar
        call bisectionlabor( aval, term0, term1, term2, term3, nupper, nval, consumption)
    else
        nupper = 1.0_rk
        nval = 0.0_rk
        consumption = yval-aval
    end if


    if (consumption .gt. 0.0_rk) then
        call gridlookup(anum, agrid, aval, aindex, aweight)
        if (ev(aindex,bindex) .le. -9999999_rk) then
            aweight = 0.0_rk
        end if
        evalf = bweight*(ev(aindex, bindex)*aweight + ev(aindex+1_ik,bindex)*(1.0_rk - aweight)) +(1.0_rk-bweight)*(ev(aindex, bindex+1)*aweight + ev(aindex+1_ik,bindex+1)*(1.0_rk - aweight))

        if (j.lt.jc .and. estat.eq.1.0_rk) then
            fval = ((consumption)**(1.0_rk - sigma))/(1.0_rk - sigma)-psi_c*((nval)**(1.0_rk+eta)/(1+eta)) + beta*evalf
        else
            fval = ((consumption)**(1.0_rk - sigma))/(1.0_rk - sigma)-psi*((nval)**(1.0_rk+eta)/(1+eta)) + beta*evalf
        end if
    else
        fval = -9999999_rk
        !write(*,*) 'error1'
    end if

    fc = -fval

    !Do it again for d

    d = a + (1-rg)*(b-a)
    aval = d

    if ( j.le. jrnum .and. j.ge.jc) then !working
        nupper = 1.0_rk
        if (ilaborfix.eq.1_ik) then
            if (estat.eq.0.0_rk) then !skilled
                nval = nskilled
            elseif (estat.eq. 0.5_rk) then !unskilled
                nval = nunskilled
            end if
            consumption =term3+term0*nval-aval
        else
            call bisectionlabor( aval, term0, term1, term2, term3, nupper, nval, consumption)
        end if
    elseif (j.lt.jc .and. j.ge.1_ik) then !studying
        nupper = 1.0_rk - tbar
        call bisectionlabor( aval, term0, term1, term2, term3, nupper, nval, consumption)

    else
        nupper = 1.0_rk
        nval = 0.0_rk
        consumption = yval-aval
    end if


    if (consumption .gt. 0.0_rk) then
        call gridlookup(anum, agrid, aval, aindex, aweight)
        if (ev(aindex,bindex) .le. -9999999_rk) then
            aweight = 0.0_rk
        end if

        evalf = bweight*(ev(aindex, bindex)*aweight + ev(aindex+1_ik,bindex)*(1.0_rk - aweight)) +(1.0_rk-bweight)*(ev(aindex, bindex+1)*aweight + ev(aindex+1_ik,bindex+1)*(1.0_rk - aweight))
        if (j.lt.jc .and. estat.eq.1.0_rk) then
            fval = ((consumption)**(1.0_rk - sigma))/(1.0_rk - sigma)-psi_c*((nval)**(1.0_rk+eta)/(1+eta)) + beta*evalf
        else
            fval = ((consumption)**(1.0_rk - sigma))/(1.0_rk - sigma)-psi*((nval)**(1.0_rk+eta)/(1+eta)) + beta*evalf
        end if
    else
        fval = -9999999_rk
        ! write(*,*) 'error2'
    end if
    fd = -fval

    !npdf
    s1 = 0;

    !As long as c and d differ, resample an intermediate value z

    do while ((d-c) > precision_gss)
        s1 = s1 + 1
        if (fc.gt.fd) then
            z = c + (1.0_rk-rg)*(b-c)
            !Case 1: [a c d b ] <--- [c d z b]
            a = c
            c = d
            fc = fd
            d = z

            aval = d

            if ( j.le. jrnum .and. j.ge.jc) then !working
                nupper = 1.0_rk
                if (ilaborfix.eq.1_ik) then
                    if (estat.eq.0.0_rk) then !skilled
                        nval = nskilled
                    elseif (estat.eq. 0.5_rk) then !unskilled
                        nval = nunskilled
                    end if
                    consumption =term3+term0*nval-aval
                else
                    call bisectionlabor( aval, term0, term1, term2, term3, nupper, nval, consumption)
                end if
            elseif (j.lt.jc .and. j.ge.1_ik) then !studying
                nupper = 1.0_rk -tbar
                call bisectionlabor( aval, term0, term1, term2, term3, nupper, nval, consumption)

            else
                nupper = 1.0_rk
                nval = 0.0_rk
                consumption = yval-aval
            end if


            if (consumption .gt. 0.0_rk) then

                call gridlookup(anum, agrid, aval, aindex, aweight)
                if (ev(aindex,bindex) .le. -9999999_rk) then
                    aweight = 0.0_rk
                end if
                evalf = bweight*(ev(aindex, bindex)*aweight + ev(aindex+1_ik,bindex)*(1.0_rk - aweight)) +(1.0_rk-bweight)*(ev(aindex, bindex+1)*aweight + ev(aindex+1_ik,bindex+1)*(1.0_rk - aweight))

                if (j.lt.jc .and. estat.eq.1.0_rk) then
                    fval = ((consumption)**(1.0_rk - sigma))/(1.0_rk - sigma)-psi_c*((nval)**(1.0_rk+eta)/(1+eta)) + beta*evalf
                else
                    fval = ((consumption)**(1.0_rk - sigma))/(1.0_rk - sigma)-psi*((nval)**(1.0_rk+eta)/(1+eta)) + beta*evalf
                end if
            else
                fval = -9999999_rk
            end if
            fd = -fval
        else
            z = a + rg*(d-a)
            !Case 2: [a c d b] <--- [a z c d]
            b = d
            d = c
            fd = fc
            c = z

            aval = c


            if ( j.le. jrnum .and. j.ge.jc) then !working
                nupper = 1.0_rk
                if (ilaborfix.eq.1_ik) then
                    if (estat.eq.0.0_rk) then !skilled
                        nval = nskilled
                    elseif (estat.eq. 0.5_rk) then !unskilled
                        nval = nunskilled
                    end if
                    consumption =term3+term0*nval-aval
                else
                    call bisectionlabor( aval, term0, term1, term2, term3, nupper, nval, consumption)
                end if
            elseif (j.lt.jc .and. j.ge.1_ik) then !studying
                nupper = 1.0_rk -tbar
                call bisectionlabor( aval, term0, term1, term2, term3, nupper, nval, consumption)

            else
                nupper = 1.0_rk
                nval = 0.0_rk
                consumption = yval-aval
            end if


            if (consumption .gt. 0.0_rk) then

                call gridlookup(anum, agrid, aval, aindex, aweight)
                if (ev(aindex,bindex) .le. -9999999_rk) then
                    aweight = 0.0_rk
                end if
                evalf = bweight*(ev(aindex, bindex)*aweight + ev(aindex+1_ik,bindex)*(1.0_rk - aweight)) +(1.0_rk-bweight)*(ev(aindex, bindex+1)*aweight + ev(aindex+1_ik,bindex+1)*(1.0_rk - aweight))

                if (j.lt.jc .and. estat.eq.1.0_rk) then
                    fval = ((consumption)**(1.0_rk - sigma))/(1.0_rk - sigma)-psi_c*((nval)**(1.0_rk+eta)/(1+eta)) + beta*evalf
                else
                    fval = ((consumption)**(1.0_rk - sigma))/(1.0_rk - sigma)-psi*((nval)**(1.0_rk+eta)/(1+eta)) + beta*evalf
                end if
            else
                fval = -9999999_rk
            end if
            fc = -fval
        end if
        if ( s1 .gt. 100_ik) then
            write(*,*) 'd-c',d-c
        end if
    end do

    if (fval .gt. -9999999_rk) then
        objective = fval
        afval = aval
        cval = consumption
    else
        objective = fval
        afval = 0.0
        cval =0.0
    end if


    end subroutine goldensection

    subroutine goldensectionb( tbar, chib, phi, term0, term1, term2, term3, blow, bhigh, bnum, bgrid, &
        evhjc,  objective, bval, cval,nval)

    integer(ik)::   bindex, bnum
    real(rk)::  phi, bvaltemp, a, blow, b,  bhigh, rg, c, bgrid(bnum), evalf, evhjc(bnum), fval,  fc, d, s1, fd, z, &
        objective, bval, cval, bweight, term0, term1, term2, term3, consumption,nval, chib, nupper, tbar

    intent(in)::term0, term1, term2, term3,  blow, bhigh,bnum, bgrid, &
        evhjc,  phi, chib, tbar
    intent(out):: objective, bval, cval, nval



    nupper = 1.0_rk- tbar
    a = blow;  b = bhigh
    rg = (3.0_rk - dsqrt(5.0_rk))/2.0_rk

    !First c evaluation
    c = a + rg*(b-a)

    bval = c
    bvaltemp =dmax1(phi + bval/(dble(jc)-1.0_rk), 0.0_rk)
    call bisectionlabor( bvaltemp, term0, term1, term2, term3+chib*bval, nupper, nval, consumption)
    !nval = ncollege
    !consumption =term3+term0*nval-bvaltemp
    if (consumption .gt. 0.0_rk) then
        call gridlookup(bnum, bgrid, bval, bindex, bweight)
        if( evhjc(bindex) .gt. -9999999_rk) then
            evalf = bweight*evhjc(bindex) +(1.0_rk-bweight)*evhjc(bindex+1)
        else
            evalf = evhjc(bindex+1)
            if (bindex+1 .eq. bnum) then
                bval = 0.0_rk
                bvaltemp =dmax1(phi + bval/(dble(jc)-1.0_rk), 0.0_rk)
                call bisectionlabor( bvaltemp, term0, term1, term2, term3+chib*bval, nupper, nval, consumption)
                !consumption =term3+term0*nval-bvaltemp
            end if
        end if
        if ( consumption .gt. 0.0_rk) then
            fval = ((consumption)**(1.0_rk - sigma))/(1.0_rk - sigma)-psi_c*((nval)**(1.0_rk+eta)/(1+eta))+ beta*evalf
            if ( bval .lt. 0.0_rk) then
                fval = fval!-chib
            end if
        else
            fval = -9999999_rk
        end if

    else
        fval = -9999999_rk
    end if

    fc = -fval

    !Do it again for d

    d = a + (1-rg)*(b-a)
    bval = d
    bvaltemp =dmax1(phi + bval/(dble(jc)-1.0_rk), 0.0_rk)
    call bisectionlabor( bvaltemp, term0, term1, term2, term3+chib*bval, nupper, nval, consumption)
    !nval = ncollege
    !consumption =term3+term0*nval-bvaltemp

    if (consumption .gt. 0.0_rk) then
        call gridlookup(bnum, bgrid, bval, bindex, bweight)
        if( evhjc(bindex) .gt. -9999999_rk) then
            evalf = bweight*evhjc(bindex) +(1.0_rk-bweight)*evhjc(bindex+1)
        else
            evalf = evhjc(bindex+1)
            if (bindex+1 .eq. bnum) then
                bval = 0.0_rk
                bvaltemp =dmax1(phi + bval/(dble(jc)-1.0_rk), 0.0_rk)
                call bisectionlabor( bvaltemp, term0, term1, term2, term3+chib*bval, nupper, nval, consumption)
                !consumption =term3+term0*nval-bvaltemp
            end if
        end if
        if ( consumption .gt. 0.0_rk) then
            fval = ((consumption)**(1.0_rk - sigma))/(1.0_rk - sigma)-psi_c*((nval)**(1.0_rk+eta)/(1+eta))+ beta*evalf
            if ( bval .lt. 0.0_rk) then
                fval = fval!-chib
            end if
        else
            fval = -9999999_rk
        end if

    else
        fval = -9999999_rk
    end if


    fd = -fval

    !npdf
    s1 = 0;

    !As long as c and d differ, resample an intermediate value z
    do while ((d-c) > precision_gss)
        s1 = s1 + 1
        if (fc.ge.fd) then
            z = c + (1.0_rk-rg)*(b-c)
            !Case 1: [a c d b ] <--- [c d z b]
            a = c
            c = d
            fc = fd
            d = z

            bval = d
            bvaltemp =dmax1(phi + bval/(dble(jc)-1.0_rk), 0.0_rk)
            call bisectionlabor( bvaltemp, term0, term1, term2, term3+chib*bval, nupper, nval, consumption)
            !nval = ncollege
            !consumption =term3+term0*nval-bvaltemp

            if (consumption .gt. 0.0_rk) then
                call gridlookup(bnum, bgrid, bval, bindex, bweight)
                if( evhjc(bindex) .gt. -9999999_rk) then
                    evalf = bweight*evhjc(bindex) +(1.0_rk-bweight)*evhjc(bindex+1)
                else
                    evalf = evhjc(bindex+1)
                    if (bindex+1 .eq. bnum) then
                        bval = 0.0_rk
                        bvaltemp =dmax1(phi + bval/(dble(jc)-1.0_rk), 0.0_rk)
                        call bisectionlabor( bvaltemp, term0, term1, term2, term3+chib*bval, nupper, nval, consumption)
                        !consumption =term3+term0*nval-bvaltemp
                    end if
                end if
                if ( consumption .gt. 0.0_rk) then
                    fval = ((consumption)**(1.0_rk - sigma))/(1.0_rk - sigma)-psi_c*((nval)**(1.0_rk+eta)/(1+eta)) + beta*evalf
                    if ( bval .lt. 0.0_rk) then
                        fval = fval!-chib
                    end if
                else
                    fval = -9999999_rk
                end if

            else
                fval = -9999999_rk
            end if


            fd = -fval
        else
            z = a + rg*(d-a)
            !Case 2: [a c d b] <--- [a z c d]
            b = d
            d = c
            fd = fc
            c = z

            bval = c
            bvaltemp =dmax1(phi + bval/(dble(jc)-1.0_rk), 0.0_rk)
            call bisectionlabor( bvaltemp, term0, term1, term2, term3+chib*bval, nupper, nval, consumption)
            !nval = ncollege
            !consumption =term3+term0*nval-bvaltemp

            if (consumption .gt. 0.0_rk) then
                call gridlookup(bnum, bgrid, bval, bindex, bweight)
                if( evhjc(bindex) .gt. -9999999_rk) then
                    evalf = bweight*evhjc(bindex) +(1.0_rk-bweight)*evhjc(bindex+1)
                else
                    evalf = evhjc(bindex+1)
                    if (bindex+1 .eq. bnum) then
                        bval = 0.0_rk
                        bvaltemp =dmax1(phi + bval/(dble(jc)-1.0_rk), 0.0_rk)
                        call bisectionlabor( bvaltemp, term0, term1, term2, term3+chib*bval, nupper, nval, consumption)
                        !consumption =term3+term0*nval-bvaltemp
                    end if
                end if
                if ( consumption .gt. 0.0_rk) then
                    fval = ((consumption)**(1.0_rk - sigma))/(1.0_rk - sigma)-psi_c*((nval)**(1.0_rk+eta)/(1+eta)) + beta*evalf
                    if ( bval .lt. 0.0_rk) then
                        fval = fval!-chib
                    end if
                else
                    fval = -9999999_rk
                end if

            else
                fval = -9999999_rk
            end if

            fc = -fval
        end if
    end do

    if (fval .gt. -9999999_rk) then
        objective = fval
        cval = consumption
    else
        objective = fval
        cval =0.0
        bval = 0.0
    end if

    if (bval .gt. -1.0E-003_rk) then
        bval = 0.0_rk
    end if


    end subroutine goldensectionb

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Subroutine: Initial Distribution !!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine distini(t, chibprob, nabil, anum, adnum, bdnum, adgrid, probcost,  pop, ad0index, ad0weight, &
        edudecision,  bfchoice, bdgrid, mulini, muhcollege, muhcollege_cur, partrans_c, partrans_nc, partrans)

    integer(ik):: nabil, anum, adnum, bdnum, iab, ia, ic, ilocad,  ilocbd, ad0index(nabil, anum), ichib,j, t

    real(rk)::  adgrid(adnum), pop(jnum), edudecision(nabil,anum, ncost,chibnum), bfchoice(nabil,anum, ncost,chibnum),   wad,  &
        mulini(nabil,adnum,bdnum), bdgrid(bdnum), wbd, bval, ad0weight(nabil, anum), &
        tempc, tempnc, partrans_c, partrans_nc, partrans, probcost(nabil, anum, ncost), chibprob(chibnum), &
        muhcollege(jc-1, nabil, adnum, bdnum, chibnum), muhcollege_cur(jc-1, nabil, adnum, bdnum, chibnum)

    intent(in):: chibprob, anum,  adnum, bdnum, adgrid, pop, ad0index, ad0weight, &
        edudecision,  bfchoice, bdgrid,  probcost, nabil, t, muhcollege_cur
    intent(out):: mulini,partrans_c, partrans_nc, partrans, muhcollege

    mulini = 0.0_rk; partrans_c=0.0_rk; tempc = 0.0_rk; partrans_nc=0.0_rk; tempnc = 0.0_rk; muhcollege = 0.0_rk

    do ic = 1,ncost
        do ia = 1, anum
            do iab = 1, nabil
                ilocad = ad0index(iab, ia); wad = ad0weight(iab, ia)
                do ichib = 1, chibnum
                    if (edudecision(iab,ia,ic,ichib) .eq. 1.0_rk) then
                        bval = bfchoice(iab,ia,ic,ichib)
                        call gridlookup(bdnum, bdgrid, bval, ilocbd, wbd)
                        muhcollege(1, iab,ilocad,ilocbd, ichib) = muhcollege(1, iab,ilocad,ilocbd, ichib) + pop(1)*wad*wbd*probcost(iab,ia, ic)*chibprob(ichib)
                        muhcollege(1, iab,ilocad+1,ilocbd, ichib) = muhcollege(1, iab,ilocad+1,ilocbd, ichib) + pop(1)*(1.0_rk - wad)*wbd*probcost(iab,ia, ic)*chibprob(ichib)
                        muhcollege(1, iab,ilocad,ilocbd+1, ichib) = muhcollege(1, iab,ilocad,ilocbd+1, ichib) + pop(1)*wad*(1.0_rk-wbd)*probcost(iab,ia, ic)*chibprob(ichib)
                        muhcollege(1, iab,ilocad+1,ilocbd+1, ichib) = muhcollege(1, iab,ilocad+1,ilocbd+1, ichib) + pop(1)*(1.0_rk - wad)*(1.0_rk-wbd)*probcost(iab,ia, ic)*chibprob(ichib)

                    else !Those who decides to enter workforce
                        mulini(iab,ilocad,bdnum) = mulini(iab,ilocad,bdnum) + pop(1)*wad*probcost(iab,ia, ic)*chibprob(ichib)
                        mulini(iab,ilocad + 1,bdnum) = mulini(iab,ilocad + 1,bdnum) + pop(1)*(1.0_rk - wad)*probcost(iab,ia, ic)*chibprob(ichib)
                    end if
                end do
            end do
        end do
    end do

    !do j = 2, jc-1
    !    if (t .eq. 0_ik) then
    !        muhcollege(j, :,:,:,:) = muhcollege(j-1,:,:,:,:)
    !    else
    !        muhcollege(j, :,:,:,:) = muhcollege_cur(j-1,:,:,:,:)
    !    end if
    !end do


    do ia = 1, adnum
        tempc = tempc + adgrid(ia)*sum(muhcollege(1,1:nabil,ia, 1:bdnum,1:chibnum))
        tempnc = tempnc +adgrid(ia)*sum(mulini(1:nabil,ia,1:bdnum))
    end do

    !parental transfers for college/ noncollege
    partrans_c = tempc/sum(muhcollege(1,:,:,:,:))
    partrans_nc = tempnc/sum(mulini)
    partrans = (tempc + tempnc)/(sum(muhcollege(1,:,:,:,:)) + sum(mulini))

    end subroutine distini

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!                                       !!!!
    !!!!             distribution              !!!!
    !!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine entiredist(pabil, phid_abil, chibprob, nabil, eagridj, anum, bnum, adnum, bdnum, wp, probcost, pop, ad0index, ad0weight,  adgrid, edudecision, &
        piel0, piec0, bfchoice, bdgrid, laborhcollege, laborh, t,  r_b, lambda, pchoicec, pchoicenc, ainiindex, ainiweight, ghcollege, gh, gl,   &
        aindex, aweight, bindex, bweight, piec, piel,  laborl, muhcollege_cur, muh_cur, mul_cur,  mud_cur,  &
        kagg, lagg, yagg,  debtout, navg, educated, noneducated, stddebtholder, parentaltransferc, parentaltransfernc, parentaltransfer, debtbyage, ddebt, &
        fracborrower, avgdebt, avgdebt_c, avgdebt_d, totaldebtbyage, ddebtpercent, totalborrower, ddebtbyage, ddebtmubyage, ddebtpercentmu, colrate, drateh, drated)

    integer(ik)::   anum, bnum, adnum, bdnum,   t, j, ip,  ia, ib, iab, nabil
    integer(ik):: ainiindex(nabil, adnum), aindex(adnum), bindex(bdnum), pchoicec(jpaynum,nT,nabil,e_num,anum,bnum), &
        pchoicenc(jpaynum,nT,nabil,e_num,anum,bnum), ad0index(nabil,anum)

    real(rk)::  kagg, laggh, laggl, laggd, lagg, yagg, r_b, ntotald, debtout, &
        ntotalh, ntotall, avgdebt, avgdebt_c, avgdebt_d,  educated, noneducated, stddebtholder,laggh_jc, &
        navg, parentaltransferc, parentaltransfernc, parentaltransfer, debtbyage(5), fracborrower(5), lambdaval, wp, &
        muj(5), totaldebtbyage(5), ddebtpercent, totalborrower, ddebtpercentmu, debtout_age(jnum), &
        popcollege, poph, popl, popd, drateh, drated

    real(rk)::  pop(jnum), adgrid(adnum), bdgrid(bdnum), piel0(e_num), ainiweight(nabil,adnum), aweight(adnum), bweight(bdnum),  &
        lambda(nT), edudecision(nabil, anum, ncost,chibnum), eagridj(jnum, nabil,e_num, 2), ad0weight(nabil,anum), &
        bfchoice(nabil,anum, ncost,chibnum),  piec(e_num,e_num), piel(e_num,e_num),  mulini(nabil,adnum,bdnum),     &
        muhcollege_cur(jc-1, nabil, adnum, bdnum, chibnum), phid_abil(nabil), pabil(nabil), &
        gl(jnum,nT,nabil,e_num,anum,bnum), laborh(jrnum,nT,nabil,e_num,anum,bnum), gh(jnum,nT,nabil,e_num,anum,bnum), laborl(jrnum,nT,nabil,e_num,anum,bnum),&
        muh_cur(jnum,nT,nabil,e_num,adnum,bdnum), mul_cur(jnum,nT,nabil,e_num,adnum,bdnum), mud_cur(jnum,nT,nabil,e_num,adnum,bdnum),  probcost(nabil, anum, ncost), &
        ddebth(jnum), ddebtd(jnum), ddebt(jnum), ddebtbyage(5), ddebtmu(jnum), ddebtmuh(jnum), ddebtmud(jnum),&
        ddebtmubyage(5), ntotaleduc, piec0(e_num), chibprob(chibnum), colrate, ghcollege(jc-1_ik,nabil,anum,bnum,chibnum), laborhcollege(jc-1_ik,nabil,anum,bnum,chibnum)


    real(rk), allocatable:: mud(:,:,:,:,:,:), muh(:,:,:,:,:,:), mul(:,:,:,:,:,:), muhcollege(:,:,:,:,:)

    intent(in):: nabil, anum, bnum, adnum, bdnum,  pop, ad0index, ad0weight,  adgrid, edudecision, piel0, bfchoice, bdgrid, &
        laborh, t,  r_b, lambda, pchoicec, pchoicenc, ainiindex, ainiweight, gh, gl,  aindex, aweight, bindex, bweight, &
        piec, piel,  laborl, probcost, eagridj, piec0, chibprob, ghcollege, laborhcollege, phid_abil, pabil

    intent(out):: totaldebtbyage, ddebtpercent, kagg, lagg, yagg,  debtout, &
        navg, totalborrower, parentaltransferc, parentaltransfernc, parentaltransfer, &
        debtbyage, ddebt, fracborrower, avgdebt, avgdebt_c, avgdebt_d, ddebtbyage, ddebtmubyage, ddebtpercentmu, &
        stddebtholder, colrate, drateh, drated

    intent(inout):: muh_cur, mul_cur, mud_cur, muhcollege_cur


    if (t.eq.0_ik) then
        write(*,*) 'entered entiredist for steady state'
    end if
    allocate( muhcollege(jc-1, nabil, adnum, bdnum, chibnum), mud(jnum,nT,nabil,e_num,adnum,bdnum), muh(jnum,nT,nabil,e_num,adnum,bdnum), &
        mul(jnum,nT,nabil,e_num,adnum,bdnum))

    mud = 0.0_rk; muh = 0.0_rk; mul =0.0_rk; muhcollege = 0.0_rk

    call distini(t, chibprob, nabil, anum,  adnum, bdnum, adgrid, probcost,  pop, ad0index, ad0weight,  edudecision,  &
        bfchoice, bdgrid, mulini, muhcollege, muhcollege_cur, parentaltransferc, parentaltransfernc, parentaltransfer)

    call disteduc(nabil, anum, bnum, adnum, bdnum, phid_abil, piel0,  t, muhcollege_cur, ainiindex,ainiweight,   &
        laborhcollege,ghcollege, muhcollege, adgrid, bindex, bweight,  eagridj(:,:,:,1), piec0, muh, laggh_jc, mud,  ntotaleduc)

    !skilled workers
    call disth(jc-1_ik,nabil, anum, bnum, adnum, bdnum, laborh, t, muh_cur, aindex,  aweight,  r_b, lambda, pchoicec, gh, &
        adgrid, bdgrid, bindex, bweight, eagridj(:,:,:,1),  piec, muh, laggh, ddebth, ntotalh, ddebtmuh)

    !drop-outs
    call disth(jd-1_ik, nabil, anum, bnum, adnum, bdnum, laborl,  t, mud_cur, aindex,  aweight,  r_b, lambda, pchoicenc, gl, &
        adgrid, bdgrid, bindex, bweight,  eagridj(:,:,:,2), piel, mud, laggd, ddebtd, ntotald, ddebtmud)

    !unskilled workers
    call distl(nabil, anum, bnum, adnum, bdnum, bindex, bweight, laborl, t, mul_cur, aindex, ainiindex, aweight, ainiweight, &
        gl, mulini,  adgrid, eagridj(:,:,:,2), piel,piel0, mul, laggl, ntotall)

    if ( t .eq. 0_ik) then
        mul_cur = mul
        muh_cur = muh
        mud_cur = mud
        muhcollege_cur = muhcollege
    end if

    ! debtout_c: College graduates' total outstanding student loan debt
    ! debtout_nc: Dropout's total outstanding student loan debt
    ! debtout: Total outstanding student loan debt
    ! avgdebt: Average student loan debt among graduating seniors with positive student debt

    avgdebt = 0.0_rk; kagg = 0.0_rk; ddebt = 0.0_rk; ddebtmu = 0.0_rk; debtout_age = 0.0_rk;
    debtout = 0.0_rk; avgdebt_c = 0.0_rk; avgdebt_d=0.0_rk


    do ia = 1, adnum
        do ib = 1, bdnum
            popcollege = sum(muhcollege_cur(:,:,ia,ib,:))
            poph = sum(muh_cur(:,:,:,:,ia,ib))
            popl = sum(mul_cur(:,:,:,:,ia,ib))
            popd = sum(mud_cur(:,:,:,:,ia,ib))

            kagg = kagg + adgrid(ia)*(popcollege+poph+popl+popd)
            debtout = debtout + bdgrid(ib)*(popcollege+poph+popd)

            avgdebt = avgdebt + bdgrid(ib)*(sum(muhcollege_cur(jd-1,:,ia,ib,:)))
            avgdebt_c = avgdebt_c + bdgrid(ib)*(sum(muh_cur(jc,1,:,:,ia,ib)))
            avgdebt_d = avgdebt_d + bdgrid(ib)*(sum(mud_cur(jd,1,:,:,ia,ib)))
        end do
    end do

    avgdebt = (avgdebt_c+avgdebt_d)/(sum(muh_cur(jc,1,:,:,:,:))+(sum(mud_cur(jd,1,:,:,:,:))))!average debt among graduating seniors.
    avgdebt_c = avgdebt_c/(sum(muh_cur(jc,1,:,:,:,:)))!average debt among non-college graduates.
    avgdebt_d = avgdebt_d/(sum(mud_cur(jd,1,:,:,:,:))) !average debt among college dropouts.
    do j = 1, jnum
        do ib = 1, bdnum
            if ( j .lt. jc) then
                debtout_age(j) = debtout_age(j) + (bdgrid(ib)*(sum(muhcollege_cur(j,:,:,ib,:))+sum(mud_cur(j,:,:,:,:,ib))))
            else
                debtout_age(j) = debtout_age(j) + (bdgrid(ib)*sum(muh_cur(j,:,:,:,:,ib)) + bdgrid(ib)*sum(mud_cur(j,:,:,:,:,ib)))
            end if
        end do
    end do


    noneducated = (sum(mul_cur(jc:jnum,:,:,:,:,:)) + sum(mud_cur(jc:jnum,:,:,:,:,:)))/(sum(mul_cur(jc:jnum,:,:,:,:,:)) + sum(mud_cur(jc:jnum,:,:,:,:,:)) + sum(muh_cur(jc:jnum,:,:,:,:,:)))
    educated = 1.0_rk- noneducated

    colrate =0.0_rk
    do iab=1, nabil
        colrate = colrate + (1-phid_abil(iab))*(sum(muhcollege_cur(1:jd-1_ik,iab,:,:,:)))/(sum(mul_cur(1:jd-1_ik,:,iab,:,:,:)) + sum(muhcollege_cur(1:jd-1_ik,iab,:,:,:)))*pabil(iab)
    end do
    
    
    stddebtholder = sum(muhcollege_cur(1:jc-1_ik,1:nabil,1:adnum,1:bdnum-1,1:chibnum))/sum(muhcollege_cur(1:jc-1_ik,1:nabil,1:adnum,1:bdnum,1:chibnum))!!calculate the fraction of BA holders with student debt.
    navg =ntotaleduc+ ntotall + ntotalh + ntotald; navg = navg/(sum(muhcollege_cur)+sum(mul_cur(1:jrnum,:,:,:,:,:)) + sum(mud_cur(1:jrnum,:,:,:,:,:)) + sum(muh_cur(1:jrnum,:,:,:,:,:)))

    !!!!!Calculate the prices!!!!!
    laggl = laggl + laggd
    laggh = laggh + laggh_jc

    lambdaval = wp*((laggh/laggl)**(1/theta))/(1 + wp*((laggh/laggl)**(1/theta)))
    lagg = (lambdaval*(laggh**((theta-1)/theta)) + (1-lambdaval)*(laggl**((theta-1)/theta)))**(theta/(theta-1))
    yagg = (kagg**alpha)*(lagg**(1-alpha))

    !debt by age group
    !age below 30 are those with j < 12.
    !age between 30 and 39 are those with 12 <= j < 22.
    !age between 40 and 49 are those with 22 <= j < 32
    !age between 50 and 59 are those with 32 <= j < 42
    !age above 60 are those with j <= 42
    ! % balance 90+ days delinquent student debt
    ddebt(1:jnum) = ddebth(1:jnum) + ddebtd(1:jnum)
    ddebtmu(1:jnum) = ddebtmuh(1:jnum) + ddebtmud(1:jnum) !population who are delinquent
    ddebtpercent = sum(ddebt)/debtout

    debtbyage = 0.0_rk; fracborrower = 0.0_rk; ddebtbyage = 0.0_rk; ddebtmubyage = 0.0_rk
    do j = 1,jnum
        if (j.lt.12_ik) then
            ddebtbyage(1) = ddebtbyage(1) + ddebt(j);
            ddebtmubyage(1) = ddebtmubyage(1) + ddebtmu(j);
        elseif (j.ge.12 .and. j.lt.22) then
            ddebtbyage(2) = ddebtbyage(2) + ddebt(j);
            ddebtmubyage(2) = ddebtmubyage(2) + ddebtmu(j);
        elseif (j.ge.22 .and. j.lt.32) then
            ddebtbyage(3) = ddebtbyage(3) + ddebt(j);
            ddebtmubyage(3) = ddebtmubyage(3) + ddebtmu(j);
        elseif (j.ge.32 .and. j.lt.42) then
            ddebtbyage(4) = ddebtbyage(4) + ddebt(j);
            ddebtmubyage(4) = ddebtmubyage(4) + ddebtmu(j);
        else !(j.ge.42 .and. j.lt. jnum)
            ddebtbyage(5) = ddebtbyage(5) + ddebt(j);
            ddebtmubyage(5) = ddebtmubyage(5) + ddebtmu(j);
        end if

        do ib= 1, bdnum - 1_ik !bdnum will result in zero.

            if (j.lt.12_ik) then
                if ( j .lt. jc) then
                    debtbyage(1) = debtbyage(1) + bdgrid(ib)*sum(muhcollege_cur(j, 1:nabil, 1:adnum, ib, 1:chibnum))
                    fracborrower(1) = fracborrower(1) + sum(muhcollege_cur(j, 1:nabil, 1:adnum, ib, 1:chibnum))
                else
                    debtbyage(1) = debtbyage(1) + bdgrid(ib)*sum(muh_cur(j, 1:nT, 1:nabil,1:e_num, 1:adnum, ib)) + bdgrid(ib)*sum(mud_cur(j, 1:nT, 1:nabil,1:e_num, 1:adnum, ib))
                    fracborrower(1) = fracborrower(1) + sum(muh_cur(j, 1:nT, 1:nabil,1:e_num, 1:adnum, ib)) + sum(mud_cur(j, 1:nT, 1:nabil,1:e_num, 1:adnum, ib))
                end if

            elseif (j.ge.12 .and. j.lt.22) then
                debtbyage(2) = debtbyage(2) + bdgrid(ib)*sum(muh_cur(j, 1:nT, 1:nabil,1:e_num, 1:adnum, ib)) + bdgrid(ib)*sum(mud_cur(j, 1:nT, 1:nabil,1:e_num, 1:adnum, ib))
                fracborrower(2) = fracborrower(2) + sum(muh_cur(j, 1:nT,1:nabil, 1:e_num, 1:adnum, ib)) +sum(mud_cur(j, 1:nT,1:nabil, 1:e_num, 1:adnum, ib))

            elseif (j.ge.22 .and. j.lt.32) then
                debtbyage(3) = debtbyage(3) + bdgrid(ib)*sum(muh_cur(j, 1:nT,1:nabil, 1:e_num, 1:adnum, ib)) + bdgrid(ib)*sum(mud_cur(j, 1:nT, 1:nabil,1:e_num, 1:adnum, ib))
                fracborrower(3) = fracborrower(3) + sum(muh_cur(j, 1:nT,1:nabil, 1:e_num, 1:adnum, ib)) + sum(mud_cur(j, 1:nT,1:nabil, 1:e_num, 1:adnum, ib))

            elseif (j.ge.32 .and. j.lt.42) then
                debtbyage(4) = debtbyage(4) + bdgrid(ib)*sum(muh_cur(j, 1:nT, 1:nabil,1:e_num, 1:adnum, ib)) + bdgrid(ib)*sum(mud_cur(j, 1:nT,1:nabil, 1:e_num, 1:adnum, ib))
                fracborrower(4) = fracborrower(4) + sum(muh_cur(j, 1:nT, 1:nabil,1:e_num, 1:adnum, ib))+ sum(mud_cur(j, 1:nT,1:nabil, 1:e_num, 1:adnum, ib))

            else !(j.ge.42 .and. j.lt. jnum)
                debtbyage(5) = debtbyage(5) + bdgrid(ib)*sum(muh_cur(j, 1:nT,1:nabil, 1:e_num, 1:adnum, ib)) + bdgrid(ib)*sum(mud_cur(j, 1:nT, 1:nabil,1:e_num, 1:adnum, ib))
                fracborrower(5) = fracborrower(5) + sum(muh_cur(j, 1:nT, 1:nabil,1:e_num, 1:adnum, ib))+ sum(mud_cur(j, 1:nT,1:nabil, 1:e_num, 1:adnum, ib))
            end if
        end do
    end do
    
    ! delinquency moments
    ddebtpercentmu = (sum(ddebtmuh(jc:jrnum))+sum(ddebtmud(jd:jrnum)))/(sum(muh_cur(jc:jrnum,:,:,:,:,1:bnum-1))+sum(mud_cur(jd:jrnum,:,:,:,:,1:bnum-1)))
    drateh = sum(ddebtmuh(jc:jrnum))/(sum(muh_cur(jc:jrnum,:,:,:,:,1:bnum-1)))
    drated = sum(ddebtmud(jd:jrnum))/(sum(mud_cur(jd:jrnum,:,:,:,:,1:bnum-1)))
    
    ddebtmubyage(1:5) = ddebtmubyage(1:5)/fracborrower(1:5);
    ddebtbyage(1:5) = ddebtbyage(1:5)/debtbyage(1:5);
    totalborrower = sum(fracborrower)
    fracborrower = fracborrower/sum(fracborrower)

    totaldebtbyage = 0.0_rk;muj =0.0_rk
    muj(1:4)= 11.0*(1.0_rk/jnum)
    muj(5) = (jnum-42)*(1.0_rk/jnum)
    do ip = 1,5
        totaldebtbyage(ip) = (debtbyage(ip))/muj(ip)
    end do
    totaldebtbyage = -totaldebtbyage/yagg

    mul_cur = mul
    muh_cur = muh
    mud_cur = mud
    muhcollege_cur = muhcollege
    end subroutine entiredist




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Subroutine: Dist for College-Grad !!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine disteduc( nabil, anum, bnum, adnum, bdnum,  phid_abil,  piel0, t, muhcollege_cur, &
        ainiindex,ainiweight,laborhcollege,ghcollege,muhcollege,adgrid,bindex, bweight, eagridj, piec0,muh,laggh_jc, mud, nvalh)

    integer(ik)::  nabil, anum, bnum, adnum, bdnum,  j, iloca, ia, ib, ilocadh, ip, ilocb, t, iab,ichib
    integer(ik):: ainiindex(nabil, adnum), bindex(bdnum)

    real(rk):: laggh_jc, muvalh, wa, wb, afvalh, wadh, nvalh, nvaltemp,  phid
    real(rk):: adgrid(adnum), piec0(e_num), muh(jnum,nT,nabil, e_num,adnum,bdnum), &
        muhcollege(jc-1, nabil, adnum, bdnum, chibnum),muhcollege_cur(jc-1, nabil, adnum, bdnum, chibnum), mud(jnum,nT,nabil, e_num,adnum,bdnum), &
        ainiweight(nabil, adnum), bweight(bdnum),  ghcollege(jc-1,nabil,anum,bnum,chibnum),laborhcollege(jc-1,nabil,anum,bnum,chibnum),&
        piel0(e_num), eagridj(jnum, nabil, e_num), phid_abil(nabil)

    intent(in)::  nabil, anum, bnum, adnum, bdnum,  piel0,  adgrid,  piec0, phid_abil, &
        ainiindex,ainiweight, bindex, bweight, t, muhcollege_cur,   eagridj,ghcollege,laborhcollege
    intent(out)::   laggh_jc, nvalh
    intent(inout):: muh, mud, muhcollege

    laggh_jc = 0.0_rk; muh =0.0_rk; nvalh = 0.0_rk; mud = 0.0_rk;
    do j = 1, jc-1
        do iab = 1, nabil
            phid =phid_abil(iab)
            do ia = 1, adnum
                do ib = 1, bdnum
                    do ichib = 1, chibnum

                        if (t .eq. 0_ik) then
                            muvalh = muhcollege(j,iab,ia,ib, ichib);
                        else
                            muvalh = muhcollege_cur(j,iab,ia,ib, ichib);
                        end if

                        if (muvalh.gt.0.0_rk) then

                            iloca =ainiindex(iab, ia); wa = ainiweight(iab, ia);
                            ilocb = bindex(ib); wb = bweight(ib)

                            afvalh = wb*(wa*ghcollege(j,iab, iloca,ilocb, ichib) + (1.0_rk-wa)*ghcollege(j,iab, iloca+1,ilocb, ichib))+ &
                                (1.0-wb)*(wa*ghcollege(j,iab, iloca,ilocb+1, ichib) + (1.0_rk-wa)*ghcollege(j,iab, iloca+1,ilocb+1, ichib))
                            nvaltemp = wb*(wa*laborhcollege(j,iab, iloca,ilocb, ichib) + (1.0_rk-wa)*laborhcollege(j,iab, iloca+1,ilocb, ichib))+ &
                                (1.0-wb)*(wa*laborhcollege(j,iab, iloca,ilocb+1, ichib) + (1.0_rk-wa)*laborhcollege(j,iab, iloca+1,ilocb+1, ichib))

                            laggh_jc = laggh_jc + eagridj(j,iab,1)*muvalh*nvaltemp
                            nvalh = nvalh + nvaltemp*muvalh
                            if ( jd .eq. jc) then
                                if ( j.eq. jc-1_ik) then
                                    call gridlookup(adnum, adgrid, afvalh, ilocadh, wadh)
                                    ip = 1_ik

                                    if (ilocadh .lt. adnum) then
                                        muh(j+1,ip,iab,:,ilocadh,ib) = muh(j+1,ip,iab,:,ilocadh,ib) + wadh*piec0(:)*muvalh*(1.0_rk-phid)
                                        muh(j+1,ip,iab,:,ilocadh+1,ib) = muh(j+1,ip,iab,:,ilocadh + 1,ib) + (1.0_rk-wadh)*piec0(:)*muvalh*(1.0_rk-phid)

                                        mud(j+1,ip,iab,:,ilocadh,ib) = mud(j+1,ip,iab,:,ilocadh,ib) + wadh*muvalh*phid*piel0(:)
                                        mud(j+1,ip,iab,:,ilocadh+1,ib) = mud(j+1,ip,iab,:,ilocadh + 1,ib) + (1.0_rk-wadh)*muvalh*phid*piel0(:)
                                    else
                                        muh(j+1,ip,iab,:,ilocadh,ib) = muh(j+1,ip,iab,:,ilocadh,ib) + wadh*muvalh*(1.0_rk-phid)*piec0(:)
                                        mud(j+1,ip,iab,:,ilocadh,ib) = mud(j+1,ip,iab,:,ilocadh,ib) + wadh*muvalh*phid*piel0(:)
                                    end if
                                else
                                    muhcollege(j+1,iab,ia,ib, ichib) = muvalh
                                end if
                            else

                                if ( j.eq. jc-1_ik) then
                                    call gridlookup(adnum, adgrid, afvalh, ilocadh, wadh)
                                    ip = 1_ik

                                    if (ilocadh .lt. adnum) then
                                        muh(j+1,ip,iab,:,ilocadh,ib) = muh(j+1,ip,iab,:,ilocadh,ib) + wadh*piec0(:)*muvalh
                                        muh(j+1,ip,iab,:,ilocadh+1,ib) = muh(j+1,ip,iab,:,ilocadh + 1,ib) + (1.0_rk-wadh)*piec0(:)*muvalh
                                    else
                                        muh(j+1,ip,iab,:,ilocadh,ib) = muh(j+1,ip,iab,:,ilocadh,ib) + wadh*muvalh*piec0(:)
                                    end if
                                elseif (j .eq. jd-1_ik) then
                                    ip = 1_ik
                                    muhcollege(j+1,iab,ia,ib, ichib) = muhcollege(j+1,iab,ia,ib, ichib) + muvalh*(1.0_rk-phid)
                                    mud(j+1,ip,iab,:,ia,ib) = mud(j+1,ip,iab,:,ia,ib) + muvalh*phid*piel0(:)
                                else
                                    muhcollege(j+1,iab,ia,ib, ichib) = muvalh
                                end if
                            end if

                        end if
                    end do
                end do
            end do

        end do
    end do

    
    end subroutine disteduc

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Subroutine: Dist for College-Grad !!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine disth(jin, nabil, anum, bnum, adnum, bdnum, laborh,  t, muh_cur,aindex,  aweight,  r_b, lambda, pchoicec,gh,adgrid,&
        bdgrid, bindex, bweight, eagridj,piec,muh,laggh, ddebt,nvalh, ddebtmu)

    integer(ik):: nabil, anum, bnum, adnum, bdnum,    j,jin, iloca, ie, ief, ia, ib, ibf, ilocadh, pchoicec(jpaynum,nT,nabil, e_num,anum,bnum),ip,ipf, &
        aindex(adnum),  bindex(bdnum), ilocb, t, iab

    real(rk)::  gh(jnum,nT,nabil, e_num,anum,bnum), adgrid(adnum), bdgrid(bdnum), bval, bfval,r_b,  &
        piec(e_num,e_num), laggh, muvalh, wa, wb, afvalh, wadh, muh(jnum,nT,nabil, e_num,adnum,bdnum), &
        piectemp(e_num),pchoicetemp,lambda(nT), aweight(adnum),  bweight(bdnum),muh_cur(jnum,nT,nabil, e_num,adnum,bdnum),  &
        laborh(jrnum,nT,nabil, e_num,anum,bnum), nvalh,nvaltemp, ddebt(jnum), ddebtmu(jnum),eagridj(jnum, nabil, e_num)

    intent(in):: jin, nabil, anum, bnum, adnum, bdnum,  laborh, pchoicec,gh,adgrid,bdgrid,piec,lambda, r_b, &
        aindex, aweight,  bindex, bweight, t, muh_cur,eagridj

    intent(out):: laggh,ddebt, nvalh, ddebtmu
    intent(inout):: muh

    laggh = 0; nvalh = 0.0_rk ;  ddebt = 0.0_rk ;  ddebtmu = 0.0_rk;



    do j = jin, jnum-1
        do iab = 1, nabil
            do ie = 1, e_num
                do ia = 1, adnum
                    do ib = 1, bdnum

                        if (j.gt.jrnum) then! for retirees, no shock
                            piectemp = 0.0_rk;
                            piectemp(ie) = 1.0_rk;
                        else
                            piectemp = 0.0_rk;
                            piectemp(1:e_num) = piec(ie,1:e_num);
                        end if

                        do ip = 1,nT
                            if (t .eq. 0_ik) then
                                muvalh = muh(j,ip,iab, ie,ia,ib)
                            else
                                muvalh = muh_cur(j,ip,iab, ie,ia,ib)
                            end if

                            if (muvalh.gt.0.0_rk) then
                                iloca  = aindex(ia); wa = aweight(ia); ilocb = bindex(ib); wb = bweight(ib)
                                afvalh = wb*(wa*gh(j,ip,iab, ie,iloca,ilocb) + (1.0_rk-wa)*gh(j,ip,iab, ie,iloca+1,ilocb))+(1.0-wb)*(wa*gh(j,ip,iab, ie,iloca,ilocb+1) + (1.0_rk-wa)*gh(j,ip,iab, ie,iloca+1,ilocb+1))
                                if( j .le. jrnum) then
                                    nvaltemp = wb*(wa*laborh(j,ip,iab, ie,iloca,ilocb) + (1.0_rk-wa)*laborh(j,ip,iab, ie,iloca+1,ilocb))+(1.0-wb)*(wa*laborh(j,ip,iab, ie,iloca,ilocb+1) + (1.0_rk-wa)*laborh(j,ip,iab, ie,iloca+1,ilocb+1))
                                    laggh = laggh + eagridj(j,iab, ie)*muvalh*nvaltemp
                                    nvalh = nvalh + nvaltemp*muvalh
                                end if

                                if (j.le.jpaynum) then !between 5 and jpaynum -1. They can choose to pay or not.

                                    !For college graduates
                                    pchoicetemp = wb*(wa*pchoicec(j,ip,iab, ie,iloca,ilocb) + (1.0_rk - wa)*pchoicec(j,ip,iab, ie,iloca+1,ilocb))+ (1.0-wb)*(wa*pchoicec(j,ip,iab, ie,iloca,ilocb+1) + (1.0_rk - wa)*pchoicec(j,ip,iab, ie,iloca+1,ilocb+1))
                                    if (pchoicetemp .gt. 0.5_rk) then
                                        pchoicetemp = 1.0_rk
                                        bval = bdgrid(ib)
                                        bfval = (1.0_rk + r_b)*(1.0_rk - lambda(ip))*bval
                                        call gridlookup(bdnum, bdgrid, bfval, ibf, wb)
                                        ipf = min(ip+1,nT) !to avoid the index problem
                                    else
                                        pchoicetemp = 0.0_rk
                                        bval = bdgrid(ib)
                                        if (ib.ne.bdnum) then
                                            ddebt(j) = ddebt(j) + muvalh*bval
                                            ddebtmu(j) = ddebtmu(j) + muvalh
                                        end if
                                        bfval = (1.0_rk + r_b)*bval !j ranges from 5 to 14.
                                        call gridlookup(bdnum, bdgrid, bfval, ibf, wb)
                                        ipf = ip
                                    end if

                                else   !From age jpaynum+1, they need to pay.

                                    !check if there are case when retirees can't pay.
                                    bval = bdgrid(ib)
                                    bfval = (1.0_rk + r_b)*(1.0_rk - lambda(ip))*bval
                                    call gridlookup(bdnum, bdgrid, bfval, ibf, wb)
                                    ipf = min(ip+1,nT)
                                end if

                                call gridlookup(adnum, adgrid, afvalh, ilocadh, wadh)


                                do ief = 1, e_num
                                    if (ilocadh.lt.adnum .and. ibf .lt. bdnum) then
                                        muh(j+1,ipf,iab,ief,ilocadh,ibf) = muh(j+1,ipf,iab,ief,ilocadh,ibf) + wb*wadh*piectemp(ief)*muvalh
                                        muh(j+1,ipf,iab,ief,ilocadh+1,ibf) = muh(j+1,ipf,iab,ief,ilocadh+1,ibf) + wb*(1.0_rk-wadh)*piectemp(ief)*muvalh
                                        muh(j+1,ipf,iab,ief,ilocadh,ibf+1) = muh(j+1,ipf,iab,ief,ilocadh,ibf+1) + (1.0-wb)*wadh*piectemp(ief)*muvalh
                                        muh(j+1,ipf,iab,ief,ilocadh+1,ibf+1) = muh(j+1,ipf,iab,ief,ilocadh+1,ibf+1) + (1.0-wb)*(1.0_rk-wadh)*piectemp(ief)*muvalh
                                    elseif (ilocadh .eq. adnum .and. ibf .lt. bdnum) then
                                        muh(j+1,ipf,iab,ief,ilocadh,ibf) = muh(j+1,ipf,iab,ief,ilocadh,ibf) + wb*wadh*piectemp(ief)*muvalh
                                        muh(j+1,ipf,iab,ief,ilocadh,ibf+1) = muh(j+1,ipf,iab,ief,ilocadh,ibf+1) +(1.0-wb)*wadh*piectemp(ief)*muvalh
                                    elseif (ilocadh .lt. adnum .and. ibf .eq. bdnum) then
                                        muh(j+1,ipf,iab,ief,ilocadh,ibf) = muh(j+1,ipf,iab,ief,ilocadh,ibf) + wadh*wb*piectemp(ief)*muvalh
                                        muh(j+1,ipf,iab,ief,ilocadh+1,ibf) = muh(j+1,ipf,iab,ief,ilocadh+1,ibf) +(1.0-wadh)*wb*piectemp(ief)*muvalh
                                    else
                                        muh(j+1,ipf,iab,ief,ilocadh,ibf) = muh(j+1,ipf,iab,ief,ilocadh,ibf) +wb*wadh*piectemp(ief)*muvalh
                                    end if
                                end do
                                
                                !update delinquent amount and borrowers for the last age
                                if (j.eq.jnum-1_ik) then
                                    if (t .eq. 0_ik) then
                                        muvalh = muh(j+1,ip,iab, ie,ia,ib)
                                    else
                                        muvalh = muh_cur(j+1,ip,iab, ie,ia,ib)
                                    end if
                                    if (ib.ne.bdnum) then
                                        ddebt(j+1) = ddebt(j+1) + muvalh*bval
                                        ddebtmu(j+1) = ddebtmu(j+1) + muvalh
                                    end if
                                end if
                                
                            end if

                        end do

                    end do

                end do
            end do
        end do

    end do


    end subroutine disth

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Subroutine: Dist for Non-College !!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine distl(nabil, anum, bnum, adnum, bdnum, bindex,bweight,laborl,t,  mul_cur, aindex, ainiindex, aweight, ainiweight, gl,mulini,adgrid,&
        eagridj,piel,piel0, mul,laggl, nval)

    integer(ik):: nabil, anum, bnum, adnum, bdnum,  j, ie, ia, ilocadl, ief, iloca, t,ilocb,ip,ib,iab, &
        aindex(adnum), ainiindex(nabil, adnum),bindex(bdnum)
    real(rk)::  gl(jnum, nT,nabil,  e_num, anum, bnum), mulini(nabil,adnum,bdnum),  adgrid(adnum), piel(e_num,e_num),laborl(jrnum,nT,nabil, e_num,anum,bnum),nval,nvaltemp, &
        laggl, wa, afvall, wadl, muvall, mul(jnum,nT,nabil, e_num,adnum,bdnum), mul_cur(jnum,nT,nabil, e_num,adnum,bdnum), pieltemp(e_num), aweight(adnum), &
        ainiweight(nabil, adnum), wb,bweight(bdnum),eagridj(jnum, nabil,e_num), piel0(e_num)

    intent(in):: nabil, anum, bnum, adnum, bdnum, bindex,bweight,laborl,t, mul_cur,gl,mulini,adgrid,piel, &
        aindex, ainiindex, aweight, ainiweight,eagridj, piel0
    intent(out):: mul,laggl, nval

    laggl = 0.0_rk; mul = 0.0_rk; nval = 0.0_rk; nvaltemp = 0.0_rk;

    do ie = 1, e_num
        mul(1, 1, 1:nabil, ie, 1:adnum, 1:bdnum) = mulini(1:nabil, 1:adnum, 1:bdnum)*piel0(ie)
    end do


    do j = 1, jnum-1
        do iab = 1, nabil
            do ie = 1, e_num
                do ia = 1, adnum
                    do ip = 1, nT
                        do ib = 1, bdnum
                            if (t.eq. 0_ik) then
                                muvall = mul(j,ip,iab, ie,ia,ib)
                            else
                                muvall = mul_cur(j,ip,iab,ie,ia,ib)
                            end if

                            if (muvall .gt. 0.0_rk) then

                                if (j .eq. 1_ik) then
                                    iloca = ainiindex(iab, ia); wa = ainiweight(iab, ia)
                                else
                                    iloca = aindex(ia); wa = aweight(ia)
                                end if
                                ilocb = bindex(ib); wb = bweight(ib)

                                !Calculate average time worked and aggregate labor.
                                if (j.le.jrnum) then
                                    nvaltemp =  wb*(wa*laborl(j,ip,iab,ie,iloca,ilocb) + (1.0_rk-wa)*laborl(j,ip,iab,ie,iloca+1,ilocb))+(1.0-wb)*(wa*laborl(j,ip,iab,ie,iloca,ilocb+1) + (1.0_rk-wa)*laborl(j,ip,iab,ie,iloca+1,ilocb+1))
                                    laggl = laggl + eagridj(j,iab,ie)*muvall*nvaltemp
                                    nval = nval + nvaltemp*muvall
                                end if

                                afvall = wb*(wa*gl(j,ip,iab,ie,iloca,ilocb) + (1.0_rk-wa)*gl(j,ip,iab,ie,iloca+1,ilocb))+(1.0-wb)*(wa*gl(j,ip,iab,ie,iloca,ilocb+1) + (1.0_rk-wa)*gl(j,ip,iab,ie,iloca+1,ilocb+1))
                                call gridlookup(adnum, adgrid, afvall, ilocadl, wadl)


                                if (j.gt.jrnum) then !income shocks on workers only. Here, no schooling.
                                    pieltemp =0.0_rk
                                    pieltemp(ie) = 1.0_rk
                                else
                                    pieltemp =0.0_rk
                                    pieltemp(1:e_num) = piel(ie,1:e_num)
                                end if

                                do ief = 1, e_num
                                    if (ilocadl .lt. adnum) then !For these non-college, they don't have debt, so they stay at ip = 1 and ib = bdnum.
                                        mul(j+1,ip,iab,ief,ilocadl,ib) = mul(j+1,ip,iab,ief,ilocadl,ib) + wadl*pieltemp(ief)*muvall
                                        mul(j+1,ip,iab,ief,ilocadl + 1,ib) = mul(j+1,ip,iab,ief,ilocadl + 1,ib) + (1.0_rk-wadl)*pieltemp(ief)*muvall
                                    else
                                        mul(j+1,ip,iab, ief, ilocadl,ib) = mul(j+1,ip,iab,ief,ilocadl,ib) + wadl*pieltemp(ief)*muvall
                                    end if
                                end do

                            end if

                        end do
                    end do
                end do
            end do

        end do
    end do


    end subroutine distl


    subroutine indexcal(nabil, anum, bnum, adnum, bdnum, r_b,lambda,  bgrid, bdgrid, adgrid, a0grid, agrid, &
        aindex, ainiindex, aweight, ainiweight, bindex, bweight,payindex,payweight,nopayindex,nopayweight, ad0index, ad0weight)

    integer(ik):: ia, ie, iloca, iloca2, ilocad, ilocb, ib, ip,  bindex1, anum, bnum, adnum, bdnum, nabil
    integer(ik):: aindex(adnum), ainiindex(nabil, adnum), bindex(bdnum), nopayindex(bnum), payindex(bnum,nT),&
        ad0index(nabil, anum)

    real(rk):: adgrid(adnum), aweight(adnum), ainiweight(nabil, adnum),  a0grid(nabil, anum), agrid(anum), &
        bgrid(bnum), bdgrid(bdnum), bweight(bdnum), payweight(bnum,nT),nopayweight(bnum),lambda(nT), ad0weight(nabil,anum)
    real(rk):: aval, wa, wa2, wb, bval, bfval, wad,  r_b, bweight1

    intent(in)::    nabil, r_b, lambda,  adgrid, a0grid, agrid, bgrid, bdgrid,  anum, bnum, adnum, bdnum
    intent(out)::   aindex, ainiindex, aweight, ainiweight, bindex, bweight, payindex, payweight, nopayindex, nopayweight, ad0index, ad0weight

    aindex = 0_ik; ainiindex = 0_ik; aweight = 0.0_rk; ainiweight=0.0_rk;
    bindex = 0_ik; bweight = 0.0_rk; payindex = 0_ik; payweight = 0.0_rk;
    nopayindex = 0_ik; nopayweight = 0.0_rk

    !save index and weight in advance to save time
    do ia = 1,adnum
        aval = adgrid(ia)
        call gridlookup(anum, agrid, aval, iloca, wa)

        aindex(ia) = iloca
        aweight(ia) =wa

        do ie = 1, nabil
            aval = adgrid(ia)
            call gridlookup(anum, a0grid(ie, 1:anum), aval, iloca2, wa2)
            ainiindex(ie, ia) = iloca2
            ainiweight(ie, ia) = wa2
        end do
    end do

    do ie = 1,nabil
        do ia = 1,anum
            aval = a0grid(ie, ia)
            call gridlookup(adnum, adgrid, aval, ilocad, wad)
            ad0index(ie, ia) = ilocad
            ad0weight(ie, ia) =wad
        end do
    end do

    do ib =1,bdnum
        bval = bdgrid(ib)
        call gridlookup(bnum, bgrid, bval, ilocb, wb)
        bindex(ib) = ilocb
        bweight(ib) = wb
    end do


    do ib = 1,bnum
        bval = bgrid(ib)

        do ip = 1,nT
            bfval = (1.0_rk + r_b)*(1.0_rk - lambda(ip))*bval
            call gridlookup(bnum, bgrid, bfval, bindex1, bweight1)
            payindex(ib,ip) = bindex1
            payweight(ib,ip) = bweight1
        end do

        bfval = (1.0_rk + r_b)*bval
        call gridlookup(bnum, bgrid, bfval, bindex1, bweight1)
        nopayindex(ib) = bindex1
        nopayweight(ib) = bweight1
    end do

    end subroutine indexcal

 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! simulate the data in the model
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine simulate(chibshock, chibprob, nabil, eagridj, pabil, abilgrid, phid_abil, anum, bnum, yagg, probcost,  pi,   wh, wl,  psycost, &
        nsim, egridc, piec, piel, piec0, piel0, piea0, a0grid,edudecision, bfchoice, gl, gh, ghcollege, laborl, laborh, laborhcollege, &
        pchoicec, pchoicenc, agrid, bgrid, r_b, lambda, netphi, asim, bsim, idxprod, edusim, prodsim, paysim,       &
        idxability, abilsim,  iss,  idxchib)

    integer(ik):: nabil, anum, bnum,  ic, ia, ix,  nsim,  j, minplace(1), e,  &
        iss, olgo,  iab, t, ichib

    integer(ik):: pchoicec(jpaynum,nT,nabil, e_num, anum, bnum), pchoicenc(jpaynum,nT,nabil, e_num, anum, bnum), pop(jnum,3)

    integer(ik), allocatable::idxprodh(:,:), idxprodl(:,:), idxprod(:,:), idxiniwealth(:,:), ind_HK(:),  &
        idxcost(:,:,:), paysim(:,:), idxtemp(:,:),  idxability(:), wealthcut(:), abilitycut(:), idxchib(:)

    real(rk):: egridc(e_num), piec(e_num, e_num), piec0(e_num), piel(e_num, e_num), piel0(e_num), &
        a0grid(nabil,anum), edusim(jnum,nsim), edudecision(nabil,anum, ncost,chibnum),    &
        psycost(nabil,anum, ncost),  phidist(2), phigrid(2), bfchoice(nabil, anum, ncost,chibnum),  &
        gl(jnum,nT,nabil,e_num,anum,bnum), gh(jnum,nT,nabil,e_num,anum,bnum), laborl(jrnum,nT,nabil,e_num,anum,bnum),     &
        laborh(jrnum,nT,nabil,e_num,anum,bnum),    agrid(anum), bgrid(bnum),        &
        lambda(nT), r_b,   lclaborsim(jnum), lcbsim(jnum), lcasim(jnum), &
        edulclaborsim(jnum,3), edulcbsim(jnum,3), edulcasim(jnum,3), edulccsim(jnum,3), defaultsim(jnum,3), &
        dchoicesim(jnum,nsim), lcearnsim(jnum), edulcearnsim(jnum,3),  wh, wl, pi, &
        avglaborini(3,4), piea0(nabil, anum), piea0temp(anum),  probcost(nabil, anum, ncost), probcosttemp(ncost), yagg,  &
        netphi(nabil,anum), tuisim(jnum,nsim), ddebt(jnum,nsim), ddebtmu(jnum,nsim),  abilgrid(nabil),  &
        debtbyagesim(5),colrate_abil(3), colrate_wealth(4),abilitycutval(3), earncutval_abil(3), drateearn(4), &
        avgearnini(3,4), earnabil(3), earnpt(4), earncutval_wealth(4), wealthcutval(4),&
        eagridj(jnum,nabil,e_num,2), pabil(nabil), phid_abil(nabil),  avgability, colrate(3,4), colrate_d(3,4), abilratio, wealthratio, avgdebt,  avgpt, &
        stddebtholder, debtout, colrateagg, colrateaggd, debtpercent,ddebtpercentmu, ticagg, excagg, avgnettui, const, const1, totborrower, chibprob(chibnum), chibshock(chibnum), &
        ghcollege(jc-1_ik,nabil,anum,bnum,chibnum), laborhcollege(jc-1_ik,nabil,anum,bnum,chibnum)


    real(rk), allocatable::bsim(:,:), asim(:,:), laborsim(:,:), simprodl(:,:), simprodh(:,:), &
        simcost(:,:,:), iniwealth(:,:), prodsim(:,:), earnsim(:,:), simtemp(:,:), csim(:,:), &
        abilsim(:,:), simchib(:)

    integer(hid_t):: fileid
    character(30):: datafile,  datafile1, datestring, timestring

    intent(in):: nabil, anum, bnum, yagg,  probcost, nsim,  egridc, piec, piel, piec0, piel0,  a0grid, edudecision,  &
        psycost,   bfchoice,gl, gh, laborl, laborh,  pchoicec, pchoicenc, eagridj, phid_abil,  &
        agrid, bgrid, r_b, lambda,  wh, wl, pi,  piea0, netphi, iss, pabil, abilgrid, chibprob, chibshock, laborhcollege, ghcollege
    intent(out):: asim, bsim, idxprod, edusim, prodsim, paysim,  idxability, abilsim, idxchib

    call date_and_time(datestring, timestring)
    datafile = 'KK_stddebt_'
    olgo = len_trim(datafile)
    datafile(olgo + 1:olgo+4) = timestring(1:4)
    datafile1 = datafile(1:olgo+4)
    datafile1(olgo+5:olgo+13) = '_simul.h5'


    allocate(idxprodh(jnum, nsim), idxprodl(jnum, nsim), idxprod(jnum,  nsim), idxiniwealth(nabil, nsim), ind_HK(nsim), &
        idxcost(nabil, anum, nsim), paysim(jnum, nsim), idxability(nsim), abilsim(jnum, nsim), idxchib(nsim), simchib(nsim))
    allocate(bsim(jnum, nsim), asim(jnum, nsim), laborsim(jnum,nsim), simprodl(jnum, nsim), &
        simprodh(jnum, nsim), simcost(nabil, anum, nsim), iniwealth(nabil, nsim), prodsim(jnum, nsim), &
        earnsim(jnum, nsim), idxtemp(jnum,nsim), simtemp(jnum, nsim), csim(jnum, nsim))


    asim= 0.0_rk; bsim = 0.0_rk; paysim =1_ik; laborsim = 0.0_rk; lclaborsim= 0.0_rk; lcbsim= 0.0_rk; lcasim= 0.0_rk;
    edulclaborsim= 0.0_rk; edulcbsim= 0.0_rk;  edulcasim= 0.0_rk; pop =0_ik; defaultsim =0.0_rk; dchoicesim = 0.0_rk;
    earnsim = 0.0_rk; lcearnsim = 0.0_rk; edulcearnsim = 0.0_rk; avglaborini=0.0_rk;
    csim = 0.0_rk; edulccsim=0.0_rk; totborrower = 0.0_rk; tuisim = 0.0_rk

    !simulate ability
    call iniwealthsim(nabil, nsim, abilgrid, pabil, idxability, abilsim(1,1:nsim))

    !simulate the intial productivity - ! HK???
    e = 1_ik;
    call prodsimulate(e, jnum,  e_num, nsim, egridc, piec, piel, piec0, piel0,  idxprodh, simprodh, idxprodl, simprodl)

    !simulate initial wealth and education cost
    do iab = 1,nabil
        piea0temp(1:anum) =piea0(iab, 1:anum)/sum(piea0(iab, 1:anum))
        call iniwealthsim(anum, nsim, a0grid(iab, 1:anum), piea0temp, idxiniwealth(iab, 1:nsim), iniwealth(iab, 1:nsim))
    end do


    do iab = 1,nabil
        do ia = 1,anum
            probcosttemp(1:ncost) = probcost(iab,ia,1:ncost)/sum(probcost(iab,ia,1:ncost))
            if ( ncost .gt. 1_ik) then
                call iniwealthsim(ncost, nsim, psycost(iab,ia,1:ncost), probcosttemp, idxcost(iab,ia,1:nsim), simcost(iab,ia,1:nsim))
            else
                idxcost(iab,ia, 1:nsim)= ncost
            end if
        end do
    end do


    call iniwealthsim(chibnum, nsim, chibshock(1:chibnum), chibprob, idxchib(1:nsim), simchib(1:nsim))


    !decide education level
    do ic = 1,nsim
        iab = idxability(ic);
        ia = idxiniwealth(iab, ic);
        ix = idxcost(iab, ia, ic)
        ichib = idxchib(ic)

        edusim(1,ic) = edudecision(iab,ia, ix,ichib);
        abilsim(2:jnum, ic) = abilsim(1,ic)
        if (edusim(1,ic) .eq. 0.0_rk) then
            !they should draw initial productivity again from the non-college shock process
            idxprod(1:jnum, ic) = idxprodl(1:jnum,ic)
            prodsim(1:jnum, ic) = simprodl(1:jnum,ic)
            asim(1,ic) = iniwealth(iab,ic)
            edusim(2:jnum, ic) = edusim(1,ic)
            paysim(1,ic) = nT
        else !
            phigrid(1) = 0.5_rk; phigrid(2) = 1.0_rk
            idxprod(1:jnum, ic) = idxprodh(1:jnum, ic)
            prodsim(1:jnum, ic) = simprodh(1:jnum, ic)
            bsim(1,ic) = bfchoice(iab,ia,ix,ichib)
            if (bsim(1,ic) .gt. -1.0E-003_rk) then
                bsim(1,ic) = 0.0_rk
            end if
            asim(1,ic) = iniwealth(iab,ic)
            edusim(2:jnum, ic) = edusim(1,ic)
            !simulation for drop-out
            phidist =0.0_rk; phidist(1) = phid_abil(iab); phidist(2) = 1.0_rk-phid_abil(iab);            
            call iniwealthsim(2_ik, 1, phigrid, phidist, minplace(1), edusim(jd,ic)) ! drop out case
            edusim(jd+1:jnum,ic) = edusim(jd,ic)
            paysim(1,ic) = 1_ik
        end if
       
    end do

    ! HK??
    do ic =1,nsim
        if (edusim(jd,ic) .eq. 0.5_rk) then
            e = 0_ik
            call prodsimulate(e,  jnum-jd+1,  e_num, 1, egridc, piec, piel, piec0, piel0,  idxtemp, simtemp, &
                idxprod(jd:jnum,ic),prodsim(jd:jnum,ic))
        end if
    end do

    call simul_sample(anum, bnum, nabil, nsim, wh, wl, r_b, idxprod, idxability, idxchib, pchoicenc, pchoicec, simchib,      &
                     netphi, a0grid, agrid, bgrid, lambda, eagridj, ghcollege, gh, gl, laborhcollege, laborh, laborl,                      &
                     laborsim, dchoicesim, csim, earnsim, asim, bsim, paysim, edusim, tuisim, ddebt, ddebtmu)
    
    pop=0_ik
    do j = 1,jnum
        lclaborsim(j) = sum(laborsim(j,1:nsim))/nsim
        lcbsim(j) = sum(bsim(j,1:nsim))/nsim
        lcasim(j) = sum(asim(j,1:nsim))/nsim
        lcearnsim(j) = sum(earnsim(j,1:nsim))/nsim
        do ic = 1,nsim
            if (edusim(j,ic) .eq. 0.0_rk) then
                edulclaborsim(j,1) = edulclaborsim(j,1) +laborsim(j,ic)
                edulcbsim(j,1) = edulcbsim(j,1) + bsim(j,ic)
                edulcasim(j,1) = edulcasim(j,1) + asim(j,ic)
                edulcearnsim(j,1) =edulcearnsim(j,1)+earnsim(j,ic)
                edulccsim(j,1) = edulccsim(j,1) + csim(j, ic)
                pop(j,1) = pop(j,1)+ 1_ik
            elseif( edusim(j,ic) .eq. 0.5_rk) then
                edulclaborsim(j,2) = edulclaborsim(j,2) +laborsim(j,ic)
                edulcbsim(j,2) = edulcbsim(j,2) + bsim(j,ic)
                edulcasim(j,2) = edulcasim(j,2) + asim(j,ic)
                defaultsim(j,2) = defaultsim(j,2) +dchoicesim(j,ic)
                edulcearnsim(j,2) =edulcearnsim(j,2)+earnsim(j,ic)
                edulccsim(j,2) = edulccsim(j,2) + csim(j, ic)
                pop(j,2) = pop(j,2)+ 1_ik
            else
                edulclaborsim(j,3) = edulclaborsim(j,3) +laborsim(j,ic)
                edulcbsim(j,3) = edulcbsim(j,3) + bsim(j,ic)
                edulcasim(j,3) = edulcasim(j,3) + asim(j,ic)
                defaultsim(j,3) = defaultsim(j,3) +dchoicesim(j,ic)
                edulcearnsim(j,3) =edulcearnsim(j,3)+earnsim(j,ic)
                edulccsim(j,3) = edulccsim(j,3) + csim(j, ic)
                pop(j,3) = pop(j,3)+ 1_ik
            end if
        end do

    end do
    edulclaborsim = edulclaborsim/pop
    edulcbsim = edulcbsim/pop
    edulcasim =edulcasim/pop
    edulccsim =edulccsim/pop
    edulcearnsim = edulcearnsim/pop
    defaultsim = defaultsim/pop
    allocate(wealthcut(nsim), abilitycut(nsim))

    t = 0_ik !steady states
    call sim_moments(phid_abil, totborrower, t, iss,  nabil, nsim, yagg, ddebt, ddebtmu,                      &
        tuisim, bsim, asim, laborsim, earnsim, idxability, abilsim, edusim, dchoicesim,     &
        abilitycut, wealthcut, earncutval_abil, wealthcutval, earncutval_wealth, avglaborini, avgearnini, earnpt,    &
        earnabil, avgability, drateearn, debtbyagesim, colrate, colrate_d, colrate_abil, colrate_wealth, abilratio, wealthratio,       &
        avgdebt,  avgpt, stddebtholder, debtout, colrateagg, colrateaggd, debtpercent,ddebtpercentmu, ticagg, excagg, avgnettui, const, const1, abilitycutval  )

    if (iarraysim .eq. 1_ik) then
        call hdf5_openf(datafile1, fileid)
        call hdf5_write(prodsim(1,1:nsim), fileid, 'prodsim')
        call hdf5_write(edusim(1,1:nsim), fileid, 'edusim')
        call hdf5_write(bsim(1,1:nsim), fileid, 'bsim')
        call hdf5_write(asim(1,1:nsim), fileid, 'asim')
        call hdf5_write(tuisim(1,1:nsim), fileid, 'tuisim')
        call hdf5_write(laborsim(1,1:nsim), fileid, 'laborsim')
        call hdf5_write(tau, fileid, 'tau')
        call hdf5_write(wl, fileid, 'wl')
        call hdf5_write(pi, fileid, 'pi')
        call hdf5_write(egridc, fileid, 'egridc')
        call hdf5_write(abilitycut, fileid, 'abilitycut')
        call hdf5_write(wealthcut, fileid, 'wealthcut')
        call hdf5_write(earnsim(1,1:nsim), fileid, 'earnsim')
        call hdf5_write(netphi, fileid, 'netphi')
        call hdf5_write(a0grid, fileid, 'a0grid')
        call hdf5_write(piea0, fileid, 'piea0')
        call hdf5_write(laborh(1,1,1:nabil,1,1:anum,1) , fileid, 'laborhini')
        call hdf5_write(bfchoice, fileid, 'bfchoice')
        call hdf5_write(abilitycutval, fileid, 'abilitycutval')
        call hdf5_write(earncutval_abil, fileid, 'earncutval_abil')
        call hdf5_write(wealthcutval, fileid, 'wealthcutval')
        call hdf5_write(earncutval_wealth, fileid, 'earncutval_wealth')
        call hdf5_write(avglaborini, fileid, 'avglaborini')
        call hdf5_write(avgearnini, fileid, 'avgearnini')
        call hdf5_write(earnpt, fileid, 'earnpt')
        call hdf5_write(earnabil, fileid, 'earnabil')
        call hdf5_write(yagg, fileid, 'yagg')
        call hdf5_write(gdp_data, fileid, 'gdp_data')
        call hdf5_write(bbar, fileid, 'bbar')
        call hdf5_write(ncost, fileid, 'ncost')
        call hdf5_write(chibnum, fileid, 'chibnum')
        call hdf5_write(abilsim, fileid, 'abilsim')

        call hdf5_closef(fileid)
    end if

    end subroutine simulate
   
    
    
    subroutine simul_sample(anum, bnum, nabil, nsim, wh, wl, r_b, idxprod, idxability, idxchib, pchoicenc, pchoicec, simchib,    &
                          netphi, a0grid, agrid, bgrid, lambda, eagridj, ghcollege, gh, gl, laborhcollege, laborh, laborl,&
                          laborsim, dchoicesim, csim, earnsim, asim, bsim, paysim, edusim, tuisim, ddebt, ddebtmu)

    integer(ik):: ic, j, nsim, ie, ip, iab, ichib, ipf, anum, bnum, nabil
    integer(ik):: idxprod(jnum,nsim), paysim(jnum,nsim), idxability(nsim), idxchib(nsim), pchoicec(jpaynum,nT,nabil, e_num, anum, bnum), &
                  pchoicenc(jpaynum,nT,nabil, e_num, anum, bnum), pchoicecin(anum,bnum), pchoicencin(anum,bnum)
    
    real(rk):: aval, bval, elevel, chib, wh, wl, r_b, afval, bfval
    real(rk):: asim(jnum,nsim), bsim(jnum,nsim), edusim(jnum,nsim), simchib(nsim), ghtemp(anum,bnum), lhtemp(anum,bnum),                 &
               lltemp(anum,bnum), ghcollege(jc-1_ik,nabil,anum,bnum,chibnum), laborhcollege(jc-1_ik,nabil,anum,bnum,chibnum),            &
               gh(jnum,nT,nabil,e_num,anum,bnum), laborl(jrnum,nT,nabil,e_num,anum,bnum), laborh(jrnum,nT,nabil,e_num,anum,bnum),        &
               gl(jnum,nT,nabil,e_num,anum,bnum), ddebt(jnum,nsim), ddebtmu(jnum,nsim), lambda(nT), a0grid(nabil,anum),                  &
               agrid(anum), bgrid(bnum), tuisim(jnum,nsim), eagridj(jnum,nabil,e_num,2), earnsim(jnum,nsim), laborsim(jnum,nsim),        &
               dchoicesim(jnum,nsim), csim(jnum,nsim), netphi(nabil,anum)


    intent(in):: anum, bnum, nabil, nsim, wh, wl, r_b, idxprod, idxability, idxchib, pchoicenc, pchoicec, simchib,               &
                 a0grid, agrid, bgrid, lambda, eagridj, ghcollege, gh, gl, laborhcollege, laborh, laborl, netphi

    intent(out):: laborsim, dchoicesim, csim, earnsim, tuisim
    
    intent(inout)::  asim, bsim, paysim, ddebt, ddebtmu, edusim

    
    ddebt = 0.0_rk; ddebtmu = 0.0_rk;

    do ic = 1,nsim
        do j= 1,jnum
            aval = asim(j,ic); bval = bsim(j,ic); elevel = edusim(j,ic)
            ie = idxprod(j, ic); ip = paysim(j,ic); iab = idxability(ic)
            ichib =  idxchib(ic); chib = simchib(ichib)

            if (elevel.eq.1.0_rk) then
                if (j.le.jc-1_ik) then
                    ghtemp = ghcollege(j,iab,:,:,ichib)
                    lhtemp = laborhcollege(j,iab,:,:,ichib)
                    pchoicecin =1_ik!
                elseif (j.le.jpaynum) then
                    ghtemp = gh(j,ip,iab,ie,:,:)
                    lhtemp = laborh(j,ip,iab,ie,:,:)
                    pchoicecin = pchoicec(j,ip,iab,ie,:,:)
                elseif (j.le.jrnum) then
                    ghtemp = gh(j,ip,iab,ie,:,:)
                    lhtemp = laborh(j,ip,iab,ie,:,:)
                    pchoicecin =1_ik
                else
                    ghtemp = gh(j,ip,iab,ie,:,:)
                    lhtemp = 0.0_rk
                    pchoicecin = 1_ik
                end if               
            else
                
                if (j.le.jpaynum) then
                    lltemp = laborl(j,ip,iab,ie,:,:)
                    pchoicencin = pchoicenc(j,ip,iab,ie,:,:)
                elseif (j.le.jrnum) then
                    lltemp = laborl(j,ip,iab,ie,:,:)
                    pchoicencin = 1_ik
                else
                    lltemp = 0.0_rk
                    pchoicencin = 1_ik
                end if
            end if



            call simcal( j, ip, anum, bnum, netphi(iab,:), pchoicencin, pchoicecin, chib, aval, bval, elevel, wl, wh, r_b, &
                lambda(ip), a0grid(iab,:), agrid, bgrid, tuisim(j,ic), gl(j,ip,iab,ie,:,:), lltemp, ghtemp, &
                lhtemp, eagridj(j,iab,ie,:), ipf, earnsim(j,ic), laborsim(j,ic), afval, dchoicesim(j,ic),  &
                csim(j,ic), bfval, ddebt(j,ic), ddebtmu(j,ic))

            if ( j.lt. jnum ) then
                asim(j+1,ic) = afval
                paysim(j+1,ic) =ipf!1_ik
                bsim(j+1,ic) = bfval
            end if
        end do
    end do

    
    end subroutine simul_sample   

    

    subroutine simcal( j, ip, anum, bnum, netphi, pchoicenc, pchoicec, chib, aval, bval, elevel, wl, wh, r_b, lambda, a0grid, agrid, bgrid, &
        tuisim, gl, laborl, gh, laborh, eagridj, paysim, earnsim, nval, afval,  dchoicesim, csim,  &
        bsim, ddebt, ddebtmu)
    integer(ik)             :: iloca, ilocb
    integer(ik)   , intent(in) :: j, ip, anum, bnum, pchoicenc(anum,bnum), pchoicec(anum,bnum)
    integer(ik)   , intent(out):: paysim
    real(rk)                :: wa, wb, paychoice
    real(rk) ,    intent(in) :: aval, bval, elevel, wl, wh, r_b, lambda, a0grid(anum), agrid(anum), bgrid(bnum), &
       gl(anum,bnum), laborl(anum,bnum), gh(anum,bnum), laborh(anum,bnum), eagridj(2), chib, netphi(anum)
    real(rk) ,    intent(out):: earnsim, nval, afval, dchoicesim, csim,  bsim,  tuisim
    real(rk) ,    intent(out):: ddebt, ddebtmu
    

    if ( elevel .eq. 0.0_rk) then

        if (j .eq. 1_ik) then
            call gridlookup(anum, a0grid, aval, iloca, wa)
        else
            call gridlookup(anum, agrid, aval, iloca, wa)
        end if

        afval = wa*gl(iloca, bnum) +(1.0_rk-wa)*gl(iloca+1,bnum)
        if( j.le. jrnum) then
            nval = wa*laborl(iloca, bnum) +(1.0_rk-wa)*laborl(iloca+1,bnum)
            earnsim = wl*eagridj(2)*nval
        else
            nval = 0.0_rk
            earnsim = wl*ssrr*eagridj(2)
        end if

        paysim = nT
        dchoicesim=0.0_rk
        ddebt = 0.0_rk
        ddebtmu = 0.0_rk
        bsim = 0.0_rk
        if ( j .eq. 1_ik) then
            csim = (1.0_rk-tau)*earnsim +aval - afval;
        else
            csim = (1.0_rk-tau)*earnsim+ (1+r)*aval-afval;
        end if

    elseif (elevel .eq. 0.5_rk) then ! this only happense from age jD

        call gridlookup(anum, agrid, aval, iloca, wa)
        call gridlookup(bnum, bgrid, bval, ilocb, wb)

        afval = wb*(wa*gl(iloca, ilocb) +(1.0_rk-wa)*gl(iloca+1,ilocb))+(1.0_rk-wb)*(wa*gl(iloca, ilocb+1) +(1.0_rk-wa)*gl(iloca+1,ilocb+1))

        if ( j .le. jrnum) then
            nval = wb*(wa*laborl(iloca, ilocb) +(1.0_rk-wa)*laborl(iloca+1,ilocb))+(1.0_rk-wb)*(wa*laborl(iloca, ilocb+1) +(1.0_rk-wa)*laborl(iloca+1,ilocb+1))
            earnsim = wl*eagridj(2)*nval
        else
            nval = 0.0_rk
            earnsim = wl*ssrr*eagridj(2)
        end if

        if ( j .le. jpaynum .and. bval .lt. 0.0_rk) then
            if ( ilocb .eq. bnum-1) then
                wb = 1.0_rk
            end if
            paychoice = wb*(wa*(dble(pchoicenc(iloca, ilocb))) +(1.0_rk-wa)*(dble(pchoicenc(iloca+1,ilocb))))+(1.0_rk-wb)*(wa*(dble(pchoicenc(iloca, ilocb+1))) +(1.0_rk-wa)*(dble(pchoicenc(iloca+1,ilocb+1))))
        else
            paychoice = 1.0_rk
        end if

        if ( paychoice .ge. 0.5_rk .or. j .gt. jpaynum) then
            paysim = min(ip+1_ik, nT)
            if ( paysim .eq. nT) then
                bsim=0.0_rk
            else
                bsim = (1.0_rk + r_b)*(1.0_rk - lambda)*bval
            end if
            csim= (1+r)*aval + (1.0_rk-tau)*earnsim +lambda*bval - afval;
            dchoicesim=0.0_rk
            ddebt = 0.0_rk
            ddebtmu = 0.0_rk
        else
            paysim = ip
            bsim =  (1.0_rk + r_b)*bval
            dchoicesim =1.0_rk
            ddebt = lambda*bval
            ddebtmu = 1.0_rk
            csim = (1+r)*aval + (1.0_rk-tau)*(1.0_rk-wgar)*earnsim- afval;
            earnsim = (1.0_rk-wgar)*earnsim
        end if

    else !college educated

        call gridlookup(bnum, bgrid, bval, ilocb, wb)

        if (j .lt. jc) then
            call gridlookup(anum, a0grid, aval, iloca, wa)
            afval = aval
            if ( j .eq. jc-1) then
                afval = wb*(wa*gh(iloca, ilocb) +(1.0_rk-wa)*gh(iloca+1,ilocb))+(1.0_rk-wb)*(wa*gh(iloca, ilocb+1) +(1.0_rk-wa)*gh(iloca+1,ilocb+1))
            end if
            tuisim = wa*netphi(iloca)+(1.0_rk-wa)*netphi(iloca+1)
        else
            call gridlookup(anum, agrid, aval, iloca, wa)
            afval = wb*(wa*gh(iloca, ilocb) +(1.0_rk-wa)*gh(iloca+1,ilocb))+(1.0_rk-wb)*(wa*gh(iloca, ilocb+1) +(1.0_rk-wa)*gh(iloca+1,ilocb+1))
        end if

        if ( j .le. jrnum) then
            if (laborh(iloca,ilocb) .eq. 0.0_rk) then
                wb = 0.0_rk;
            end if

            nval = wb*(wa*laborh(iloca, ilocb) +(1.0_rk-wa)*laborh(iloca+1,ilocb))+(1.0_rk-wb)*(wa*laborh(iloca, ilocb+1) +(1.0_rk-wa)*laborh(iloca+1,ilocb+1))
        else
            nval = 0.0_rk
        end if
        
        if ( j .lt. jc) then
            paysim = ip
            bsim = bval
            earnsim =wl*eagridj(1)*nval
            csim= aval + (1.0_rk-tau)*earnsim -tuisim - bval/(jc-1)+chib*bval
            dchoicesim=0.0_rk
        else
            if ( j .le. jrnum) then
                earnsim= wh*eagridj(1)*nval
            else
                earnsim= wh*ssrr*eagridj(1)
            end if

            if ( j .le. jpaynum .and. bval .lt. 0.0_rk) then
                if ( ilocb .eq. bnum-1) then
                    wb = 1.0_rk
                end if
                paychoice = wb*(wa*(dble(pchoicec(iloca, ilocb))) +(1.0_rk-wa)*(dble(pchoicec(iloca+1,ilocb))))+(1.0_rk-wb)*(wa*(dble(pchoicec(iloca, ilocb+1))) +(1.0_rk-wa)*(dble(pchoicec(iloca+1,ilocb+1))))
            else
                paychoice =1.0_rk
            end if

            if ( paychoice .ge. 0.5_rk .or. j .gt. jpaynum) then
                paysim = min(ip +1_ik, nT)
                if ( paysim.eq. nT) then
                    bsim =0.0_rk
                else
                    bsim = (1.0_rk + r_b)*(1.0_rk - lambda)*bval
                end if
                csim= (1+r)*aval + (1.0_rk-tau)*earnsim +lambda*bval - afval;
                dchoicesim=0.0_rk                
                ddebt = 0.0_rk
                ddebtmu = 0.0_rk
            else
                paysim = ip
                bsim=  (1.0_rk + r_b)*bval
                dchoicesim =1.0_rk
                ddebt = lambda*bval
                ddebtmu = 1.0_rk
                csim = (1+r)*aval + (1.0_rk-tau)*(1.0_rk-wgar)*earnsim- afval;
                earnsim = (1.0_rk-wgar)*earnsim                
            end if
        end if
    end if

    end subroutine simcal


    subroutine sim_moments(phid_abil, totborrower, t, iss, nabil, nsim,  yagg, ddebt, ddebtmu,                              &
        tuisim, bsim, asim, laborsim, earnsim, idxability, abilsim, edusim, dchoicesim, &
        abilitycut, wealthcut, earncutval_abil, wealthcutval, earncutval_wealth, avglaborini, avgearnini, earnpt, &
        earnabil, avgability, drateearn, debtbyagesim, colrate_ini, colrate_ini_d, colrate_abil_ini, colrate_wealth_ini,abilratio, wealthratio,     &
        avgdebt,  avgpt, stddebtholder, debtout, colrateagg, colrateaggd, debtpercent,ddebtpercentmu, ticagg, excagg, avgnettui, const, const1, abilitycutval  )

    integer(ik)::  ic, ia,  nsim, i, j, ncut, ind, k, k2, nabil, iss,  ng, count, iab, t, jin
    integer(ik):: pop(jnum,3), wealthcut(nsim), abilitycut(nsim), nyearpay(nsim), idxability(nsim)

    integer(ik), allocatable::ind_HK(:), indearn(:), earncut(:)

    real(rk)::  yagg, edusim(jnum,nsim), colrate_ini(3,4), colrateagg, total_ini_d(3,4), dchoicesim(jnum,nsim), colrate_ini_d(3,4),&
        avglaborini(3,4),avgasimini(3,4),avgbsimini(3,4), debtbyagesim(5), laborabil(3), laborpt(4),  &
        stdloanabil(3), stdloanpt(4), tuisim(jnum,nsim), avgnetphi(3,4), tuipt(4), tuiabil(3), temp3(3,4), &
        colrate_abil_ini(3), colrate_wealth_ini(4), abilratio, wealthratio, avgdebt, avgdebt_c, avgdebt_d,stddebtholder, total_bor, total_borj1,  &
        ddebt(jnum,nsim), ddebtmu(jnum,nsim),  debtout, avgnettui,  tic(3,4), ticagg, tic_abil(3), tic_wealth(4), tic_abil_wealth(3,4),  &
        exc(3,4), exc_abil(3), exc_wealth(4), exc_abil_wealth(3,4) , excagg, abilitycutval(3), earncutval_abil(3), &
        avgearnini(3,4), earnabil(3), earnpt(4), earncutval_wealth(4), wealthcutval(4), phid_abil(nabil), &
        const, drateearn(4), dsample(4), bearnq(4), bsim(jnum, nsim), asim(jnum, nsim), laborsim(jnum,nsim),  &
        earnsim(jnum, nsim), abilsim(jnum, nsim), avgability,debtpercent, count1,ddebtpercentmu, totborrower, borrower(3,4),  &
        drate, drateearn1(4), avglaborhini, phid, educated, noneducated, colrateaggd, total, total_col, const1, colh, cold, &
        drateh, drated, count1_h, count1_d, avgpt, avgpt_c, avgpt_nc, avgpttotal, avgpttotal_c, avgpttotal_nc,              &
        avgearn, earn_w, earn_wc, earn_wnc, earn_wd, count_w,  count_wc,  count_wnc,  count_wd, cdrop_count(4), c_count(4), cdrop_earn(4)

    real(rk), allocatable:: assetHK(:), prodHK(:),   earning(:), default(:), bearn(:), bsimtemp(:), edusimtemp(:)

    intent(in):: phid_abil, nabil, nsim, iss,  yagg, bsim, asim, laborsim,  earnsim,idxability, abilsim, edusim, dchoicesim, ddebt, ddebtmu, tuisim, t
    intent(out):: abilitycut, wealthcut, earncutval_abil, wealthcutval, earncutval_wealth, avglaborini, avgearnini, &
        earnpt, earnabil, avgability, drateearn, debtbyagesim, colrate_ini, colrate_ini_d, abilratio, wealthratio, avgdebt,  avgpt, &
        colrate_abil_ini, colrate_wealth_ini, stddebtholder, debtout, colrateagg, colrateaggd, debtpercent,ddebtpercentmu, ticagg, excagg, &
        avgnettui, const, totborrower, abilitycutval, const1

    if ( t.eq. 0_ik) then
        jin = jd-1
    else
        jin = 1_ik
    end if
    
    allocate(ind_HK(nsim))
    pop =0_ik;  avglaborini=0.0_rk; avgasimini=0.0_rk; avgbsimini=0.0_rk; avgnetphi=0.0_rk;
    tuipt=0.0_rk; tuiabil=0.0_rk; nyearpay = 0.0_rk;

    write(*,*) 'sim_moments starts'
    !calculate deafult rate over earnings
    count = 0
    do ic = 1, nsim
        do j = jc, jrnum
            count = count +1
        end do
    end do
    ng = count

    allocate(earning(ng), default(ng), indearn(ng), earncut(ng), bearn(ng), bsimtemp(ng), edusimtemp(ng))
    earning = 0.0_rk; default = 0.0_rk; indearn = 0_ik; earncut = 0_ik; bearn = 0.0_rk; bsimtemp = 0.0_rk; edusimtemp = 0.0_rk   
    count = 0
    do ic = 1, nsim
        do j = jc, jrnum
            count = count +1
            earning(count) = earnsim(j,ic)
            default(count) = dchoicesim(j,ic)
            
            if ( edusim(j,ic) .eq. 1.0_rk) then   
            bearn(count) = bsim(jc,ic)
            elseif ( edusim(j,ic) .eq. 0.5_rk) then
            bearn(count) = bsim(jd,ic)
            else
            bearn(count) = 0.0_rk
            end if
            bsimtemp(count) = bsim(j,ic)
            edusimtemp(count) = edusim(j,ic)
        end do
    end do
    call qsortd(earning, indearn, ng)

    ncut = int(ng/4);
    do k =1,4
        if (k.eq. 4_ik) then
            k2 = ng
        else
            k2 = k*ncut
        end if
        do i = (k-1)*ncut+1, k2
            ind = indearn(i)
            earncut(ind)= k
        end do
    end do

    drateearn = 0.0_rk; dsample = 0.0_rk; bearnq =0.0_rk; drateearn1 = 0.0_rk
    c_count = 0.0_rk; cdrop_count = 0.0_rk; cdrop_earn = 0.0_rk
    do ic = 1, ng
        k = earncut(ic)
        if (edusimtemp(ic).ge.0.5_rk  .and. bsimtemp(ic) .lt. 0.0_rk) then
            drateearn1(k) = drateearn1(k) + default(ic)
            dsample(k) = dsample(k) +1.0_rk
            bearnq(k) = bearnq(k) + bearn(ic)            
        end if
        if (edusimtemp(ic).ge.0.5_rk) then !college goers
            c_count(k) = c_count(k) + 1.0_rk
            if (edusimtemp(ic).eq.0.5_rk) then !college dropouts
                cdrop_count(k) = cdrop_count(k) + 1.0_rk
            end if
        end if      
    end do

    drate = sum(drateearn1)/sum(dsample)
    drateearn = drateearn1/dsample
    bearnq = (bearnq/yagg)*gdp_data
    bearnq= bearnq/dsample
    cdrop_earn = cdrop_count/c_count

    !calculate the college enrollment rate over initial wealth distribution and intial persistent shock
    allocate(assetHK(nsim), prodHK(nsim))

    call qsortd(asim(1,1:nsim),ind_HK, nsim)
    assetHK(1:nsim) = asim(1,ind_HK)
    ncut = int(nsim/4);
    earncutval_wealth = 9999999.0_rk; earncutval_abil = 9999999.0_rk
    do k =1,4
        if (k.eq. 4_ik) then
            k2 = nsim
        else
            k2 = k*ncut
        end if
        do i = (k-1)*ncut+1, k2
            ind = ind_HK(i)
            wealthcut(ind)= k
            earncutval_wealth(k) = dmin1(earncutval_wealth(k),earnsim(1,ind))
            if ( i .eq. k2) then
                wealthcutval(k) = asim(1,ind)
            end if
        end do
    end do
    write(*,*) 'wealthcutval', wealthcutval
    call qsortd(abilsim(1,1:nsim),ind_HK, nsim)
    ncut = int(nsim/3);
    do k =1,3
        if (k.eq. 3_ik) then
            k2 = nsim
        else
            k2 = k*ncut
        end if
        do i = (k-1)*ncut+1, k2
            ind = ind_HK(i)
            abilitycut(ind)= k
            earncutval_abil(k) = dmin1(earncutval_abil(k),earnsim(1,ind))
            if ( i .eq. k2) then
                abilitycutval(k) = abilsim(1,ind)
            end if

        end do
    end do

    write(*,*) 'abilitycutval', abilitycutval 
    !calculate college enrollment rate over ability and parental wealth
    noneducated = 0.0_rk; educated = 0.0_rk
    colrateagg = 0.0_rk; colrateaggd=0.0_rk; stddebtholder = 0.0_rk;
    debtbyagesim = 0.0_rk; avgnettui = 0.0_rk;
    temp3 = 0.0_rk; avgdebt =0.0_rk;  
    avgability = 0.0_rk; totborrower = 0.0_rk; debtout = 0.0_rk; 
    colrate_ini=0.0_rk; colrate_ini_d =0.0_rk; colrate_abil_ini=0.0_rk; colrate_wealth_ini = 0.0_rk; total_ini_d = 0.0_rk; 
    debtpercent = 0.0_rk; ddebtpercentmu=0.0_rk;
    tic = 0.0_rk; ticagg = 0.0_rk; tic_abil = 0.0_rk; tic_wealth = 0.0_rk; tic_abil_wealth = 0.0_rk; 
    exc = 0.0_rk; excagg = 0.0_rk; exc_abil = 0.0_rk; exc_wealth = 0.0_rk; exc_abil_wealth = 0.0_rk; 
    avgearnini = 0.0_rk; const = 0.0_rk; const1 =0.0_rk      
    borrower =0.0_rk; colh = 0.0_rk; cold = 0.0_rk 
    avgpt = 0.0_rk;  avgpt_c = 0.0_rk; avgpt_nc = 0.0_rk;
    total_bor = 0.0_rk; count1=0.0_rk; total_borj1 = 0.0_rk;
    total_col = 0.0_rk; 
    total=0.0_rk; 
    avgpttotal = 0.0_rk;avgpttotal_c = 0.0_rk;avgpttotal_nc = 0.0_rk;

    avgearn = 0.0_rk; earn_w= 0.0_rk;  earn_wc= 0.0_rk;  earn_wnc= 0.0_rk;  earn_wd= 0.0_rk; 
    avglaborhini = 0.0_rk; count_w= 0.0_rk;  count_wc= 0.0_rk;  count_wnc= 0.0_rk;  count_wd= 0.0_rk; 
    avgdebt_c =0.0_rk; avgdebt_d =0.0_rk; colh = 0.0_rk; cold = 0.0_rk 
    drateh = 0.0_rk; drated = 0.0_rk; count1_h = 0.0_rk; count1_d= 0.0_rk
    do ic = 1, nsim

        iab  = abilitycut(ic); ia = wealthcut(ic)

            !Moments across distribution using j = 1:jd-1  
            do j =1, jin
            total_ini_d(iab,ia) = total_ini_d(iab,ia) + 1.0_rk
                if (edusim(j,ic) .eq. 1.0_rk) then
                    phid  = phid_abil(idxability(ic))
                    colrate_ini_d(iab,ia) = colrate_ini_d(iab,ia) +1.0_rk-phid
                    colrate_ini(iab,ia) = colrate_ini(iab,ia) + 1.0_rk
                    colrate_abil_ini(iab) = colrate_abil_ini(iab) + 1.0_rk-phid
                    colrate_wealth_ini(ia) = colrate_wealth_ini(ia)+1.0_rk-phid

                    avgnetphi(iab,ia) = avgnetphi(iab,ia) + tuisim(j,ic)
                    avgearnini(iab,ia) = avgearnini(iab,ia) + earnsim(j,ic)
                    !if ( bsim(j,ic) .lt. 0.0_rk) then
                    avgbsimini(iab,ia) = avgbsimini(iab,ia) + bsim(j,ic)
                    borrower(iab,ia) = borrower(iab,ia) +1.0_rk
                    !end if
                    avglaborini(iab,ia) = avglaborini(iab,ia) + laborsim(j,ic)
            
                    ! Initial Parental Tranfers distribution using j = 1:jc-1                
                    avgasimini(iab,ia) = avgasimini(iab,ia) + asim(j,ic)
                    temp3(iab,ia) = temp3(iab,ia) + 1.0_rk
                end if
            end do
            
            do j =1, jd-1_ik
                !Aggregate count using j = 1, jc-1
                total = total + 1.0_rk
                if (edusim(j,ic) .eq. 1.0_rk) then !including dropout
                    phid  = phid_abil(idxability(ic))
                    colrateaggd = colrateaggd + 1.0_rk-phid     
                end if
            end do
            
                    
            do j = 1, jc-1_ik

                if (edusim(j,ic) .eq. 0.0_rk) then
                    if (j.eq.1_ik) then
                        avgpt_nc = avgpt_nc + asim(j,ic)
                        avgpttotal_nc = avgpttotal_nc + 1.0_rk
                    end if
                elseif (edusim(j,ic) .eq. 1.0_rk) then
                    avgpt_c = avgpt_c + asim(j,ic)
                    avgpttotal_c = avgpttotal_c + 1.0_rk

                end if

                if (edusim(j,ic) .eq. 1.0_rk) then !including dropout
                    total_col = total_col + 1.0_rk
                    !Calculating aggregate moments using age from 1 to jc-1

                    avgnettui = avgnettui + tuisim(j,ic)
                    avgability = avgability + abilsim(j, ic)
                    avgdebt = avgdebt + bsim(j,ic)


                    if (bsim(j,ic) .le. dmax1(bbar + 1.0E-003_rk, -fam*tuisim(j,ic)*(dble(jc)-1.0_rk)+1.0E-003_rk )) then
                        const = const + 1.0_rk
                    end if

                    !fraction of constrained people(max(tic,exc) at j=1:jc-1)
                    if (j.le.jin) then
                        if (bsim(j,ic) .le. dmax1(bbar + 1.0E-003_rk, -fam*tuisim(j,ic)*(dble(jc)-1.0_rk)+1.0E-003_rk )) then
                            const1 = const1 + 1.0_rk
                        end if
                        if ( bsim(j,ic) .le. -fam*tuisim(j,ic)*(dble(jc)-1.0_rk) +1.0E-003_rk) then
                            tic(iab,ia) = tic(iab,ia) +1.0_rk
                            tic_abil(iab) = tic_abil(iab) + 1.0_rk
                            tic_wealth(ia) = tic_wealth(ia) + 1.0_rk
                        end if
                        if ( bsim(j,ic) .le. bbar+1.0E-003_rk) then
                            exc(iab,ia) = exc(iab,ia) +1.0_rk
                            exc_abil(iab) = exc_abil(iab) + 1.0_rk
                            exc_wealth(ia) = exc_wealth(ia) + 1.0_rk
                        end if
                    end if

                    !Calculating student debt holder ratio
                    if (bsim(j,ic) .lt. 0.0_rk) then
                        total_bor = total_bor + 1.0_rk
                        if (j.le.jin) then
                            total_borj1 = total_borj1 + 1.0_rk
                        end if
                    end if

                    avgearn = avgearn + earnsim(j,ic)
                    avglaborhini = avglaborhini+ laborsim(j,ic)
                end if
            end do

        ! calculate non-educated between jc and jnum     

        if ( edusim(jc,ic) .eq. 1.0_rk) then
            avgdebt_c = avgdebt_c + bsim(jc,ic)
            colh = colh + 1.0_rk
        elseif ( edusim(jd,ic) .eq. 0.5_rk) then
            avgdebt_d = avgdebt_d + bsim(jd,ic)
            cold = cold +1.0_rk
        end if

        
        do j= jc,jnum

            if (edusim(j,ic) .le. 0.5_rk) then !including those who dropped out or high school graduates
                noneducated = noneducated + 1.0_rk
            else
                educated = educated + 1.0_rk
            end if
        end do        

        do j= jd, jrnum

             if ( edusim(j,ic) .ge. 0.5_rk) then
                 if (bsim(j,ic).lt.0.0_rk  ) then
                     count1 = count1 + 1.0_rk
                 end if
             end if

             if ( edusim(j,ic) .eq. 0.5_rk .and. bsim(j,ic).lt.0.0_rk) then
                count1_d = count1_d + 1.0_rk
                if (dchoicesim(j,ic) .eq. 1.0_rk) then
                    drated = drated +1.0_rk
                end if
             end if
        end do
            
        !calculate labor income from age 25-64
        do j= jc, jrnum
             earn_w= earn_w + earnsim(j,ic) 
             count_w = count_w + 1.0_rk            
              if ( edusim(j,ic) .eq. 0.5_rk) then
                 earn_wd= earn_wd + earnsim(j,ic) 
                 count_wd = count_wd + 1.0_rk           
             elseif ( edusim(j,ic) .eq. 1.0_rk) then
                 earn_wc= earn_wc + earnsim(j,ic) 
                 count_wc = count_wc + 1.0_rk
                 if (bsim(j,ic).lt.0.0_rk  ) then
                     count1_h = count1_h + 1.0_rk
                     if (dchoicesim(j,ic) .eq. 1.0_rk) then
                         drateh = drateh +1.0_rk
                     end if
                 end if
             else
                 earn_wnc= earn_wnc + earnsim(j,ic)
                 count_wnc = count_wnc + 1.0_rk
             end if
        end do
        
        
        do j = 1,jnum
            if (bsim(j,ic).lt.0.0_rk) then
                totborrower = totborrower + 1_rk
            end if


            if (j.lt.12_ik) then
                debtbyagesim(1) = debtbyagesim(1) + bsim(j,ic)
            elseif (j.ge.12 .and. j.lt.22) then
                debtbyagesim(2) = debtbyagesim(2) + bsim(j,ic)
            elseif (j.ge.22 .and. j.lt.32) then
                debtbyagesim(3) = debtbyagesim(3) + bsim(j,ic)
            elseif (j.ge.32 .and. j.lt.42) then
                debtbyagesim(4) = debtbyagesim(4) + bsim(j,ic)
            else
                debtbyagesim(5) = debtbyagesim(5) + bsim(j,ic)
            end if
        end do

    end do

    do ia = 1,4
        laborpt(ia) =sum(avglaborini(1:3,ia))/sum(colrate_ini(1:3,ia))
        stdloanpt(ia) = sum(avgbsimini(1:3,ia))/sum(borrower(1:3,ia))!sum(colrate(1:3,ia))
        tuipt(ia) = sum(avgnetphi(1:3,ia))/sum(colrate_ini(1:3,ia))
        earnpt(ia) = sum(avgearnini(1:3,ia))/sum(colrate_ini(1:3,ia)) !this is only for j = 1
    end do


    do iab = 1,3
        laborabil(iab) = sum(avglaborini(iab,1:4))/sum(colrate_ini(iab,1:4))
        stdloanabil(iab) = sum(avgbsimini(iab,1:4))/sum(borrower(iab,1:4))!sum(colrate(iab,1:4))
        tuiabil(iab) = sum(avgnetphi(iab,1:4))/sum(colrate_ini(iab,1:4))
        earnabil(iab) = sum(avgearnini(iab,1:4))/sum(colrate_ini(iab,1:4)) !this is only for j = 1
    end do

    stdloanpt = (stdloanpt/yagg)*gdp_data/4
    stdloanabil =(stdloanabil/yagg)*gdp_data/4
    tuipt = (tuipt/yagg)*gdp_data
    tuiabil = (tuiabil/yagg)*gdp_data
    earnabil = (earnabil/yagg)*gdp_data
    earnpt = (earnpt/yagg)*gdp_data    
    
    avgdebt = -(avgdebt/total_col)/yagg
    avgdebt_c = -(avgdebt_c/colh)/yagg
    avgdebt_d = -(avgdebt_d/cold)/yagg    
    avgability = avgability/total_col    
    
    avgpt = ((avgpt_c+avgpt_nc)/(avgpttotal_c+avgpttotal_nc))/yagg
    avgpt_c = (avgpt_c/avgpttotal_c)/yagg
    avgpt_nc = (avgpt_nc/avgpttotal_nc)/yagg
    stddebtholder = total_bor/total_col    
    avgnettui = (avgnettui/total_col)/yagg    
    
    debtbyagesim =((debtbyagesim/yagg)/(nsim*jnum))*yaggdata   
    totborrower = (totborrower/(nsim*jnum))*totalpopdata

    avglaborini = avglaborini/colrate_ini*18*5
    avgasimini = ((avgasimini/temp3)/yagg)*gdp_data !parental transfer
    avgbsimini = ((avgbsimini/borrower)/yagg)*gdp_data/4 !student loan
    avgearnini = ((avgearnini/colrate_ini)/yagg)*gdp_data !average earning at j=1    
    avgnetphi  = ((avgnetphi/colrate_ini)/yagg)*gdp_data !tuition

    earn_w = (earn_w/count_w)/yagg 
    earn_wc = (earn_wc/count_wc)/yagg 
    earn_wd = (earn_wd/count_wd)/yagg 
    earn_wnc = (earn_wnc/count_wnc)/yagg 
    
    ticagg = sum(tic)/total_borj1
    excagg = sum(exc)/total_borj1
    const = const/total_bor
    const1 = const1/total_borj1
    
    tic_abil_wealth(1:3,1:4) = tic(1:3,1:4)/sum(tic)
    exc_abil_wealth(1:3,1:4) = exc(1:3,1:4)/sum(exc)

    do iab = 1, 3
        colrate_abil_ini(iab) = (colrate_abil_ini(iab))/sum(total_ini_d(iab, 1:4))
        tic_abil(iab) = tic_abil(iab)/sum(tic)
        exc_abil(iab) = exc_abil(iab)/sum(exc)
    end do

    do iab = 1, 4
        colrate_wealth_ini(iab) = (colrate_wealth_ini(iab))/sum(total_ini_d(1:3, iab))
        tic_wealth(iab) = tic_wealth(iab)/sum(tic)
        exc_wealth(iab) = exc_wealth(iab)/sum(exc)
    end do    
    
    abilratio =   tuiabil(3)/tuiabil(1)
    wealthratio = tuipt(4)/tuipt(1)

    colrateaggd = colrateaggd/total
    colrateagg =  educated/(educated+noneducated)

    colrate_ini_d = colrate_ini_d(1:3,1:4)/total_ini_d(1:3,1:4)

    debtpercent =  sum(ddebt(1:jnum,1:nsim))/sum(bsim)
    ddebtpercentmu = (drateh+drated)/(count1_h+count1_d)! sum(ddebtmu(jd:jrnum,1:nsim))/count1
    drateh = drateh/count1_h
    drated = drated/count1_d    
    debtout = -(sum(bsim)/(nsim*jnum))/yagg

    !aggregate labor supply and earning while in college 
    avgearn = (avgearn/total_col)/yagg    
    avglaborhini = (avglaborhini/total_col)*18*5
    
    laborpt = laborpt*18*5
    laborabil = laborabil*18*5
    

    
    
    if (iss .eq. 1_ik) then
        if (t.ne. 0_ik) then
        write(*,*)'--------------------------->>> simulation results for time: ', 1983 + t
        end if
        write(*,*)
        write(*,'(1x,3(a,f8.4),a)') 'Graduation rate (25-85 pop)                                0.1887 :', colrateagg
        write(*,'(1x,3(a,f8.4),a)') 'Graduation rate (j=1:jc-1)                                        :', colrateaggd
        write(*,'(1x,3(a,f8.4),a)') 'students graduating with student debt in simulation(1997)    0.60 :', stddebtholder
        write(*,'(1x,3(a,f8.4),a)') 'average student debt in simulation                         0.1680 :', avgdebt
        write(*,'(1x,3(a,f8.4),a)') 'outstanding std debt / gdp  in simulation                  0.0070 :', debtout
        write(*,'(1x,3(a,f8.4),a)') 'delinquent debt in simulation                                     :', debtpercent
        write(*,'(1x, a, f8.4)')    'Average net college tuition in simulation(w/ selection)    0.1809 :', avgnettui
        write(*,'(1x,3(a,f8.4),a)') 'averate productivity of college goers                             :', avgability
        !write(*,'(1x,3(a,f8.4),a)') 'fraction of constrained people(tic)                               :', ticagg
        !write(*,'(1x,3(a,f8.4),a)') 'fraction of constrained people(exc)                               :', excagg
        !write(*,'(1x,3(a,f8.4),a)') 'fraction of constrained people(max(tic,exc))                      :', const
        !write(*,'(1x,3(a,f8.4),a)') 'fraction of constrained people(max(tic,exc) at j=1)               :', const1
        write(*,'(1x,3(a,f8.4),a)') 'average earnings while in college in sim(8884/31123)        0.2854 :', avgearn        
        write(*,'(1x,3(a,f8.4),a)') 'Average labor supply weekly in hours                        23.61 :', avglaborhini   
        write(*,*)        
        write(*,'(1x, a, f8.4)')    ' Average parental transfer: $4,719.56 (0.1516)                         : ', avgpt
        write(*,'(1x, a, f8.4)')    ' Average parental transfer for college: $9,722.83 (0.3124)             : ', avgpt_c
        write(*,'(1x, a, f8.4)')    ' Average parental transfer for non-college: $3,230.14 (0.1038)         : ', avgpt_nc   
        write(*,*)
        write(*,'(1x, a, f8.4)')    ' Average debt at graduation for all              : $5,229 (0.1680)     : ', avgdebt
        write(*,'(1x, a, f8.4)')    ' Average debt at graduation for college          : $9,493 (0.3050)     : ', avgdebt_c
        write(*,'(1x, a, f8.4)')    ' Average debt at graduation for dropout          : $2,618 (0.0841)     : ', avgdebt_d    
        write(*,*)   
        write(*,'(1x, a, f8.4)')    ' delinqunecy rate for all workers (SCF1998)      : 0.1863              : ',  ddebtpercentmu
        write(*,'(1x, a, f8.4)')    ' delinqunecy rate for college graduate workers   : 0.1388              : ',  drateh
        write(*,'(1x, a, f8.4)')    ' delinqunecy rate for college dropouts workers   : 0.2523              : ',  drated    
        write(*,*)  
        write(*,'(1x, a, f8.4)')    ' earnings for all workers (SCF1998)          : $51,625.15 (1.6587)     : ',  earn_w     
        write(*,'(1x, a, f8.4)')    ' earnings for college graduate workers       : $63,756.82 (2.0485)     : ',  earn_wc
        write(*,'(1x, a, f8.4)')    ' earnings for college dropouts workers       : $51,567.80 (1.6569)     : ',  earn_wd 
        write(*,'(1x, a, f8.4)')    ' earnings for college non-college workers    : $43,681.75 (1.4035)     : ',  earn_wnc          
        write(*,*)           
        write(*,'(1x,3(a,f8.4),a)') 'college completion data (1T, 2T, 3T) = (', 0.0425, ',', 0.1524, ',', 0.3772, ')'
        write(*,'(1x,3(a,f8.4),a)') 'college completion mode (1T, 2T, 3T) = (', colrate_abil_ini(1), ',', colrate_abil_ini(2), ',', colrate_abil_ini(3), ')'
        write(*,'(1x,4(a,f8.4),a)') 'college completion data (1Q, 2Q, 3Q, 4Q) = (', 0.0970, ',', 0.1396, ',', 0.2230, ',', 0.3119,')'
        write(*,'(1x,4(a,f8.4),a)') 'college completion mode (1Q, 2Q, 3Q, 4Q) = (', colrate_wealth_ini(1), ',', colrate_wealth_ini(2), ',', colrate_wealth_ini(3), ',', colrate_wealth_ini(4),')'
        write(*,*)
        write(*,*) '---------------------------over earnings quartile 25-65--------------------------'
        write(*,'(1x,4(a,f8.4),a)') 'default rate among borrowers SCF1998  (1Q, 2Q, 3Q, 4Q) = (', 0.2680, ',', 0.2208, ',', 0.2223,  ',', 0.0967,')'
        write(*,'(1x,4(a,f8.4),a)') 'default rate among borrowers          (1Q, 2Q, 3Q, 4Q) = (', drateearn(1), ',', drateearn(2), ',', drateearn(3),  ',', drateearn(4),')'
        write(*,'(1x,4(a,f8.0),a)') 'bval among borrowers                  (1Q, 2Q, 3Q, 4Q) = (', bearnq(1), ',', bearnq(2), ',', bearnq(3),  ',', bearnq(4),')'
        write(*,'(1x,4(a,f8.4),a)') 'dropouts among college goers          (1Q, 2Q, 3Q, 4Q) = (', cdrop_earn(1), ',', cdrop_earn(2), ',', cdrop_earn(3),  ',', cdrop_earn(4),')'
        write(*,*)        
        do i= 1,4
            write(*,'(1x,3(a,f8.0),a)') 'total simulated sample (1Q, 2Q, 3Q) = (', total_ini_d(1,i), ',', total_ini_d(2,i), ',', total_ini_d(3,i), ')'
        end do
        do i=1,4
            write(*,'(1x,3(a,f8.4),a)') 'college completion rate (1Q, 2Q, 3Q) = (', colrate_ini_d(1,i), ',', colrate_ini_d(2,i), ',', colrate_ini_d(3,i), ')'
        end do
        write(*,*)
        write(*,*) 'college completion rate in data (1Q, 2Q, 3Q) = (0.0099 , 0.0500, 0.3212)'
        write(*,*) 'college completion rate in data (1Q, 2Q, 3Q) = (0.0137 , 0.0560, 0.2042)'
        write(*,*) 'college completion rate in data (1Q, 2Q, 3Q) = (0.0059 , 0.0698, 0.2517)'
        write(*,*) 'college completion rate in data (1Q, 2Q, 3Q) = (0.0130 , 0.1103, 0.4053)'
        write(*,*)
        write(*,*)
        !do i=1,4
        !    write(*,'(1x,3(a,f15.4),a)') 'initial parental transfer (1Q, 2Q, 3Q) = (', avgasimini(1,i), ',', avgasimini(2,i), ',', avgasimini(3,i), ')'
        !end do
        !write(*,*)
        !do i=1,4
        !    write(*,'(1x,3(a,f15.4),a)') '    average studetn debt (1Q, 2Q, 3Q) = (', avgbsimini(1,i), ',', avgbsimini(2,i), ',', avgbsimini(3,i), ')'
        !end do
        !write(*,*)
        !do i=1,4
        !    write(*,'(1x,3(a,f15.4),a)') '     average net tuition (1Q, 2Q, 3Q) = (', avgnetphi(1,i), ',', avgnetphi(2,i), ',', avgnetphi(3,i), ')'
        !end do
        !write(*,*)
        !do i=1,4
        !    write(*,'(1x,3(a,f15.4),a)') '    average labor supply (1Q, 2Q, 3Q) = (', avglaborini(1,i), ',', avglaborini(2,i), ',', avglaborini(3,i), ')'
        !end do
        !write(*,*)
        !do i=1,4
        !    write(*,'(1x,3(a,f15.4),a)') '    average earnings (1Q, 2Q, 3Q) = (', avgearnini(1,i), ',', avgearnini(2,i), ',', avgearnini(3,i), ')'
        !end do
        write(*,*)
        write(*,*) ' labor supply over ability tercile'
        write(*,*) '    data          (1T, 2T, 3T)   = (   17.11,   17.25,   18.22)'
        write(*,'(1x,3(a,f8.2),a)') 'average labor supply (1T, 2T, 3T)     = (', laborabil(1), ',', laborabil(2), ',', laborabil(3), ')'
        write(*,*) '    data          (1Q, 2Q, 3Q, 4Q) = (   16.08,   18.48,   18.17,   17.96)'
        write(*,'(1x,4(a,f8.2),a)') 'average labor supply (1Q, 2Q, 3Q, 4Q) = (', laborpt(1), ',', laborpt(2), ',', laborpt(3),  ',', laborpt(4),')'
        write(*,*)
        write(*,*)
        write(*,*) '    data          (1T, 2T, 3T)   = (  345,   833,   972)'
        write(*,'(1x,3(a,f8.0),a)') 'annual avg student debt (1T, 2T, 3T)     = (', -stdloanabil(1), ',', -stdloanabil(2), ',', -stdloanabil(3), ')'
        write(*,*) '    data          (1Q, 2Q, 3Q, 4Q) = (  558,   662,   809,   1047)'
        write(*,'(1x,4(a,f8.0),a)') 'annual avg student debt (1Q, 2Q, 3Q, 4Q) = (', -stdloanpt(1), ',', -stdloanpt(2), ',', -stdloanpt(3),  ',', -stdloanpt(4),')'
        write(*,*)
        write(*,*)
        write(*,*) '   data(1979)     (1T, 2T, 3T)   = (  2802,   4408,   5684)'
        write(*,'(1x,3(a,f8.0),a)') 'average net tuition (1T, 2T, 3T)     = (', tuiabil(1), ',', tuiabil(2), ',', tuiabil(3), ')'
        write(*,*) '   data(1979)     (1Q, 2Q, 3Q, 4Q) = (  5223,   4152,   3752,   5984)'
        write(*,'(1x,4(a,f8.0),a)') 'average net tuition (1Q, 2Q, 3Q, 4Q) = (', tuipt(1), ',', tuipt(2), ',', tuipt(3),  ',', tuipt(4),')'
        write(*,*)
        write(*,*)
        write(*,*) '   data(1979)     (1T, 2T, 3T)   = (  5640,   7100,   8058)'
        write(*,'(1x,3(a,f8.0),a)') 'avg earning during college (1T, 2T, 3T)     = (', earnabil(1), ',', earnabil(2), ',', earnabil(3), ')'
        write(*,'(1x,4(a,f8.0),a)') 'avg earning during college (1Q, 2Q, 3Q, 4Q) = (', earnpt(1), ',', earnpt(2), ',', earnpt(3),  ',', earnpt(4),')'

        write(*,'(1x, a, f8.4)')    '  A7.1 Tuition ability-ratio      :    9318/4593 = 2.0287   ',   abilratio
        write(*,'(1x, a, f8.4)')    '  A7.1 Tuition PT-ratio           :    9810/8562 = 1.1458   ',   wealthratio

        !write(*,*)
        !write(*,*)
        !write(*,'(1x, a, f12.2)')    'debt_under30 / tot_debt age in simulation (data: 0.36)', (debtbyagesim(1))/sum(debtbyagesim)
        !write(*,'(1x, a, f12.2)')    'debt_30-39   / tot_debt age in simulation (data: 0.33)', (debtbyagesim(2))/sum(debtbyagesim)
        !write(*,'(1x, a, f12.2)')    'debt_40-49   / tot_debt age in simulation (data: 0.16)', (debtbyagesim(3))/sum(debtbyagesim)
        !write(*,'(1x, a, f12.2)')    'debt_50-59   / tot_debt age in simulation (data: 0.11)', (debtbyagesim(4))/sum(debtbyagesim)
        !write(*,'(1x, a, f12.2)')    'debt_60+     / tot_debt age in simulation (data: 0.04)', (debtbyagesim(5))/sum(debtbyagesim)


    end if

    if (t.eq.0_ik) then

        open(unit = 14, file = 'debtbyagesim_ss.txt', status = 'replace', action = 'write')
        do i=1,5
            write(14,*),(debtbyagesim(i))/sum(debtbyagesim(:))
        end do
        close(14)

        open(unit = 22, file = 'colratesim_ss.txt', status = 'replace', action = 'write')
        do i=1,3
            do j=1,4
                write(22,*),colrate_ini(i,j)
            end do
        end do
        close(22)
    end if


    if (itext.eq.1_ik .and. t .eq. 0_ik) then
        open(unit = 77, file = 'ss_sim_result.txt', status = 'replace', action = 'write')
        write(77,'(1x,3(a,f8.0),a)') 'Average labor supply weekly in hours                        18.29 :', avglaborhini
        write(77,'(1x,3(a,f8.0),a)') 'average earnings while in college in simulation            0.2971 :', avgearn*gdp_data
        write(77,'(1x, a, f8.0)')    'Average net college tuition in simulation(w/ selection)    0.1909 :', avgnettui*gdp_data
        write(77,'(1x, a, f8.0)')    ' Average parental transfer: $4.719.56                      0.1907 :', avgpt*gdp_data
        write(77,'(1x,3(a,f8.4),a)') 'outstanding std debt / gdp  in simulation                  0.0092 :', debtout
        write(77,'(1x,3(a,f8.4),a)') 'delinquent population in simulation                       17%-21% :', ddebtpercentmu
        write(77,'(1x,3(a,f8.4),a)') 'college completion mode 1T                                 0.0103 :', colrate_abil_ini(1)
        write(77,'(1x,3(a,f8.4),a)') 'college completion mode 2T                                 0.0810 :', colrate_abil_ini(2)
        write(77,'(1x,3(a,f8.4),a)') 'college completion mode 3T                                 0.3315 :', colrate_abil_ini(3)
        write(77,'(1x,4(a,f8.4),a)') 'college completion mode 1Q                                 0.0903 :', colrate_wealth_ini(1)
        write(77,'(1x,4(a,f8.4),a)') 'college completion mode 2Q                                 0.0783 :', colrate_wealth_ini(2)
        write(77,'(1x,4(a,f8.4),a)') 'college completion mode 3Q                                 0.1165 :', colrate_wealth_ini(3)
        write(77,'(1x,4(a,f8.4),a)') 'college completion mode 4Q                                 0.2466 :', colrate_wealth_ini(4)
        write(77,'(1x,3(a,f8.4),a)') 'Graduation rate in 1979 in simulation                      0.1563 :', colrateagg
        close(77)
    end if
    
        if (t .eq. 19_ik) then !saving results for year 1997 from transition
            open(unit = 19, file = 'laborsim_97.txt', status = 'replace', action = 'write')
            do i=1,4
                write(19,*),laborpt(i)
            end do
            do i=1,3
                write(19,*),laborabil(i)
            end do
            close(19)
            open(unit = 20, file = 'stdloansim_97.txt', status = 'replace', action = 'write')
            do i=1,4
                write(20,*),stdloanpt(i)
            end do
            do i=1,3
                write(20,*),stdloanabil(i)
            end do
            close(20)
            open(unit = 21, file = 'tuisim_97.txt', status = 'replace', action = 'write')
            do i=1,4
                write(21,*),tuipt(i)
            end do
            do i=1,3
                write(21,*),tuiabil(t)
            end do
            close(21)
            open(unit = 23, file = 'colratesim_97.txt', status = 'replace', action = 'write')
            do i=1,3
                do j=1,4
                    write(23,*),colrate_ini_d(i,j)
                end do
            end do
            close(23)
            open(unit = 24, file = 'debtbyagesim_97.txt', status = 'replace', action = 'write')
            do i=1,5
                write(24,*),(debtbyagesim(i))/sum(debtbyagesim(:))
            end do
            close(24)
            open(unit = 25, file = 'defaultearns_97.txt', status = 'replace', action = 'write')
            write(*,'(1x,4(a,f8.4),a)') 'default rate  (1Q, 2Q, 3Q, 4Q) = (', drateearn(1), ',', drateearn(2), ',', drateearn(3),  ',', drateearn(4),')'
            close(25)
        end if

    end subroutine sim_moments


    subroutine prodsimulate(e,  jnum1, num, nsim, grid, piz1, piz2, piz10, piz20, isimh, simh, isiml, siml)

    integer(ik)::num, isimh(jnum1, nsim), isiml(jnum1, nsim),nsim, i, minplace(1), j, jnum1, ind0,  e
    !integer(ik), allocatable:: seed(:)
    real(rk):: grid(num), simh(jnum1, nsim), siml(jnum1, nsim), choicevec(jnum1,  nsim), cumpiz(num, num), &
        piz1(num, num), piz2(num, num), cumpiz0(num), piz10(num), piz20(num),cumpizi0(num), temp(num)

    intent(in):: e,  num, nsim, grid, piz1, piz2,  piz10, piz20, jnum1
    intent(out):: isimh, simh, isiml, siml

    !!fix the randome number seed in the beginning
    !!allocate(seed(nsim))
    !!seed(:) = 86934 !fix the seed
    !seed(:) = 432
    !!random number generator
    !call random_seed(put=seed)
    
    !call set_seed(432)
    !call random_seed(size=n)
    !allocate(seed(n))
    !seed = 432    ! putting arbitrary seed to all elements
    !call random_seed(put=seed) 
    call random_number(choicevec)

    
    if (e.eq.1_ik) then
        ! to draw the productivity ie =1 for the college education period
        isimh(1:jc-1,1:nsim) = 1_ik
        do j = 1, jc-1
            do i = 1, nsim
                simh(j,i) = grid(isimh(j,i))
            end do
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! simulate for the college-educated!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !random number generator
        !call random_seed(put=seed)
        call random_number(choicevec)

        !calculate cumulative distribution
        cumpiz0 = 0.0_rk
        cumpiz0(1) = piz10(1) !initial distribution
        do j = 2,num
            cumpiz0(j) = cumpiz0(j-1) + piz10(j)
        end do
        call cumsum(piz1, num,cumpiz) !for the markov chain

        !simulate data
        do i = 1,nsim
            cumpizi0 =0.0_rk
            cumpizi0(1:num) = cumpiz0(1:num)-choicevec(1, i)
            minplace = minloc(cumpizi0, mask = cumpizi0 .gt. 0.0_rk)
            isimh(jc,i) = minplace(1)
            simh(jc,i) = grid(isimh(jc,i))
        end do


        do i = 1, nsim
            do j =jc+1,jnum1
                if ( j .le. jrnum+1) then
                    temp =0.0_rk
                    ind0 = isimh(j-1,i)
                    temp(1:num) = cumpiz(ind0,1:num)-choicevec(j,i)
                    minplace = minloc(temp, mask = temp .gt. 0.0_rk)
                    isimh(j,i) = minplace(1)
                    simh(j,i) = grid(isimh(j,i))
                else
                    isimh(j,i) = isimh(j-1,i)
                    simh(j,i) = simh(j-1,i)
                end if
            end do
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !simulate for non-college-educated!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !random number generator
        !call random_seed(put=seed)
        call random_number(choicevec)

        !calculate cumulative distribution
        cumpiz0=0.0_rk
        cumpiz0(1) = piz20(1) !initial distribution
        do j = 2,num
            cumpiz0(j) = cumpiz0(j-1) + piz20(j)
        end do
        call cumsum(piz2, num,cumpiz) !for the markov chain

        !simulate data
        do i = 1,nsim
            cumpizi0 =0.0_rk
            cumpizi0(1:num) = cumpiz0(1:num)-choicevec(1, i)
            minplace = minloc(cumpizi0, mask = cumpizi0 .gt. 0.0_rk)
            isiml(1,i) = minplace(1)
            siml(1,i) = grid(isiml(1,i))
        end do

        !calculate cumulative distribution

        do j =2,jnum1
            do i =1,nsim
                if (j .le. jrnum+1 ) then
                    temp =0.0_rk
                    ind0 = isiml(j-1,i)
                    temp(1:num) = cumpiz(ind0,1:num)-choicevec(j,i)
                    minplace = minloc(temp, mask = temp .gt. 0.0_rk)
                    isiml(j,i) = minplace(1)
                    siml(j,i) = grid(isiml(j,i))
                else
                    isiml(j,i) = isiml(j-1,i)
                    siml(j,i) = siml(j-1,i)
                end if
            end do
        end do

    else

        !college- dropouts: here j=1 corresponds to jc

        call random_number(choicevec)

        ! simulate for the college-educated
        !calculate cumulative distribution
        cumpiz0=0.0_rk
        cumpiz0(1) = piz20(1) !initial distribution
        do j = 2,num
            cumpiz0(j) = cumpiz0(j-1) + piz20(j)
        end do
        call cumsum(piz2, num,cumpiz) !for the markov chain

        !simulate data
        do i = 1,nsim
            cumpizi0 =0.0_rk
            cumpizi0(1:num) = cumpiz0(1:num)-choicevec(1, i)
            minplace = minloc(cumpizi0, mask = cumpizi0 .gt. 0.0_rk)
            isiml(1,i) = minplace(1)
            siml(1, i) = grid(isiml(1,i))
        end do


        do i = 1, nsim
            do j =2,jnum1
                if (j .le. jrnum-jc+1) then
                    temp =0.0_rk
                    ind0 = isiml(j-1,i)
                    temp(1:num) = cumpiz(ind0,1:num)-choicevec(j,i)
                    minplace = minloc(temp, mask = temp .gt. 0.0_rk)
                    isiml(j,i) = minplace(1)
                    siml(j,i) = grid(isiml(j,i))
                else
                    isiml(j,i) = isiml(j-1,i)
                    siml(j,i) = siml(j-1,i)
                end if

            end do
        end do
    end if

    end subroutine prodsimulate


    subroutine iniwealthsim(num1, nsim, grid, pia0, idxsim, sim)
    integer(ik)::  j, num1, nsim, idxsim(nsim), minplace(1),i
    !integer(ik),allocatable:: seed(:)
    real(rk):: choicevec(nsim), cumpia0(num1), pia0(num1), cumpiai0(num1), sim(nsim), grid(num1)
    intent(in):: num1, nsim, grid, pia0
    intent(out):: idxsim, sim

    !allocate(seed(nsim))
    !seed(:) = 86934 !fix the seed
   ! seed(:) = 432
    !random number generator
    !call random_seed(put=seed)
    !call set_seed(432)
    !call random_seed(size=n)
    !allocate(seed(n))
    !seed = 123456789    ! putting arbitrary seed to all elements
    !call random_seed(put=seed) 
    !call set_seed(123)
    call random_number(choicevec)

    
    
    !calculate cumulative distribution
    cumpia0=0.0_rk
    cumpia0(1) = pia0(1) !initial distribution
    do j = 2,num1
        cumpia0(j) = cumpia0(j-1) + pia0(j)
    end do
    do i = 1,nsim
        cumpiai0 =0.0_rk
        cumpiai0(1:num1) = cumpia0(1:num1)-choicevec(i)
        minplace = minloc(cumpiai0, mask = cumpiai0 .gt. 0.0_rk)
        idxsim(i) = minplace(1)
        sim(i) = grid(idxsim(i))
    end do

    end subroutine iniwealthsim

   subroutine set_seed(index)
        ! Determines the dimension of the seed vector and sets the seed of
        ! the random number generator in a deterministic way.
        ! The function must be called before calling RANDOM_NUMBER
        ! Usage:
        ! call set_seed(index)
        ! Inputs:
        ! index  - choosen seed
        integer, intent(in) :: index
        integer             :: i, n
        integer, dimension(:), allocatable :: seed
        call random_seed(size = n)
        allocate(seed(n))
        seed = index + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = seed)
    end subroutine set_seed


    end module StdDebtCrisis_mod