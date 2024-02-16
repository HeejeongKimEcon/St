    module transition_solve

    use Lib_kindset
    use Lib_conshock
    use Lib_inequality
    use Lib_sort
    use Lib_grid
    use omp_lib
    use lib_rwhdf5
    use StdDebtCrisis_mod
    use get_params
    implicit none

    contains


    subroutine transtioninput(nabil, pabil, anum,  wh, netphi, probcost, psycost,   piec, piel, epgridc, etgridc, eagridj, &
        psycost_t, probcost_t, wh_t, phi_t, bbar_t, piec_t, piel_t, piel0_t, piec0_t, pia0, pic0, piea0_t, eagridj_t)

    integer(ik):: anum, nabil
    integer(ik):: ie, it, ic, ia
    real(rk):: wh, netphi(nabil, anum),  probcost(nabil,anum,ncost), psycost(nabil,anum,ncost)
    real(rk):: piec(e_num, e_num), piel(e_num, e_num), epgridc(epnum), etgridc(etnum), &
        varepc_t(tnum), varepnc_t(tnum), varetc_t(tnum), varetnc_t(tnum), pia0(anum), pic0(ncost), &
        pabil(nabil), discount, temp, eagridj(jnum, nabil, e_num, 2)

    real(rk), allocatable:: piec_t(:,:,:), piel_t(:,:,:), piel0_t(:,:), piec0_t(:,:),phi_t(:,:,:),  wh_t(:), tuition_t(:), bbar_t(:), &
        probcost_t(:,:,:,:), psycost_t(:,:,:,:), piea0_t(:,:,:), phiavg_t(:), grah_t(:), gral_t(:), eagridj_t(:,:,:,:,:)

    intent(in)::  anum,  wh, netphi, probcost, psycost,  piec, piel, epgridc, etgridc,  pia0, pic0, nabil,  pabil, eagridj
    intent(out):: psycost_t, probcost_t, wh_t, phi_t, bbar_t, piec_t, piel_t, piel0_t, piec0_t, piea0_t, eagridj_t

    allocate( piec_t(e_num,e_num,tnum), piel_t(e_num,e_num,tnum), piel0_t(e_num,tnum), piec0_t(e_num,tnum),phi_t(nabil, anum,tnum),  wh_t(tnum), &
        eagridj_t(jnum, nabil, e_num, 2, tnum))

    !Calling wage premium
    if (iwp .eq. 1_ik) then
        wh_t(1:tnum) = wh
    else
        OPEN(3, FILE='wagepremium.txt', &
            STATUS='OLD', ACTION='READ')
        do ie = 1, tnum
            if (ie .le. tmid) then
                read(3,*), wh_t(ie)
            else
                wh_t(ie) = wh_t(tmid)
            end if
        end do
        close(3)
    end if


    !tuition (1997-2015 , using until 2015)
    allocate(tuition_t(tnum), phiavg_t(tnum))
    if (ituition .eq. 1_ik) then
        tuition_t(1:tnum) = 1.0_rk
    else
        OPEN(3, FILE='phi_t.txt', STATUS='OLD', ACTION='READ')
        do ie = 1, tnum
            if(ie.le.tmid) then
                read(3,*), tuition_t(ie)
            else
                tuition_t(ie) = tuition_t(tmid)
            end if
        end do
        close(3)
    end if


    !ability premium(1997-2015 , using until 2015)
    allocate(grah_t(tnum), gral_t(tnum))
    OPEN(3, FILE='skilled_abil_HP_time.txt', STATUS='OLD', ACTION='READ')
    do ie = 1, tnum
        if (ie.le.tmid) then
            read(3,*), grah_t(ie)
        else
            grah_t(ie) = grah_t(tmid)
        end if
    end do
    close(3)
    OPEN(3, FILE='unskilled_abil_HP_time.txt', STATUS='OLD', ACTION='READ')
    do ie = 1, tnum
        if (ie.le.tmid) then
            read(3,*), gral_t(ie)
        else
            gral_t(ie) = gral_t(tmid)
        end if
    end do
    close(3)

    if (iability .eq. 1_ik) then
        grah_t(2:tnum) = grah_t(1)
        gral_t(2:tnum) = gral_t(1)
    end if

    ! earnings grid updated with time-varying ability premium
    do it = 1, tnum
        eagridj_t(:,:,:,1, it) = eagridj(:,:,:,1)*(exp(grah_t(it))/exp(grah_t(1)))
        eagridj_t(:,:,:,2, it) = eagridj(:,:,:,2)*(exp(gral_t(it))/exp(gral_t(1)))
    end do

    !phiavg = 0.0_rk
    !do iab = 1,nabil
    !    do ia = 1,anum
    !        phiavg= phiavg+ netphi(iab,ia)*pabil(iab)*pia0(ia)
    !    end do
    !end do
    !
    !
    !phiavg_t = 0.0_rk
    !do ie=1,tnum
    !
    !    do iab = 1,nabil
    !        do ia = 1,anum
    !            phi_t(iab,ia,ie) = netphi(iab,ia) + phiavg*(tuition_t(ie)-1.0_rk)
    !            phiavg_t(ie) = phiavg_t(ie) + phi_t(iab,ia,ie)*pabil(iab)*pia0(ia)
    !        end do
    !    end do
    !    write(*,*) ie, phiavg_t(ie)/phiavg
    !end do

    do ie=1,tnum
        phi_t(1:nabil,1:anum,ie) = netphi(1:nabil, 1:anum)*tuition_t(ie)
    end do


    !credit limit (1997-2015)
    allocate(bbar_t(tnum))
    bbar_t = 0.0_rk;
    bbar_t(1:tnum) = bbar

    !!!!!!Earning shocks (1993)
    call transitionshock( epgridc, etgridc, piec,piel,piec_t, piel_t, piel0_t, piec0_t,varepc_t,varepnc_t,varetc_t,varetnc_t)

    !decreasing psychic education cost
    allocate(probcost_t(nabil,anum,ncost,tnum), psycost_t(nabil,anum,ncost,tnum))
    do it = 1,tnum
        psycost_t(1:nabil, 1:anum, 1:ncost, it) =psycost(1:nabil, 1:anum, 1:ncost)
        probcost_t(1:nabil, 1:anum, 1:ncost, it) =probcost(1:nabil, 1:anum, 1:ncost)
    end do

    discount = 80_rk!45_rk

    if (ipsycost .eq. 0_ik) then
        do it = 1,tnum
            temp = 0.0_rk
            if (it.gt.1_ik .and. it .le. tmid) then
                !psycost_t(1:nabil, 1:anum, 1:ncost, it) = ((1-discount)**(it-1_ik))*psycost_t(1:nabil, 1:anum, 1:ncost, it)
                psycost_t(1:nabil, 1:anum, 1:ncost, it) = psycost_t(1:nabil, 1:anum, 1:ncost, it) - (discount/tnum)*it
                !write(*,*) 'decreasing psycost in ', it-1_ik+1979, ' by: ', ((1-discount)**(it-1_ik))
                write(*,*) 'decreasing psycost in ', it-1_ik+1983, ' by: ', (discount/tnum)*it
            elseif (it .gt. tmid) then
                psycost_t(1:nabil, 1:anum, 1:ncost, it) = psycost_t(1:nabil, 1:anum, 1:ncost, tmid)
            end if

            do ie = 1, nabil
                do ia = 1, anum
                    do ic = 1, ncost
                        probcost_t(ie,ia,ic,it) =pic0(ic)*pabil(ie)*pia0(ia)
                        temp = temp + probcost_t(ie,ia,ic,it)*psycost_t(ie,ia,ic,it)
                    end do
                end do
            end do
            write(*,*) 'avg psychost', temp
        end do
    end if

    allocate(piea0_t(nabil,anum,tnum))
    !Time varying piea0_t
    do it = 1, tnum
        do ie= 1, nabil
            do ia = 1,anum
                piea0_t(ie,ia,it) = pia0(ia)*pabil(ie)
            end do
        end do
    end do

    end subroutine transtioninput


    subroutine transitionshock( epgridc, etgridc, piec,piel, piec_t, piel_t, piel0_t, piec0_t, varepc_t,varepnc_t,varetc_t,varetnc_t)

    integer(ik)::  ie, place,ip,it,ipf,itf,placef,t
    real(rk):: varepc_t(tnum),varepnc_t(tnum), varetc_t(tnum),varetnc_t(tnum), &
        sigmacp_t(tnum),sigmact_t(tnum),sigmancp_t(tnum),sigmanct_t(tnum),piec_t(e_num,e_num,tnum), &
        piel_t(e_num,e_num,tnum), piepc(epnum,epnum), pietc(etnum,etnum), &
        piepl(epnum,epnum), pietl(etnum,etnum), piec(e_num,e_num), piel(e_num,e_num), epgridc(epnum), etgridc(etnum),&
        piel0_t(e_num,tnum), piec0_t(e_num,tnum)


    intent(in)::  epgridc, etgridc, piec,piel
    intent(out)::piec_t, piel_t, piel0_t, piec0_t,varepc_t,varepnc_t,varetc_t,varetnc_t


    if ( ishock .eq. 1_ik) then
        varepc_t(1:tnum) = varp1;
        varepnc_t(1:tnum) = varp2;
        varetc_t(1:tnum) = vart1;
        varetnc_t(1:tnum) = vart2;
        do t = 1, tnum
            piec_t(1:e_num,1:e_num,t) = piec(1:e_num,1:e_num)
            piel_t(1:e_num,1:e_num,t) = piel(1:e_num,1:e_num)
        end do
        do t = 1,tnum
            call  ergodicdist(e_num, piel_t(1:e_num,1:e_num,t), precerg, piel0_t(1:e_num,t))
            call  ergodicdist(e_num, piec_t(1:e_num,1:e_num,t), precerg, piec0_t(1:e_num,t))
        end do
    else

        OPEN(30, FILE='varhpcp_t.txt', &
            STATUS='OLD', ACTION='READ')
        do ie = 1, tnum
            if (ie .le. tmid) then
                read(30,*), varepc_t(ie)
            else
                varepc_t(ie) = varepc_t(tmid)
            end if
        end do
        close(30)

        OPEN(3, FILE='varhpcp_t.txt', &
            STATUS='OLD', ACTION='READ')
        do ie = 1, tnum
            if (ie .le. tmid) then
                read(3,*), varepnc_t(ie)
            else
                varepnc_t(ie) = varepnc_t(tmid)
            end if
        end do
        close(3)

        OPEN(5, FILE='varhpct_t.txt', &
            STATUS='OLD', ACTION='READ')
        do ie = 1, tnum
            if (ie .le. tmid) then
                read(5,*), varetc_t(ie)
            else
                varetc_t(ie) = varetc_t(tmid)
            end if
        end do
        close(5)

        OPEN(10, FILE='varhpct_t.txt', &
            STATUS='OLD', ACTION='READ')
        do ie = 1, tnum
            if (ie .le. tmid) then
                read(10,*), varetnc_t(ie)
            else
                varetnc_t(ie) = varetnc_t(tmid)
            end if
        end do
        close(10)

        !Labor productivity shock process for college-educated & non-college-educated households
        !1 is for college educated, 2 is for non-college educated


        do ie = 1,tnum
            sigmacp_t(ie) = dsqrt(varepc_t(ie))
            sigmancp_t(ie) = dsqrt(varepnc_t(ie))
            sigmact_t(ie) = dsqrt(varetc_t(ie))
            sigmanct_t(ie) = dsqrt(varetnc_t(ie))
        end do

        !!!!!Earning shocks (1984)
        piec_t(1:e_num,1:e_num,1) = piec
        piel_t(1:e_num,1:e_num,1) = piel

        !Discretize labor shock process for college-educated & non-college-educated households
        do t = 2, tnum
            !Combining persistent and transitional shock
            call tauchenonsupport(mean,sigmacp_t(t),rho1,epnum,epgridc,piepc)
            call tauchenonsupport(mean,sigmancp_t(t),rho2,epnum,epgridc,piepl)  !epgridl
            call tauchenonsupport(mean,sigmact_t(t),0.0_rk,etnum,etgridc,pietc)
            call tauchenonsupport(mean,sigmanct_t(t),0.0_rk,etnum,etgridc,pietl) !epgridlc

            do ip = 1, epnum
                do it = 1, etnum
                    place = (ip - 1_ik)*etnum + it
                    do ipf = 1, epnum
                        do itf = 1, etnum
                            placef = (ipf - 1_ik)*etnum + itf
                            piec_t(place,placef,t) = piepc(ip,ipf)*pietc(it,itf)
                            piel_t(place,placef,t) = piepl(ip,ipf)*pietl(it,itf)
                        end do
                    end do
                end do
            end do
            do ie = 1, e_num
                piec_t(ie, 1:e_num, t) = piec_t(ie, 1:e_num, t)/sum(piec_t(ie, 1:e_num, t))
                piel_t(ie, 1:e_num, t) = piel_t(ie, 1:e_num, t)/sum(piel_t(ie, 1:e_num, t))
            end do
        end do



        !Ergodic distribution for income shock and parental transfer
        do t = 1,tnum
            call  ergodicdist(e_num, piel_t(1:e_num,1:e_num,t), precerg, piel0_t(1:e_num,t))
            call  ergodicdist(e_num, piec_t(1:e_num,1:e_num,t), precerg, piec0_t(1:e_num,t))
        end do

    end if

    end subroutine transitionshock

    subroutine transition(tbar, chibprob, eagridj_t,abilgrid, pabil, phid_abil, chibshock, nabil, anum, bnum, adnum, bdnum, wp, &
        pop,   chir, datafile, olgo, r_b, chidage, pi,  psycost_t, probcost_t, wh_t, phi_t, bbar_t, piec_t, piel_t, piel0_t, piec0_t, &
        asim, bsim, idxprod, edusim, prodsim, paysim,  idxability, abilsim,  lambda, pia0, a0grid, agrid, bgrid, adgrid, bdgrid, egridl, egridc,  &
        muhcollege_ss, muh_ss, mul_ss, mud_ss, piea0_t,nsim, idxchib)

    integer(ik) :: anum, bnum, adnum, bdnum,  nsim, nabil
    integer(ik) :: ie, t, olgo, ia,  idxprod(jnum, nsim), paysim(jnum, nsim), idxability(nsim), idxchib(nsim), t1

    real(rk) :: wp,  r_b, wl,    pi, navg, tbar

    real(rk):: pop(jnum),  lambda(nT),  pia0(anum), a0grid(nabil, anum), agrid(anum), bgrid(bnum), &
        adgrid(adnum), bdgrid(bdnum), egridc(e_num), egridl(e_num), piea0_t(nabil, anum,tnum), chibshock(chibnum), chibprob(chibnum)

    real(rk):: muh_ss(jnum,nT,nabil, e_num,adnum,bdnum), mul_ss(jnum,nT,nabil,e_num,adnum,bdnum), mud_ss(jnum,nT,nabil,e_num,adnum,bdnum), &
        parentaltransferc, parentaltransfernc, parentaltransfer, debtbyage(5), ddebt(jnum), chir, fracborrower_t(5,tsim),&
        totaldebtbyage_t(5,tsim), debtbyage_t(5,tsim), ddebtbyage_t(5,tsim), ddebtmubyage_t(5,tsim), welvalue(tnum+jnum), &
        psycost_t(nabil,anum,ncost,tnum), chidage(jnum), ddebtpercentmu_t(tsim), probcost_t(nabil,anum,ncost,tnum), &
        piec_t(e_num,e_num,tnum), piel_t(e_num,e_num,tnum), piel0_t(e_num,tnum), piec0_t(e_num,tnum),phi_t(nabil, anum,tnum), wh_t(tnum), bbar_t(tnum), &
        asim(jnum, nsim), bsim(jnum, nsim), edusim(jnum, nsim), prodsim(jnum, nsim),  &
        hitblow_t(tsim), hitahigh_t(tsim), eagridj_t(jnum, nabil, e_num, 2, tnum), pabil(nabil), abilsim(jnum,nsim),abilgrid(nabil), &
        muhcollege_ss(jc-1, nabil, adnum, bdnum, chibnum), phid_abil(nabil)


    real(rk), allocatable:: aweight(:), ainiweight(:,:), bweight(:), kagg_t(:), lagg_t(:), educated_t(:), noneducated_t(:),debtoutstanding_t(:), &
        nopayweight(:), & stddebtholder_t(:), avgnettutition(:), yagg_t(:), totalstudentdebt_t(:), avgdebt_t(:), avgdebtc_t(:), avgdebtd_t(:), &
        totalborrower_t(:),  ddebtpercent_t(:), totdebt_yagg_t(:),avgdebt_yagg_t(:),  avgtuition_yagg_t(:),&
        afsim(:,:), bfsim(:,:), edufsim(:),lcbsim(:), lcasim(:), edulclaborsim(:,:), edulcbsim(:,:), &
        edulcasim(:,:), defaultsim(:,:),  debtbyagesim(:,:), avgpt(:), avgability(:), &
        lcearnsim(:), edulcearnsim(:,:), colrate_t(:,:,:), colrate_abil_t(:,:), colrate_wealth_t(:,:),&
        abilratio_t(:), wealthratio_t(:), avgdebt_sim(:), stddebtholder_sim(:), totborrower_sim(:), totdebt_sim(:),colrate_aggsim1_t(:),&
        colrate_aggsim_t(:),ddebtratio_t(:), ddebtmuratio_t(:), ticagg_t(:),excagg_t(:), avgnettui_t(:), const_t(:),const1_t(:),        &
        bhini_t(:,:,:), vhini_t(:,:,:), colrated_t(:), drate_t(:), drateh_t(:), drated_t(:)
        

    real(rk), dimension(:,:), allocatable::  payweight, ad0weight
    real(rk), dimension(:,:,:,:,:), allocatable:: vhcollegef_t, vhcollege_t, bfchoice_t, edudecision_t,muhcollege_t
    real(rk), dimension(:,:,:,:,:,:), allocatable:: vhf_t, vh_t, vlf_t, vl_t, muh_t, mul_t, mud_t, laborhcollege_t, ghcollege_t
    real(rk), dimension(:,:,:,:,:,:,:), allocatable:: laborh_t, laborl_t, gh_t, gl_t

    integer(ik), dimension(:), allocatable:: bindex, aindex, nopayindex
    integer(ik), dimension(:,:), allocatable:: payindex, ainiindex,  ad0index
    integer(ik), dimension(:,:,:,:,:,:,:), allocatable:: pchoicec_t,pchoicenc_t


    integer(hid_t):: fileid
    character(30):: datafile,  datafile1

    intent(in):: tbar, chibprob, nabil, pabil, chibshock, anum, bnum, adnum, bdnum, wp, pop,  chir,    r_b,  chidage, pi,  lambda, &
        psycost_t, probcost_t, wh_t, phi_t, bbar_t, piec_t, piel_t, piel0_t, piec0_t, phid_abil, &
        pia0, a0grid, agrid, bgrid, adgrid, bdgrid, egridc, egridl, datafile, olgo, piea0_t, nsim, asim, bsim, idxprod, &
        edusim, prodsim, paysim,  eagridj_t, idxability, abilsim, abilgrid, idxchib, muh_ss, mul_ss, mud_ss, muhcollege_ss


    wl = 1.0_rk
    datafile1 = datafile(1:olgo+4)
    datafile1(olgo+5:olgo+12) = '_tran.h5'

    write(*,*) ' Transition start'
    write(*,*)
    write(*,*) ' Transition dynamic inputs:'

    if (ituition.eq.1_ik) then
        write(*,*) '   - tuition is fixed'
    else
        write(*,*) '   - tuition is time-varying'
    end if

    if (iwp.eq.1_ik) then
        write(*,*) '   - wage premium is fixed'
    else
        write(*,*) '   - wage premium is time-varying'
    end if

    if (ishock.eq.1_ik) then
        write(*,*) '   - variance shocks are fixed'
    else
        write(*,*) '   - variance shocks are time-varying'
    end if

    if (iability.eq.1_ik) then
        write(*,*) '   - ability premia are fixed'
    else
        write(*,*) '   - ability premia are time-varying'
    end if

    if (ipsycost.eq.1_ik) then
        write(*,*) '   - psychic costs are fixed'
    else
        write(*,*) '   - psychic costs are time-varying'
    end if

    write(*,*)
    write(*,*)
    write(*,*) ' Transition starting year is ', -(tmid-1-2015)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Transition Dynamics !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    allocate( vhini_t(nabil,anum,chibnum),  pchoicec_t(jpaynum,nT,nabil,e_num,anum,bnum,tsim), pchoicenc_t(jpaynum,nT,nabil,e_num,anum,bnum,tsim), &
        ad0index(nabil, anum), ad0weight(nabil, anum), bhini_t(nabil,anum,chibnum), edudecision_t(nabil,anum,ncost,chibnum, tsim),bfchoice_t(nabil,anum, ncost,chibnum,tsim), &
        aindex(adnum), ainiindex(nabil, adnum), aweight(adnum), ainiweight(nabil, adnum), bindex(bdnum), bweight(bdnum), payindex(bnum,nT), &
        payweight(bnum,nT), nopayindex(bnum), nopayweight(bnum))


    allocate(vhf_t(jnum,nT,nabil,e_num,anum,bnum), vlf_t(jnum,nT,nabil,e_num,anum,bnum), vhcollegef_t(jc-1_ik,nabil,anum,bnum,chibnum), &
        vl_t(jnum,nT,nabil,e_num,anum,bnum),  vh_t(jnum,nT,nabil,e_num,anum,bnum), vhcollege_t(jc-1_ik,nabil,anum,bnum,chibnum), &
        laborhcollege_t(jc-1_ik,nabil,anum,bnum,chibnum, tsim),ghcollege_t(jc-1_ik,nabil,anum,bnum,chibnum, tsim))
    allocate(laborh_t(jrnum,nT,nabil,e_num,anum,bnum,tsim),  laborl_t(jrnum,nT,nabil,e_num,anum,bnum,tsim))
    allocate(gh_t(jnum,nT,nabil,e_num,anum,bnum,tsim), gl_t(jnum,nT,nabil,e_num,anum,bnum, tsim))

    vl_t =0.0_rk; vh_t = 0.0_rk; vhcollege_t =0.0_rk !

    !!!!!! calculate the index/weights in advance to save time!!!!!
    call indexcal(nabil, anum, bnum, adnum, bdnum,  r_b, lambda,   bgrid, bdgrid,  adgrid, a0grid, agrid, aindex, ainiindex, aweight, ainiweight, &
        bindex, bweight, payindex, payweight, nopayindex, nopayweight, ad0index, ad0weight)

    do t = tnum,1,-1
        if (t .eq. tnum) then
            vlf_t = 0.0_rk; vhf_t = 0.0_rk; vhcollegef_t =0.0_rk
        else
            vlf_t = vl_t; vhf_t = vh_t; vhcollegef_t = vhcollege_t
        end if


        if (t.gt.tsim) then
            t1 = tsim
        else
            t1 = t
        end if

        !!!!!!! Retirees and Workers!!!!!!!
        call decisionrule(tbar, nabil, eagridj_t(:,:,:,:,t), chibshock, anum, bnum,  chir, psycost_t(1:nabil,1:anum,1:ncost,t),  &
            piel0_t(1:e_num, t),  piec0_t(1:e_num, t),  phi_t(1:nabil, 1:anum, t), bbar_t(t), &
            a0grid,  t,  vhf_t, vlf_t, vhcollegef_t, chidage, lambda, agrid, bgrid,  wh_t(t), wl, &
            piec_t(1:e_num,1:e_num,t), piel_t(1:e_num,1:e_num,t), phid_abil, payindex, payweight, nopayindex, nopayweight,&
            gh_t(1:jnum,1:nT,1:nabil, 1:e_num,1:anum,1:bnum,t1),gl_t(1:jnum,1:nT,1:nabil,1:e_num,1:anum,1:bnum,t1), vh_t, vl_t, &
            pchoicec_t(1:jpaynum,1:nT,1:nabil,1:e_num,1:anum,1:bnum,t1), pchoicenc_t(1:jpaynum,1:nT,1:nabil,1:e_num,1:anum,1:bnum,t1),&
            laborh_t(1:jrnum,1:nT,1:nabil,1:e_num,1:anum,1:bnum,t1), laborl_t(1:jrnum,1:nT,1:nabil,1:e_num,1:anum,1:bnum,t1), &
            edudecision_t(1:nabil,1:anum,1:ncost,1:chibnum, t1), bfchoice_t(1:nabil,1:anum,1:ncost,1:chibnum, t1), bhini_t, vhini_t, &
            vhcollege_t, laborhcollege_t(1:jc-1_ik,1:nabil,1:anum,1:bnum,1:chibnum, t1), ghcollege_t(1:jc-1_ik,1:nabil,1:anum,1:bnum,1:chibnum, t1))

        write(*,*) 'time = ', tnum - t + 1_ik

    end do

    deallocate(vl_t, vh_t, vhf_t, vlf_t, vhcollegef_t, vhcollege_t)

    allocate(stddebtholder_t(tsim),  kagg_t(tsim), lagg_t(tsim),   yagg_t(tsim), totalstudentdebt_t(tsim), avgdebt_t(tsim), &
        avgdebtc_t(tsim), avgdebtd_t(tsim), debtoutstanding_t(tsim) ,educated_t(tsim), noneducated_t(tsim), &
        avgnettutition(tsim), ddebtpercent_t(tsim), totdebt_yagg_t(tsim), totalborrower_t(tsim), avgdebt_yagg_t(tsim), &
        avgtuition_yagg_t(tsim), muh_t(jnum,nT,nabil,e_num,adnum,bdnum), mul_t(jnum,nT,nabil,e_num,adnum,bdnum), mud_t(jnum,nT,nabil,e_num,adnum,bdnum),&
        afsim(jnum,nsim), bfsim(jnum,nsim),edufsim(nsim), colrated_t(tsim), muhcollege_t(jc-1, nabil, adnum, bdnum, chibnum), &
        drateh_t(tsim), drated_t(tsim))

    avgnettutition =0.0_rk
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !       find intitial distribution in 1980
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !initial guess for the distribution
    muh_t = muh_ss
    mul_t = mul_ss
    mud_t = mud_ss
    muhcollege_t = muhcollege_ss


    do t = 1, tsim, +1
        write(*,*) 'entiredistribution for t', t
        ! !!!!! Distribution !!!!!!


        call entiredist(pabil, phid_abil, chibprob, nabil, eagridj_t(:,:,:,:,t), anum, bnum, adnum, bdnum, wp,  probcost_t(1:nabil,1:anum,1:ncost,t), &
            pop, ad0index, ad0weight,  adgrid, edudecision_t(1:nabil,1:anum,1:ncost,1:chibnum,t), piel0_t(1:e_num, t), piec0_t(1:e_num, t),&
            bfchoice_t(1:nabil,1:anum,1:ncost,1:chibnum,t), bdgrid, laborhcollege_t(1:jc-1_ik,1:nabil,1:anum,1:bnum,1:chibnum, t),          &
            laborh_t(1:jrnum,1:nT,1:nabil, 1:e_num,1:anum,1:bnum,t), t, r_b, lambda, &
            pchoicec_t(1:jpaynum,1:nT,1:nabil,1:e_num,1:anum,1:bnum,t), pchoicenc_t(1:jpaynum,1:nT,1:nabil,1:e_num,1:anum,1:bnum,t), &
            ainiindex, ainiweight, ghcollege_t(1:jc-1_ik,1:nabil,1:anum,1:bnum,1:chibnum, t), gh_t(1:jnum,1:nT,1:nabil,1:e_num,1:anum,1:bnum,t), &
            gl_t(1:jnum,1:nT,1:nabil,1:e_num,1:anum,1:bnum,t),&
            aindex, aweight, bindex, bweight, piec_t(1:e_num,1:e_num,t), piel_t(1:e_num,1:e_num,t),   &
            laborl_t(1:jrnum,1:nT,1:nabil,1:e_num,1:anum,1:bnum,t), muhcollege_t, muh_t, mul_t,  mud_t, kagg_t(t),lagg_t(t), yagg_t(t), &
            debtoutstanding_t(t), navg, educated_t(t), noneducated_t(t), stddebtholder_t(t), parentaltransferc, parentaltransfernc, &
            parentaltransfer,debtbyage ,ddebt, fracborrower_t(1:5,t),avgdebt_t(t), avgdebtc_t(t), avgdebtd_t(t),totaldebtbyage_t(1:5,t), &
            ddebtpercent_t(t), totalborrower_t(t), ddebtbyage_t(1:5,t), ddebtmubyage_t(1:5,t), ddebtpercentmu_t(t), colrated_t(t), drateh_t(t), drated_t(t))


        !!!!!Calculate the prices!!!!!
        do ie  = 1, nabil
            do ia = 1, anum
                avgnettutition(t) = avgnettutition(t) + phi_t(ie,ia,t)*pia0(ia)*pabil(ie)
            end do
        end do

        !!!!!yaggdata and yaggpercapitadata is all set up for 2004 for all t!!!!!!!!!.
        totdebt_yagg_t(t) = -debtoutstanding_t(t)/yagg_t(t)
        totalstudentdebt_t(t) = totdebt_yagg_t(t)*yaggdata
        debtbyage_t(1:5,t) = debtbyage(1:5)/debtoutstanding_t(t)
        avgdebt_yagg_t(t) = -avgdebt_t(t)/yagg_t(t)
        avgtuition_yagg_t(t) = avgnettutition(t)/ yagg_t(t)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                   !
        !                targets and moments                !
        !                                                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        write(*,*) 'Solving distribution for year ', t + 1983
        if (t.lt.19) then !write only for 1997


            write(*,*) '2004 federal loan supply (college board)', 4148.72 , 'student loan borrowing limit', -0.25*(bbar_t(t)*(gdp_data))/yagg_t(t)


            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(*,*) '!                                                   !'
            write(*,*) '!      Targets Identifying 4 Parameters             !'
            write(*,*) '!                                                   !'
            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(*,*)
            write(*,*) '                                            targets                                       data            model'
            write(*,*)
            write(*,'(1x, a, f8.4)')    ' A1.  K/Y                                                           :    3.00   ',    kagg_t(t)/yagg_t(t)
            write(*,*)
            write(*,'(1x, a, f8.4)')    ' A2.  average labor supply                                          :    0.3333 ',    navg
            write(*,*)
            write(*,'(1x, a, f7.3)')    ' A3.  Frisch elasticity of male labor supply                        :   0.50    ',    (1/eta)
            write(*,*)
            write(*,'(1x, a, f8.4)')    ' A4.  Student loan borrowing limit(unsubsidized) $23,000            :    0.74   ',    -bbar_t(t)/yagg_t(t)
            write(*,*)
            write(*,'(1x, a, f8.4)')    ' A5.  Graduation rate in 1997(NLSY97 cohort/ whole pop)             :    0.2294 ',    educated_t(t)
            write(*,*)
            write(*,'(1x, a, f8.4)')    ' A6.  Graduation rate in 1997(NLSY97 cohort/ j = 1:jc-1)            :    0.2294 ',    colrated_t(t)
            write(*,*)
            write(*,'(1x, a, f8.4)')    ' A7.  % students graduating with student deb                        :    0.60   ',    stddebtholder_t(t)
            write(*,*)
            write(*,'(1x, a, f8.4)')    ' A8.  Average student debt at graduation($11,348/24739.68)          :    0.4587 ',   avgdebt_yagg_t(t)
            write(*,*)
            write(*,'(1x, a, f8.4)')    ' A9.  outstanding std debt as a fraction of gdp (179.77B/5553.87B)  :    0.0324 ',    totdebt_yagg_t(t)
            write(*,*)
            write(*,'(1x,2(a,f8.4),a)') ' A10.  fraction of delinquent student debt                          :    _____  ',   ddebtpercent_t(t)
            write(*,*)
            write(*,'(1x, a, f8.4)')    ' A11.  Average net college tuition($7,737/24739.68)                 :    0.3127 ',    avgtuition_yagg_t(t)
            write(*,*)
            write(*,*)
            write(*,'(1x, a, f8.4)')    ' Average parental transfer: $4,719.56 (0.1908)                                 ', parentaltransfer/yagg_t(t)
            write(*,'(1x, a, f8.4)')    ' Average parental transfer for college: $9,722.83 (0.3930)                     ', parentaltransferc/yagg_t(t)
            write(*,'(1x, a, f8.4)')    ' Average parental transfer for non-college: $3,230.14 (0.1306)                 ', parentaltransfernc/yagg_t(t)
            write(*,*)
            write(*,'(1x, a, f8.4)')    ' Average debt at graduation for all           : $5,229.176 (0.1680)            : ', -avgdebt_t(t)/yagg_t(t)
            write(*,'(1x, a, f8.4)')    ' Average debt at graduation for college       : $9,493.583 (0.3050)            : ', -avgdebtc_t(t)/yagg_t(t)
            write(*,'(1x, a, f8.4)')    ' Average debt at graduation for dropout       : $2,618.589 (0.0841)            : ', -avgdebtd_t(t)/yagg_t(t)    
            write(*,*)            
            write(*,'(1x, a, f8.4)')    ' fraction of delinquent debt over total debt                                   ', ddebtpercent_t(t)
            write(*,'(1x, a, f8.4)')    ' fraction of delinquent borrowers over total borrowers                         ', ddebtpercentmu_t(t)
            write(*,*)
            write(*,'(1x, a, f8.4)')    ' delinqunecy rate for all workers                                            : ',  ddebtpercentmu_t(t)
            write(*,'(1x, a, f8.4)')    ' delinqunecy rate for college graduate workers                               : ',  drateh_t(t)
            write(*,'(1x, a, f8.4)')    ' delinqunecy rate for college dropouts workers                               : ',  drated_t(t)              
            write(*,*)
            !write(*,'(1x,5(a,f8.2),a)')    'debt/gdp age under 30 (Data: 0.10, Model:', totaldebtbyage_t(1,t), ')'
            !write(*,'(1x,5(a,f8.2),a)')    'debt/gdp age 30-39    (Data: 0.08, Model:', totaldebtbyage_t(2,t), ')'
            !write(*,'(1x,5(a,f8.2),a)')    'debt/gdp age 40-49    (Data: 0.03, Model:', totaldebtbyage_t(3,t), ')'
            !write(*,'(1x,5(a,f8.2),a)')    'debt/gdp age 50-59    (Data: 0.02, Model:', totaldebtbyage_t(4,t), ')'
            !write(*,'(1x,5(a,f8.2),a)')    'debt/gdp age 60+      (Data: 0.00, Model:', totaldebtbyage_t(5,t), ')'
            !write(*,*)
            !!write(*,*)
            !write(*,'(1x, a, f8.2,a)')     'fraction of borrower age under 30 (Data: 0.50, Model:', fracborrower_t(1,t),')'
            !write(*,'(1x, a, f8.2,a)')     'fraction of borrower age 30-39    (Data: 0.25, Model:', fracborrower_t(2,t),')'
            !write(*,'(1x, a, f8.2,a)')     'fraction of borrower age 40-49    (Data: 0.14, Model:', fracborrower_t(3,t),')'
            !write(*,'(1x, a, f8.2,a)')     'fraction of borrower age 50-59    (Data: 0.09, Model:', fracborrower_t(4,t),')'
            !write(*,'(1x, a, f8.2,a)')     'fraction of borrower age 60+      (Data: 0.03, Model:', fracborrower_t(5,t),')'
            !write(*,*)
            write(*,'(1x, a, f12.0)')    'Total number of borrowers', totalborrower_t(t)*totalpopdata!
            write(*,*)
            write(*,'(1x, a, f12.6)') 'population hitting blow', sum(muhcollege_t(:,:,:,1,:))+ sum(muh_t(:,:,:,:,:,1))+ sum(mul_t(:,:,:,:,:,1))+ sum(mud_t(:,:,:,:,:,1))
            write(*,'(1x, a, f12.6)') 'population hitting ahigh', sum(muhcollege_t(:,:,adnum,:,:))+sum(muh_t(:,:,:,:,adnum,:))+sum(mul_t(:,:,:,:,adnum,:))+sum(mud_t(:,:,:,:,adnum,:))


            hitblow_t(t) = sum(muhcollege_t(:,:,:,1,:))+sum(muh_t(:,:,:,:,:,1))+ sum(mul_t(:,:,:,:,:,1))+ sum(mud_t(:,:,:,:,:,1))
            hitahigh_t(t) = sum(muhcollege_t(:,:,adnum,:,:))+sum(muh_t(:,:,:,:,adnum,:))+sum(mul_t(:,:,:,:,adnum,:))+sum(mud_t(:,:,:,:,adnum,:))
        end if

    end do




    if ( itransim .eq. 1_ik) then
    write(*,*) 'transimulation before allocation'
        allocate( lcbsim(jnum), lcasim(jnum), edulclaborsim(jnum,3), edulcbsim(jnum,3), edulcasim(jnum,3), defaultsim(jnum,3),debtbyagesim(5,tsim), &
            lcearnsim(jnum), edulcearnsim(jnum,3), colrate_t(3,4,tsim), colrate_abil_t(3,tsim), colrate_wealth_t(4,tsim), abilratio_t(tsim), wealthratio_t(tsim),&
            avgability(tsim), avgpt(tsim), avgdebt_sim(tsim), stddebtholder_sim(tsim), totborrower_sim(tsim), totdebt_sim(tsim), colrate_aggsim_t(tsim),&
            ddebtratio_t(tsim), ddebtmuratio_t(tsim), ticagg_t(tsim), excagg_t(tsim), avgnettui_t(tsim), const_t(tsim), const1_t(tsim), colrate_aggsim1_t(tsim))
    write(*,*) 'transimulation after allocation'
        call transimulate(totborrower_sim, eagridj_t(:,:,:,:,1:tsim), abilgrid, pabil, phid_abil, r_b, lambda, nabil,nsim, anum, bnum, yagg_t(1:tsim), egridc, egridl,  &
            piec0_t(1:e_num,1:tsim), piel0_t(1:e_num,1:tsim), piec_t(1:e_num,1:e_num, 1:tsim), piel_t(1:e_num,1:e_num, 1:tsim), &
            a0grid, piea0_t(1:nabil,1:anum,1:tsim),  laborl_t(1:jrnum,1:nT,1:nabil, 1:e_num,1:anum,1:bnum, 1:tsim), &
            laborh_t(1:jrnum,1:nT,1:nabil,1:e_num,1:anum,1:bnum, 1:tsim), laborhcollege_t(1:jc-1_ik,1:nabil,1:anum,1:bnum,1:chibnum, 1:tsim), &
            gl_t(1:jnum,1:nT,1:nabil,1:e_num,1:anum,1:bnum, 1:tsim), gh_t(1:jnum,1:nT,1:nabil,1:e_num,1:anum,1:bnum, 1:tsim), ghcollege_t(1:jc-1_ik,1:nabil,1:anum,1:bnum,1:chibnum, 1:tsim), &
            pchoicec_t(1:jpaynum,1:nT,1:nabil,1:e_num,1:anum,1:bnum, 1:tsim), pchoicenc_t(1:jpaynum,1:nT,1:nabil,1:e_num,1:anum,1:bnum, 1:tsim), &
            bfchoice_t(1:nabil,1:anum,1:ncost, 1:chibnum,1:tsim), edudecision_t(1:nabil,1:anum,1:ncost, 1:chibnum,1:tsim), agrid, bgrid, &
            wh_t(1:tsim), wl, pi,  phi_t(1:nabil,1:anum,1:tsim), probcost_t(1:nabil,1:anum,1:ncost, 1:tsim), &
            psycost_t(1:nabil,1:anum,1:ncost,  1:tsim),  asim, bsim, edusim, idxprod, prodsim, paysim, idxability,&
            abilsim, colrate_t, colrate_abil_t, colrate_wealth_t, abilratio_t, wealthratio_t, debtbyagesim, avgability, avgpt, &
            avgdebt_sim, stddebtholder_sim, totdebt_sim, colrate_aggsim_t, colrate_aggsim1_t, ddebtratio_t, ddebtmuratio_t, ticagg_t, excagg_t,avgnettui_t, const_t, idxchib, &
            chibprob, chibshock, const1_t)

    end if

    call hdf5_openf(datafile1, fileid)
    call hdf5_write(wh_t, fileid, 'wh_t')
    call hdf5_write(phi_t, fileid, 'phi_t')
    call hdf5_write(bbar_t, fileid, 'bbar_t')
    call hdf5_write(totalstudentdebt_t, fileid, 'totalstudentdebt_t')
    call hdf5_write(stddebtholder_t, fileid, 'stddebtholder_t')
    call hdf5_write(avgdebt_t, fileid, 'avgdebt_t')
    call hdf5_write(avgnettutition, fileid, 'avgnettutition')
    call hdf5_write(educated_t, fileid, 'educated_t')
    call hdf5_write(colrated_t, fileid, 'colrated_t')
    !call hdf5_write(totalborrower_t, fileid, 'totalborrower_t')
    !call hdf5_write(ddebtpercent_t, fileid, 'ddebtpercent_t')
    call hdf5_write(totdebt_yagg_t, fileid, 'totdebt_yagg_t')
    !call hdf5_write(avgdebt_yagg_t, fileid, 'avgdebt_yagg_t')
    call hdf5_write(yagg_t, fileid, 'yagg_t')
    call hdf5_write(fracborrower_t, fileid, 'fracborrower_t')
    call hdf5_write(totaldebtbyage_t, fileid, 'totaldebtbyage_t')
    call hdf5_write(debtbyage_t, fileid, 'debtbyage_t')
    call hdf5_write(ddebtbyage_t, fileid, 'ddebtbyage_t')
    call hdf5_write(ddebtmubyage_t, fileid, 'ddebtmubyage_t')
    !call hdf5_write(ddebtpercentmu_t, fileid, 'ddebtpercentmu_t')
    call hdf5_write(colrate_t, fileid, 'colrate_t')
    call hdf5_write(colrate_aggsim_t, fileid, 'colrate_aggsim_t')
    call hdf5_write(colrate_aggsim1_t, fileid, 'colrate_aggsim1_t')
    call hdf5_write(colrate_abil_t, fileid, 'colrate_abil_t')
    call hdf5_write(colrate_wealth_t, fileid, 'colrate_wealth_t')
    call hdf5_write(abilratio_t, fileid, 'abilratio_t')
    call hdf5_write(wealthratio_t, fileid, 'wealthratio_t')
    call hdf5_write(debtbyagesim, fileid, 'debtbyagesim')
    call hdf5_write(avgability, fileid, 'avgability')
    call hdf5_write(avgpt, fileid, 'avgpt')
    call hdf5_write(avgdebt_sim, fileid, 'avgdebt_sim')
    call hdf5_write(stddebtholder_sim, fileid, 'stddebtholder_sim')
    call hdf5_write(totborrower_sim, fileid, 'totborrower_sim')
    call hdf5_write(totdebt_sim, fileid, 'totdebt_sim')
    call hdf5_write(ddebtratio_t, fileid, 'ddebtratiosim_t')
    call hdf5_write(ddebtmuratio_t, fileid, 'ddebtmuratiosim_t')
    call hdf5_write(ticagg_t, fileid, 'ticagg_t')
    call hdf5_write(excagg_t, fileid, 'excagg_t')
    call hdf5_write(const_t, fileid, 'const_t')
    call hdf5_write(const1_t, fileid, 'const1_t')
    call hdf5_write(avgnettui_t, fileid, 'avgnettui_t')
    call hdf5_write(hitblow_t, fileid, 'hitblow_t')
    call hdf5_write(hitahigh_t, fileid, 'hitahigh_t')
    call hdf5_write(edudecision_t, fileid, 'edudecision_t')
    call hdf5_write(bfchoice_t, fileid, 'bfchoice_t')
    if ( iwelfare .eq. 1_ik) then
        call hdf5_write(welvalue, fileid, 'welvalue')
    end if

    call hdf5_closef(fileid)


    end subroutine transition


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! simulate the data in the model
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine transimulate(totborrower, eagridj_t,abilgrid, pabil, phid_abil, r_b, lambda, nabil, nsim, anum, bnum, yagg,              &
        egridc, egridl,piec0, piel0, piec, piel, a0grid,  piea0_t,  laborl, laborh, laborhcollege, gl, gh, ghcollege,                       &
        pchoicec, pchoicenc, bfchoice, edudecision, agrid, bgrid, wh, wl, pi,  netphi, probcost, psycost,                                   &
        asim_ss, bsim_ss, edusim_ss, idxprod_ss, prodsim_ss, paysim_ss, idxability_ss, abilsim_ss,                                          &
        colrate_ini_d,colrate_abil_ini, colrate_wealth_ini, abilratio, wealthratio, debtbyagesim, &
        avgability, avgpt, avgdebt, stddebtholder, totdebt, colrateagg,colrateaggd,ddebtratio, ddebtmuratio, ticagg_t, &
        excagg_t, avgnettui, const_t, idxchib_ss, chibprob, chibshock, const1_t)

    integer(ik):: nabil, nsim,anum, bnum,  ic, ie, ia,   j, ip,  ipf,  t,  ix, minplace(1), &
        olgo,  iab,  ichib,    iss, e, t1, j1

    integer(ik)::  pchoicec(jpaynum,nT, nabil, e_num, anum, bnum, tsim), pchoicenc(jpaynum,nT,nabil, e_num, anum, bnum, tsim), &
        idxprod_ss(jnum, nsim), paysim_ss(jnum, nsim), idxability_ss(nsim), idxchib_ss(nsim), &
        pchoicecin(anum,bnum), pchoicencin(anum,bnum)

    integer(ik), allocatable ::  idxprod(:,:,:), idxprodh(:,:), idxprodl(:,:), paysim(:,:,:), idxiniwealth(:, :), wealthcut(:,:), abilitycut(:,:),  &
        idxcost(:,:,:),  idxabil_n(:), idxabil(:,:,:), idxchib(:,:,:), idxchib_n(:)

    real(rk), allocatable:: bsim(:,:,:), asim(:,:,:), edusim(:,:,:),  prodsim(:,:,:), prodsimh(:,:), prodsiml(:,:),  laborsim(:,:,:),&
        earnsim(:,:,:),  iniwealth(:, :),simcost(:,:,:), abilsim_n(:), csim(:,:,:), abilsim(:,:,:),  simchib_n(:)

    real(rk)::egridc(e_num), egridl(e_num), a0grid(nabil, anum), agrid(anum), bgrid(bnum), wh(tsim), phid_abil(nabil), &
        laborl(jrnum,nT,nabil, e_num,anum,bnum, tsim), laborh(jrnum,nT,nabil, e_num,anum,bnum, tsim),  piec(e_num,e_num,tsim), piel(e_num,e_num,tsim), &
        yagg(tsim), netphi(nabil,anum, tsim), tuisim(jnum, nsim, tsim), gl(jnum,nT,nabil,e_num,anum,bnum, tsim), gh(jnum,nT,nabil,e_num,anum,bnum, tsim),  &
        piec0(e_num, tsim), piel0(e_num, tsim), piea0_t(nabil, anum,tsim), piea0temp(anum), probcosttemp(ncost),  probcost(nabil,anum,ncost, tsim),     &
        psycost(nabil,anum,ncost, tsim), bfchoice(nabil,anum,ncost,chibnum, tsim), edudecision(nabil,anum,ncost,chibnum, tsim), asim_ss(jnum, nsim), bsim_ss(jnum, nsim), &
        edusim_ss(jnum, nsim), prodsim_ss(jnum, nsim), abilsim_ss(jnum, nsim), ddebt(jnum, nsim, tsim), ddebtmu(jnum, nsim, tsim), lambda(nT), phigrid(2), phidist(2), &
        chibprob(chibnum), chibshock(chibnum), laborhcollege(jc-1_ik,nabil,anum,bnum,chibnum, tsim), ghcollege(jc-1_ik,nabil,anum,bnum,chibnum, tsim), &
        ghtemp(anum,bnum), lhtemp(anum,bnum), lltemp(anum,bnum)

    real(rk):: avgdebt(tsim), dchoicesim(jnum,nsim,tsim), avglaborini(3,4, tsim),avgasimini(3,4, tsim),avgbsimini(3,4, tsim), debtbyagesim(5, tsim), &
        avgnetphi(3,4, tsim),  colrateagg(tsim), colrateaggd(tsim), stddebtholder(tsim), tuipt(4, tsim), tuiabil(3, tsim),   avgnettui(tsim), &
        abilratio(tsim), wealthratio(tsim),  avgability(tsim),  avgpt(tsim),  totdebt(tsim), &
        colrate_ini(3,4, tsim), colrate_ini_d(3,4, tsim),   colrate_abil_ini(3, tsim), colrate_wealth_ini(4, tsim), ddebtratio(tsim), ddebtmuratio(tsim),&
        ticagg_t(tsim),  excagg_t(tsim),  avgearnini(3,4,tsim), earnabil(3,tsim), earnpt(4,tsim), earncutval_wealth(4,tsim), &
        wealthcutval(4,tsim), abilitycutval(3,tsim), earncutval_abil(3,tsim), const_t(tsim), &
        drateearn(4,tsim), dsample(4,tsim), bearnq(4,tsim), eagridj_t(jnum, nabil, e_num, 2, tsim),abilgrid(nabil), pabil(nabil),&
        totborrower(tsim), drate(tsim),  const1_t(tsim)



    real(rk):: afval, aval, bval, elevel, wl, pi, r_b,  chib, bfval


    integer(hid_t):: fileid
    character(30):: datafile,  datafile1, datestring, timestring


    intent(in)::  r_b, lambda, nsim, anum, bnum, yagg,  egridc, egridl, a0grid, laborhcollege, ghcollege, &
        laborl, laborh, gl, gh,  pchoicec, pchoicenc, agrid, bgrid,  wh, wl, pi, netphi, piec0, piel0, &
        piea0_t, probcost, psycost, bfchoice, edudecision, asim_ss, bsim_ss, edusim_ss, idxprod_ss, prodsim_ss, piec, &
        piel, paysim_ss, eagridj_t, idxability_ss, abilsim_ss, nabil,abilgrid, pabil, phid_abil, idxchib_ss,  chibprob, chibshock

    intent(out)::  colrate_ini_d, colrate_abil_ini, colrate_wealth_ini, abilratio, wealthratio, debtbyagesim, avgability, avgpt, &
        avgdebt, stddebtholder, totdebt, colrateagg,colrateaggd,ddebtratio, ticagg_t, excagg_t, avgnettui, const_t, totborrower, ddebtmuratio, const1_t


    call date_and_time(datestring, timestring)
    datafile = 'KK_stddebt_'
    olgo = len_trim(datafile)
    datafile(olgo + 1:olgo+4) = timestring(1:4)
    datafile1 = datafile(1:olgo+4)
    datafile1(olgo+5:olgo+17) = '_simultran.h5'

    allocate(bsim(jnum, nsim, tsim), asim(jnum, nsim, tsim), edusim(jnum,nsim, tsim), paysim(jnum,nsim, tsim), idxchib(jnum, nsim,tsim),&
        idxprod(jnum,nsim, tsim), idxprodl(jnum, nsim), idxprodh(jnum,  nsim), prodsim(jnum,nsim, tsim), prodsimh(jnum, nsim), prodsiml(jnum,  nsim), &
        earnsim(jnum,nsim, tsim), laborsim(jnum, nsim, tsim), csim(jnum, nsim, tsim))
    allocate( idxabil_n(nsim), abilsim_n( nsim), idxiniwealth(nabil, nsim), iniwealth(nabil, nsim), idxcost(nabil,anum,nsim), &
        simcost(nabil,anum,nsim), idxabil(jnum, nsim, tsim), abilsim(jnum, nsim, tsim), idxchib_n(nsim), simchib_n(nsim))

    write(*,*) 'transimulate starts'
    asim= 0.0_rk; bsim = 0.0_rk;  earnsim = 0.0_rk; laborsim =0.0_rk; avglaborini=0.0_rk; avgasimini=0.0_rk; avgbsimini=0.0_rk; avgnetphi=0.0_rk;
    tuipt=0.0_rk; tuiabil=0.0_rk; tuisim = 0.0_rk; edusim = -1.0_rk; ddebt = 0.0_rk; ddebtmu = 0.0_rk; csim = 0.0_rk; dchoicesim=0.0_rk;
    idxabil = 0_ik; abilsim = 0.0_rk;


    phigrid(1) = 0.5_rk; phigrid(2) = 1.0_rk

    ! initial period state values are from the steady-state
    asim(:,:,1) = asim_ss(:,:)
    bsim(:,:,1) = bsim_ss(:,:)
    prodsim(:,:,1) = prodsim_ss(:,:)
    idxprod(:,:,1) = idxprod_ss(:,:)
    edusim(:,:,1) = edusim_ss(:,:)
    paysim(:,:,1) = paysim_ss(:,:)

    do j = 1, jnum
        idxchib(j,:,1) = idxchib_ss(:)
        idxabil(j,:,1) = idxability_ss(:)
        abilsim(j,:,1) = abilsim_ss(j,:)
    end do

    do t = 1, tsim
        write(*,*) 'starting simcal for t = ', 1983 + t

        do ic = 1,nsim
            do j= 1,jnum

                aval = asim(j,ic,t); bval = bsim(j,ic,t); elevel = edusim(j,ic,t)
                ie = idxprod(j,ic,t); ip = paysim(j,ic,t); iab = idxabil(j,ic,t)
                ichib = idxchib(j,ic,t); chib = chibshock(ichib);
                
                if (elevel.eq.1.0_rk) then
                    if (j.le.jc-1_ik) then
                        ghtemp = ghcollege(j,iab,:,:,ichib, t)
                        lhtemp = laborhcollege(j,iab,:,:,ichib,t)
                        pchoicecin =1_ik!
                    elseif (j.le.jpaynum) then
                        ghtemp = gh(j,ip,iab,ie,:,:,t)
                        lhtemp = laborh(j,ip,iab,ie,:,:,t)
                        pchoicecin = pchoicec(j,ip,iab,ie,:,:,t)
                    elseif (j.le.jrnum) then
                        ghtemp = gh(j,ip,iab,ie,:,:,t)
                        lhtemp = laborh(j,ip,iab,ie,:,:,t)
                        pchoicecin =1_ik
                    else
                        ghtemp = gh(j,ip,iab,ie,:,:,t)
                        lhtemp = 0.0_rk
                        pchoicecin = 1_ik
                    end if
                else
                    if (j.le.jpaynum) then
                        lltemp = laborl(j,ip,iab,ie,:,:,t)
                        pchoicencin = pchoicenc(j,ip,iab,ie,:,:,t)
                    elseif (j.le.jrnum) then
                        lltemp = laborl(j,ip,iab,ie,:,:,t)
                        pchoicencin = 1_ik
                    else
                        lltemp = 0.0_rk
                        pchoicencin = 1_ik
                    end if
                end if

                call simcal( j, ip, anum, bnum, netphi(iab,:,t), pchoicencin, pchoicecin, chib, aval, bval, elevel, wl, wh(t), r_b, &
                    lambda(ip), a0grid(iab,:), agrid, bgrid, tuisim(j,ic,t), gl(j,ip,iab,ie,:,:, t), lltemp, ghtemp, &
                    lhtemp, eagridj_t(j,iab,ie,:, t), ipf, earnsim(j,ic,t), laborsim(j,ic,t), afval,  dchoicesim(j,ic,t),  &
                    csim(j,ic,t), bfval, ddebt(j,ic,t), ddebtmu(j,ic,t))

                if ( j.lt. jnum .and. t .lt. tsim) then
                    asim(j+1,ic,t+1) = afval
                    bsim(j+1,ic,t+1) = bfval
                    paysim(j+1,ic,t+1) = ipf
                    idxabil(j+1,ic,t+1) = iab
                    idxchib(j+1,ic,t+1) = ichib
                    abilsim(j+1,ic,t+1) = abilsim(j,ic,t)
                    edusim(j+1,ic, t+1) = elevel
                    idxprod(j+1,ic,t+1) = ie
                    prodsim(j+1,ic,t+1) = prodsim(j,ic,t)
                    tuisim(j+1,ic,t+1) = tuisim(j,ic,t)

                    !if ( j .eq. jc-1 .and. elevel .eq. 1.0_rk) then
                    !    phidist =0.0_rk; phidist(1) = phid_abil(iab); phidist(2) = 1.0_rk-phid_abil(iab);
                    !    call iniwealthsim(2_ik, 1, phigrid, phidist, minplace(1), edusim(jc,ic,t+1)) ! drop out case
                    !    if (edusim(j+1,ic, t+1) .eq. 0.5_rk) then
                    !        e = 0_ik
                    !        call prodsimulate(e,  jnum-jc+1,  e_num, 1, egridc, piec(ie,1:e_num,t+1), piel(ie,1:e_num,t+1), piec0(:,t+1), piel0(:,t+1),  &
                    !            idxprodh, prodsimh, idxprodl, prodsiml)
                    !        do t1 =t+1, tsim
                    !            j1 = t1-t+jc-1
                    !            if (j1 .le. jnum ) then
                    !                prodsim(j1,1:nsim,t1) = prodsiml(j1,1:nsim)
                    !                idxprod(j1,1:nsim,t1) = idxprodl(j1,1:nsim)
                    !            end if
                    !        end do
                    !
                    !    end if
                    !end if
                    if ( elevel .le. 0.5_rk) then
                        if ( j .le. jrnum) then
                            call iniwealthsim(e_num, 1, egridl, piel(ie,1:e_num,t+1), idxprod(j+1,ic, t+1), prodsim(j+1,ic,t+1))
                        end if
                    else
                        tuisim(j+1,ic,t+1) = tuisim(j,ic,t)
                        if ( jc .eq. jd) then
                            if ( j .eq. jc-1) then
                                phidist =0.0_rk; phidist(1) = phid_abil(iab); phidist(2) = 1.0_rk-phid_abil(iab);
                                call iniwealthsim(2_ik, 1, phigrid, phidist, minplace(1), edusim(j+1,ic,t+1)) ! drop out case
                                if (edusim(j+1,ic, t+1) .eq. 0.5_rk) then
                                    call iniwealthsim(e_num, 1, egridl, piel0(:,t+1), idxprod(j+1, ic,t+1), prodsim(j+1, ic,t+1))
                                else
                                    call iniwealthsim(e_num, 1, egridc, piec0(:,t+1), idxprod(j+1, ic,t+1), prodsim(j+1, ic,t+1))
                                end if
                            elseif (j .ge. jc .and. j.le. jrnum) then
                                call iniwealthsim(e_num, 1, egridc, piec(ie,1:e_num,t+1), idxprod(j+1,ic, t+1), prodsim(j+1,ic,t+1))
                            end if
                        else
                             if ( j .eq. jd-1) then
                                phidist =0.0_rk; phidist(1) = phid_abil(iab); phidist(2) = 1.0_rk-phid_abil(iab);
                                call iniwealthsim(2_ik, 1, phigrid, phidist, minplace(1), edusim(j+1,ic,t+1)) ! drop out case
                                if (edusim(j+1,ic, t+1) .eq. 0.5_rk) then
                                    call iniwealthsim(e_num, 1, egridl, piel0(:,t+1), idxprod(j+1, ic,t+1), prodsim(j+1, ic,t+1))
                                end if
                            elseif (j .eq. jc-1 .and. elevel .eq. 1.0_rk) then
                                 call iniwealthsim(e_num, 1, egridc, piec0(:,t+1), idxprod(j+1, ic,t+1), prodsim(j+1, ic,t+1))
                            elseif (j .ge. jc .and. j.le. jrnum) then
                                call iniwealthsim(e_num, 1, egridc, piec(ie,1:e_num,t+1), idxprod(j+1,ic, t+1), prodsim(j+1,ic,t+1))
                            end if
                        end if
                        
                    end if


                end if

            end do
        end do
        !    !$omp end parallel do
        if ( t+1 .le. tsim) then
            !simulate initial wealth and education cost
            call iniwealthsim(nabil, nsim, abilgrid, pabil, idxabil_n, abilsim_n)

            call iniwealthsim(chibnum, nsim, chibshock(1:chibnum), chibprob, idxchib_n(1:nsim), simchib_n(1:nsim))

            do ie = 1,nabil
                piea0temp(1:anum) = piea0_t(ie, 1:anum,t)/sum(piea0_t(ie, 1:anum,t))
                call iniwealthsim(anum, nsim, a0grid(ie, 1:anum), piea0temp, idxiniwealth(ie, 1:nsim), iniwealth(ie, 1:nsim))
            end do

            do ie = 1,nabil
                do ia = 1,anum
                    probcosttemp(1:ncost) = probcost(ie,ia,1:ncost,t+1)/sum(probcost(ie,ia,1:ncost,t+1))
                    call iniwealthsim(ncost, nsim, psycost(ie,ia,1:ncost,t+1), probcosttemp, idxcost(ie,ia,1:nsim), simcost(ie,ia,1:nsim))
                end do
            end do

            !decide education level
            do ic = 1,nsim
                iab = idxabil_n(ic); ia = idxiniwealth(iab, ic); ix = idxcost(iab, ia, ic)
                ichib = idxchib_n(ic);

                edusim(1,ic,t+1) = edudecision(iab,ia, ix,ichib, t+1);
                idxabil(1,ic,t+1) = iab; abilsim(1,ic,t+1) = abilgrid(iab)
                idxchib(1,ic,t+1) = ichib

                if (edusim(1,ic,t+1) .eq. 0.0_rk) then
                    asim(1,ic,t+1) = iniwealth(iab,ic)
                    bsim(1,ic,t+1) = 0.0_rk
                    paysim(1,ic,t+1) = nT
                    call iniwealthsim(e_num, 1, egridl, piel0(:,t+1), idxprod(1, ic,t+1), prodsim(1, ic,t+1))
                else !
                    bsim(1,ic,t+1) = bfchoice(iab,ia,ix,ichib,t+1)
                    if (bsim(1,ic,t+1) .gt. -1.0E-003_rk) then
                        bsim(1,ic,t+1) = 0.0_rk
                    end if
                    asim(1,ic,t+1) = iniwealth(iab,ic)
                    tuisim(1,ic,t+1) = netphi(iab,ia,t+1);
                    paysim(1,ic,t+1) = 1_ik
                    idxprod(1,ic,t+1) = 1_ik
                    prodsim(1,ic,t+1) =egridc(idxprod(1,ic,t+1))
                end if

            end do
        end if
        
        
        ! make new borns in t+1
        !if ( t+1 .le. tsim) then
        !
        !    !simulate initial wealth and education cost
        !    call iniwealthsim(nabil, nsim, abilgrid, pabil, idxabil_n, abilsim_n)
        !
        !    call iniwealthsim(chibnum, nsim, chibshock(1:chibnum), chibprob, idxchib_n(1:nsim), simchib_n(1:nsim))
        !
        !    do ie = 1,nabil
        !        piea0temp(1:anum) = piea0_t(ie, 1:anum,t)/sum(piea0_t(ie, 1:anum,t))
        !        call iniwealthsim(anum, nsim, a0grid(ie, 1:anum), piea0temp, idxiniwealth(ie, 1:nsim), iniwealth(ie, 1:nsim))
        !    end do
        !
        !    do ie = 1,nabil
        !        do ia = 1,anum
        !            probcosttemp(1:ncost) = probcost(ie,ia,1:ncost,t+1)/sum(probcost(ie,ia,1:ncost,t+1))
        !            call iniwealthsim(ncost, nsim, psycost(ie,ia,1:ncost,t+1), probcosttemp, idxcost(ie,ia,1:nsim), simcost(ie,ia,1:nsim))
        !        end do
        !    end do
        !
        !    !simulate the intial productivity - ! HK???
        !    e = 1_ik;
        !    call prodsimulate(e, jnum,  e_num, nsim, egridc, piec(ie,1:e_num,t+1), piel(ie,1:e_num,t+1), piec0(:,t+1), piel0(:,t+1),  &
        !        idxprodh, prodsimh, idxprodl, prodsiml)
        !
        !
        !    !decide education level
        !    do ic = 1,nsim
        !        iab = idxabil_n(ic); ia = idxiniwealth(iab, ic); ix = idxcost(iab, ia, ic)
        !        ichib = idxchib_n(ic);
        !
        !        edusim(1,ic,t+1) = edudecision(iab,ia, ix,ichib, t+1);
        !        idxabil(1,ic,t+1) = iab; abilsim(1,ic,t+1) = abilgrid(iab)
        !        idxchib(1,ic,t+1) = ichib
        !
        !        if (edusim(1,ic,t+1) .eq. 0.0_rk) then
        !            asim(1,ic,t+1) = iniwealth(iab,ic)
        !            bsim(1,ic,t+1) = 0.0_rk
        !            paysim(1,ic,t+1) = nT
        !            do t1 =t+1, tsim
        !                j1 = t1-t
        !                if (j1 .le. jnum ) then
        !                    prodsim(j1,1:nsim,t1) = prodsiml(j1,1:nsim)
        !                    idxprod(j1,1:nsim,t1) = idxprodl(j1,1:nsim)
        !                end if
        !            end do
        !
        !        else !
        !            bsim(1,ic,t+1) = bfchoice(iab,ia,ix,ichib,t+1)
        !            if (bsim(1,ic,t+1) .gt. -1.0E-003_rk) then
        !                bsim(1,ic,t+1) = 0.0_rk
        !            end if
        !            asim(1,ic,t+1) = iniwealth(iab,ic)
        !            tuisim(1,ic,t+1) = netphi(iab,ia,t+1);
        !            paysim(1,ic,t+1) = 1_ik
        !
        !            do t1 =t+1, tnum
        !                j1 = t1-t
        !                if (j1 .le. jnum ) then
        !                    prodsim(j1,1:nsim,t1) = prodsimh(j1,1:nsim)
        !                    idxprod(j1,1:nsim,t1) = idxprodh(j1,1:nsim)
        !                end if
        !            end do
        !
        !
        !        end if
        !
        !    end do
        !end if
    end do


    !calculate deafult rate over earnings level
    drateearn = 0.0_rk; dsample = 0.0_rk; bearnq = 0.0_rk; drate = 0.0_rk
    allocate( wealthcut(nsim, tsim), abilitycut(nsim, tsim))
    do t = 1, tsim
        iss = 1_ik
        call sim_moments(phid_abil, totborrower(t), t, iss, nabil, nsim, yagg(t), ddebt(:,:,t), ddebtmu(:,:,t),                      &
            tuisim(:,:,t), bsim(:,:,t), asim(:,:,t), laborsim(:,:,t), earnsim(:,:,t), idxabil(1,:,t), abilsim(:,:,t), edusim(:,:,t), dchoicesim(:,:,t),     &
            abilitycut(:,t), wealthcut(:,t), earncutval_abil(:,t), wealthcutval(:,t), earncutval_wealth(:,t), avglaborini(:,:,t), avgearnini(:,:,t), earnpt(:,t),    &
            earnabil(:,t), avgability(t), drateearn(:,t), debtbyagesim(:,t), colrate_ini(:,:,t), colrate_ini_d(:,:,t), colrate_abil_ini(:,t), colrate_wealth_ini(:,t), abilratio(t), wealthratio(t),       &
            avgdebt(t), avgpt(t), stddebtholder(t), totdebt(t), colrateagg(t),  colrateaggd(t), ddebtratio(t), ddebtmuratio(t), ticagg_t(t), excagg_t(t), avgnettui(t), const_t(t), const1_t(t), abilitycutval(:,t)  )
    end do


    if (itransim .eq. 1_ik .and. t .ne. 0_ik) then
        call hdf5_openf(datafile1, fileid)
        call hdf5_write(prodsim(1,1:nsim,1:tsim), fileid, 'prodsim_t')
        call hdf5_write(edusim(1,1:nsim,1:tsim), fileid, 'edusim_t')
        call hdf5_write(bsim(1,1:nsim,1:tsim), fileid, 'bsim_t')
        call hdf5_write(asim(1,1:nsim,1:tsim), fileid, 'asim_t')
        call hdf5_write(tuisim(1,1:nsim,1:tsim), fileid, 'tuisim_t')
        call hdf5_write(laborsim(1,1:nsim,1:tsim), fileid, 'laborsim_t')
        call hdf5_write(tau, fileid, 'tau')
        call hdf5_write(wl, fileid, 'wl')
        call hdf5_write(pi, fileid, 'pi')
        call hdf5_write(egridc, fileid, 'egridc')
        call hdf5_write(abilitycut, fileid, 'abilitycut_t')
        call hdf5_write(wealthcut, fileid, 'wealthcut_t')
        call hdf5_write(earnsim(1,1:nsim,1:tsim), fileid, 'earnsim_t')
        call hdf5_write(netphi, fileid, 'netphi_t')
        call hdf5_write(a0grid, fileid, 'a0grid')
        call hdf5_write(piea0_t, fileid, 'piea0_t')
        call hdf5_write(laborh(1,1,1:nabil,1,1:anum,1,1:tsim) , fileid, 'laborhini_t')
        call hdf5_write(bfchoice, fileid, 'bfchoice_t')
        call hdf5_write(abilitycutval, fileid, 'abilitycutval')
        call hdf5_write(earncutval_abil, fileid, 'earncutval_abil')
        call hdf5_write(wealthcutval, fileid, 'wealthcutval')
        call hdf5_write(earncutval_wealth, fileid, 'earncutval_wealth')
        call hdf5_write(avglaborini, fileid, 'avglaborini')
        call hdf5_write(avgearnini, fileid, 'avgearnini')
        call hdf5_write(earnpt, fileid, 'earnpt')
        call hdf5_write(earnabil, fileid, 'earnabil')
        call hdf5_write(yagg, fileid, 'yagg_t')
        call hdf5_write(gdp_data, fileid, 'gdp_data')
        call hdf5_write(drateearn, fileid, 'drateearn')
        call hdf5_write(ncost, fileid, 'ncost')
        call hdf5_write(chibnum, fileid, 'chibnum')
        call hdf5_write(abilsim(1,1:nsim,1:tsim), fileid, 'abilsim_t')
        call hdf5_closef(fileid)
    end if


    end subroutine transimulate

    subroutine bisectionwelfare( nabil,adnum, bdnum, conshd, laborhd, muh_ss, ftran, oval)
    integer(ik):: iter, adnum, bdnum,  j, ip, ie, ia, ib, iab, nabil
    real(rk):: oval,  conshd(jnum,nT,nabil,e_num,adnum,bdnum), laborhd(jnum,nT,nabil,e_num,adnum,bdnum), muh_ss(jnum,nT,nabil,e_num,adnum,bdnum), &
        utility(jnum,nT,nabil,e_num,adnum,bdnum)
    real(rk):: olow, ohigh, dist2, fval, flow, fhigh, ftran, fss

    intent(in):: nabil,adnum, bdnum, conshd, laborhd, muh_ss, ftran
    intent(out)::oval
    utility = 0.0_rk

    dist2 = 2.0_rk*precision
    olow = -0.9_rk; ohigh = 0.5_rk
    oval = olow;

    do j =1, jnum
        do ip = 1, nT
            do iab = 1, nabil
                do ie = 1, e_num
                    do ia = 1,adnum
                        do ib = 1,bdnum
                            if ( conshd(j,ip, iab,ie,ia,ib) .gt. 0.0_rk) then
                                if ( j .lt. jc) then
                                    utility(j,ip, iab,ie,ia,ib)=( ((1+oval)*conshd(j,ip, iab,ie,ia,ib))**(1-sigma))/(1-sigma) -psi_c*((laborhd(j,ip, iab,ie,ia,ib)**(1+eta))/(1+eta))
                                else                       !((consumption)**(1.0_rk - sigma))/(1.0_rk - sigma)-psi*((nval)**(1.0_rk+eta)/(1+eta)) + beta*evalf
                                    utility(j,ip, iab,ie,ia,ib)=( ((1+oval)*conshd(j,ip, iab,ie,ia,ib))**(1-sigma))/(1-sigma) -psi*((laborhd(j,ip, iab,ie,ia,ib)**(1+eta))/(1+eta))
                                end if
                            end if
                        end do

                    end do
                end do
            end do
        end do
    end do
    fss = sum(utility*muh_ss)


    fval =fss-ftran
    flow = fval;

    oval =ohigh;
    do j =1, jnum
        do ip = 1, nT
            do iab = 1, nabil
                do ie = 1, e_num
                    do ia = 1,adnum
                        do ib = 1,bdnum
                            if ( conshd(j,ip, iab,ie,ia,ib) .gt. 0.0_rk) then
                                if ( j .lt. jc) then
                                    utility(j,ip, iab,ie,ia,ib)=( ((1+oval)*conshd(j,ip, iab,ie,ia,ib))**(1-sigma))/(1-sigma) -psi_c*((laborhd(j,ip, iab,ie,ia,ib)**(1+eta))/(1+eta))
                                else
                                    utility(j,ip, iab,ie,ia,ib)=( ((1+oval)*conshd(j,ip, iab,ie,ia,ib))**(1-sigma))/(1-sigma) -psi*((laborhd(j,ip, iab,ie,ia,ib)**(1+eta))/(1+eta))
                                end if
                            end if
                        end do

                    end do
                end do
            end do
        end do
    end do

    fss = sum(utility*muh_ss)
    fval =fss-ftran
    fhigh = fval;

    if (flow*fhigh .gt. 0.0_rk) then
        write(*,'(1x,a)', advance = 'no') ' Could not bisect '
        write(*,*) flow, fhigh
    else
        iter =1_ik
        do while (dist2 .gt. precision)

            oval = (olow + ohigh)/2.0_rk;
            do j =1, jnum
                do ip = 1, nT
                    do iab =1, nabil
                        do ie = 1, e_num
                            do ia = 1,adnum
                                do ib = 1,bdnum
                                    if ( conshd(j,ip, iab,ie,ia,ib) .gt. 0.0_rk) then
                                        if ( j .lt. jc) then
                                            utility(j,ip,iab,ie,ia,ib)=( ((1+oval)*conshd(j,ip,iab,ie,ia,ib))**(1-sigma))/(1-sigma) -psi_c*((laborhd(j,ip,iab,ie,ia,ib)**(1+eta))/(1+eta))
                                        else
                                            utility(j,ip,iab,ie,ia,ib)=( ((1+oval)*conshd(j,ip,iab,ie,ia,ib))**(1-sigma))/(1-sigma) -psi*((laborhd(j,ip,iab,ie,ia,ib)**(1+eta))/(1+eta))
                                        end if
                                    end if

                                end do
                            end do

                        end do
                    end do
                end do
            end do
            fss = sum(utility*muh_ss)

            fval = fss-ftran

            if (fval*flow.gt.0) then
                olow = oval;
                flow = fval;
            else
                ohigh = oval;
                fhigh = fval;
            end if

            dist2 = ohigh - olow;
            iter = iter+1_ik
            if ( iter .gt. 100_ik) then
                write(*,*) 'dist2', dist2
            end if

        end do


    end if


    end subroutine bisectionwelfare


    subroutine welfare(nabil, tnum,  anum, bnum, adnum, bdnum, consh_t, laborh_t, conshss, laborhss, muh_ss, &
        ainiindex_c, ainiweight_c, aindex, aweight, bindex, bweight, welvalue)

    integer(ik):: nabil, anum, bnum, adnum, bdnum, j, ip, ie, ia, ib, tnum, t, t1, iab
    integer(ik):: ainiindex_c(nabil,adnum), aindex(adnum), bindex(bdnum)

    real(rk):: consh_t(jnum,nT,nabil,e_num,anum,bnum, tnum), laborh_t(jrnum,nT,nabil,e_num,anum,bnum, tnum), conshss(jnum,nT,nabil,e_num,anum,bnum), &
        laborhss(jrnum,nT,nabil,e_num,anum,bnum), muh_ss(jnum,nT,nabil, e_num,adnum,bdnum), ainiweight_c(nabil,adnum), &
        aweight(adnum), bweight(bdnum),  welvalue(tnum+jnum)
    real(rk):: valuet

    real(rk), allocatable:: conshdt(:,:,:,:,:,:), laborhdt(:,:,:,:,:,:), conshdss(:,:,:,:,:,:), laborhdss(:,:,:,:,:,:), &
        utility_t(:,:,:,:,:,:), conshpanel(:,:,:,:,:,:), laborpanel(:,:,:,:,:,:)

    intent(in):: nabil, anum, bnum, adnum, bdnum, consh_t, laborh_t, conshss, laborhss, muh_ss, &
        ainiindex_c, ainiweight_c, aindex, aweight, bindex, bweight,  tnum
    intent(out):: welvalue

    allocate( conshdss(jnum,nT,nabil, e_num,adnum,bdnum), laborhdss(jnum,nT,nabil, e_num,adnum,bdnum), conshdt(jnum,nT,nabil, e_num,adnum,bdnum), &
        laborhdt(jnum,nT,nabil, e_num,adnum,bdnum),  utility_t(jnum,nT,nabil, e_num,adnum,bdnum))

    call dist_welfare(conshss, laborhss, ainiindex_c, ainiweight_c, aindex, aweight, bindex, bweight,nabil, anum, bnum, adnum, bdnum, conshdss, laborhdss)

    allocate(conshpanel(jnum,nT,nabil, e_num,anum,bnum), laborpanel(jrnum,nT,nabil, e_num,anum,bnum))

    ! !$omp parallel do private(t, t1, conshpanel, laborpanel, conshdt, laborhdt, utility_t, valuet )
    do t =1, tnum+jnum !t =1 1913

        conshpanel =0.0_rk; laborpanel = 0.0_rk

        do j =1, jnum
            t1 = t-jnum+j-1
            t1 = min(t1, tnum)
            !write(*,*) 't1', j, t1
            if ( t1 .le. 0_ik) then  ! living in a steady-state world
                conshpanel(j,:,:,:,:,:) = conshss(j,:,:,:,:,:)
                if ( j .le. jrnum) then
                    laborpanel(j,:,:,:,:,:) = laborhss(j,:,:,:,:,:)
                end if
            else ! living in the transition world
                conshpanel(j,:,:,:,:,:) = consh_t(j,:,:,:,:,:,t1)
                if ( j .le. jrnum) then
                    laborpanel(j,:,:,:,:,:) = laborh_t(j,:,:,:,:,:,t1)
                end if
            end if

        end do

        call dist_welfare(conshpanel, laborpanel, ainiindex_c, ainiweight_c, aindex, aweight, bindex, bweight, nabil, &
            anum, bnum, adnum, bdnum, conshdt, laborhdt )

        utility_t = 0.0_rk
        do j =1, jnum
            do ip = 1, nT
                do iab = 1, nabil
                    do ie = 1, e_num
                        do ia = 1,adnum
                            do ib = 1,bdnum
                                if ( conshdt(j,ip, iab,ie,ia,ib) .gt. 0.0_rk) then
                                    if ( j .lt. jc) then
                                        utility_t(j,ip,iab,ie,ia,ib)=( conshdt(j,ip,iab,ie,ia,ib)**(1-sigma))/(1-sigma) -psi_c*(((laborhdt(j,ip,iab,ie,ia,ib))**(1+eta))/(1+eta))
                                    else
                                        utility_t(j,ip,iab,ie,ia,ib)=( conshdt(j,ip,iab,ie,ia,ib)**(1-sigma))/(1-sigma) -psi*(((laborhdt(j,ip,iab,ie,ia,ib))**(1+eta))/(1+eta))
                                    end if
                                end if

                            end do
                        end do
                    end do
                end do
            end do
        end do
        valuet = sum(utility_t*muh_ss)

        call bisectionwelfare( nabil, adnum, bdnum, conshdss, laborhdss, muh_ss, valuet, welvalue(t))
        write(*,*) 'bisectionwelfare done for year = ', t+1912, welvalue(t)

    end do


    ! !$omp end parallel do
    deallocate(conshpanel, laborpanel, conshdt, laborhdt, utility_t)

    end subroutine welfare


    subroutine dist_welfare(consh, laborh, ainiindex_c, ainiweight_c, aindex, aweight, bindex, bweight, nabil, anum, bnum, adnum, bdnum, conshd, laborhd )

    integer(ik):: j, ip, ie, ia, ib, iab, iloca, wa, ilocb, wb, anum, bnum, adnum, bdnum, nabil
    integer(ik):: ainiindex_c(nabil,adnum), aindex(adnum), bindex(bdnum)
    real(rk):: consh(jnum,nT,nabil,e_num,anum,bnum), laborh(jrnum,nT,nabil,e_num,anum,bnum), ainiweight_c(nabil,adnum), aweight(adnum), bweight(bdnum)
    real(rk):: conshd(jnum,nT,nabil,e_num,adnum,bdnum), laborhd(jnum,nT,nabil,e_num,adnum,bdnum)

    intent(in):: consh, laborh, ainiindex_c, ainiweight_c, aindex, aweight, bindex, bweight, nabil, anum, bnum, adnum, bdnum
    intent(out):: conshd, laborhd
    conshd=0.0_rk; laborhd=0.0_rk

    do j = 1,jnum
        do ip = 1, nT
            do iab = 1, nabil
                do ie = 1, e_num
                    do ia = 1, adnum
                        if (j.lt.jc) then
                            iloca = ainiindex_c(iab, ia);  wa = ainiweight_c(iab, ia);
                        else
                            iloca  = aindex(ia); wa = aweight(ia);
                        end if

                        do ib = 1, bdnum
                            ilocb = bindex(ib); wb = bweight(ib)

                            conshd(j, ip,iab,ie, ia, ib) = wb*(wa*consh(j,ip,iab,ie,iloca,ilocb) + (1.0_rk-wa)*consh(j,ip,ie,iab,iloca+1,ilocb))+(1.0-wb)*(wa*consh(j,ip,iab,ie,iloca,ilocb+1) + (1.0_rk-wa)*consh(j,ip,iab,ie,iloca+1,ilocb+1))
                            if ( j .le. jrnum) then
                                laborhd(j, ip,iab, ie, ia, ib) = wb*(wa*laborh(j,ip,iab,ie,iloca,ilocb) + (1.0_rk-wa)*laborh(j,ip,iab,ie,iloca+1,ilocb))+(1.0-wb)*(wa*laborh(j,ip,iab,ie,iloca,ilocb+1) + (1.0_rk-wa)*laborh(j,ip,iab,ie,iloca+1,ilocb+1))
                            end if
                        end do
                    end do
                end do
            end do
        end do

    end do


    end subroutine dist_welfare

    end module transition_solve
