    module steadystate_solve


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


    subroutine ss_solve(tbar, chibprob, nabil, chibshock, anum, bnum, adnum, bdnum, probcost,eagridj,pabil, abilgrid, pop, &
        chir, psycost, datafile, olgo, netphi,  pi,  phid_abil, chidage,    t,  nsim, &
        r_b,  wp, wl, wh,  lambda, agrid, adgrid, bgrid, bdgrid, egridc, egridl,&
        piec,  piel, a0grid, piel0, piec0, pia0, piea0,  mul, mud, muh,muhcollege, &
        asim, bsim, idxprod, edusim, prodsim, paysim, idxability, abilsim,  epgridc, idxchib)

    integer(ik) :: np, olgo, ie, ia, anum, bnum, adnum, bdnum, ntarg, iss,  t,  nsim, nabil

    integer(ik), allocatable:: bindex(:), pchoicec(:,:,:,:,:,:), pchoicenc(:,:,:,:,:,:),  aindex(:),ad0index(:,:), &
        ainiindex(:,:),payindex(:,:),nopayindex(:),  idxprod(:,:),  paysim(:,:), idxability(:), idxchib(:)

    real(rk) ::   r_b,  wp, wl, wh,  debtout, kagg, lagg, yagg, stddebtholder, educated, noneducated,    &
        avgdebt, avgdebt_c, avgdebt_d,  navg, phicost,  ddebtpercentmu, &
        pi,partransfer_c, partransfer_nc, partransfer, debtbyage(5), chir, totaldebtbyage(5), ddebtpercent,  &
        totalborrower, tbar, colrate, drateh, drated

    real(rk):: pop(jnum),  lambda(nT), agrid(anum), bgrid(bnum), adgrid(adnum), bdgrid(bdnum), egridc(e_num), egridl(e_num), &
        a0grid(nabil, anum),  pia0(anum), netphi(nabil,anum), psycost(nabil, anum, ncost), probcost(nabil, anum, ncost),  piec(e_num,e_num), &
        piel(e_num,e_num), piel0(e_num), piec0(e_num),  ddebt(jnum), ddebtbyage(5), ddebtmubyage(5), &
        fracborrower(5), piea0(nabil, anum), chidage(jnum), epgridc(epnum), eagridj(jnum,nabil, e_num, 2), pabil(nabil), abilgrid(nabil),    &
        chibshock(chibnum), chibprob(chibnum),muhcollege(jc-1, nabil, adnum, bdnum, chibnum), phid_abil(nabil)


    real(rk), allocatable :: gl(:,:,:,:,:,:), vl(:,:,:,:,:,:),  bhini(:,:,:),  gh(:,:,:,:,:,:), vhini(:,:,:),  vh(:,:,:,:,:,:), edudecision(:,:,:,:),&
        bfchoice(:,:,:,:), para(:), avglabord(:), ad0weight(:,:), aweight(:), ainiweight(:,:), &
        bweight(:),  nopayweight(:),payweight(:,:), laborl(:,:,:,:,:,:), laborh(:,:,:,:,:,:),avglaborl(:),avglaborh(:), &
        muh(:,:,:,:,:,:), mul(:,:,:,:,:,:), mud(:,:,:,:,:,:), asim(:,:), bsim(:,:),  &
        edusim(:,:), prodsim(:,:), model(:), abilsim(:,:), &
        vhcollege(:,:,:,:,:), laborhcollege(:,:,:,:,:), ghcollege(:,:,:,:,:)


    integer(hid_t):: fileid
    character(30):: datafile,  datafile1, datafile2, datafile3, datestring, timestring

    intent(in):: tbar, nabil, pop,  anum, bnum, adnum, bdnum, pabil, abilgrid,&
        chir, netphi,   pi,   chidage,  t,  nsim,  r_b,  wp, wl, wh, lambda, agrid, adgrid, bgrid, &
        bdgrid, egridc,  egridl, piec,  piel, a0grid, piel0, piec0, pia0, datafile, olgo, psycost, piea0, &
        probcost, chibshock, epgridc, eagridj, chibprob, phid_abil
    intent(out)::  muhcollege, mul, mud, muh, asim, bsim, idxprod, edusim, prodsim, paysim,  idxability, abilsim, idxchib

    iss = 1_ik
    ntarg = 13_ik
    !!!!!!!!!!!!!!!!!!!!!!!!
    !   reporting values   !
    !!!!!!!!!!!!!!!!!!!!!!!!
    datafile1 = datafile(1:olgo+4)
    datafile1(olgo+5:olgo+12) = '_mom.txt'

    datafile2 = datafile(1:olgo+4)
    datafile2(olgo+5:olgo+7) = '.h5'

    datafile3 = datafile(1:olgo+4)
    datafile3(olgo+5:olgo+12) = '_sim.h5'

    phicost = 0.0_rk
    do ie  = 1, nabil
        do ia = 1, anum
            phicost = phicost + netphi(ie,ia)*pia0(ia)*pabil(ie)
        end do
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Decision Rules !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!

    allocate(vl(jnum,nT,nabil, e_num,anum,bnum),gl(jnum,nT,nabil,e_num,anum,bnum),  vh(jnum,nT,nabil, e_num,anum,bnum), gh(jnum,nT,nabil,e_num,anum,bnum),  & !college-educated retirees & workers
        laborl(jrnum,nT,nabil,e_num,anum,bnum), laborh(jrnum,nT,nabil,e_num,anum,bnum), vhini(nabil,anum,chibnum), bhini(nabil,anum,chibnum), bindex(bdnum), &
        edudecision(nabil,anum,ncost,chibnum), bfchoice(nabil,anum,ncost,chibnum), pchoicec(jpaynum,nT,nabil, e_num, anum, bnum), pchoicenc(jpaynum,nT,nabil, e_num, anum, bnum), &
        aweight(adnum), ainiweight(nabil,adnum),aindex(adnum), ainiindex(nabil,adnum), ad0index(nabil,anum), ad0weight(nabil, anum), bweight(bdnum), &
        nopayindex(bnum), nopayweight(bnum), payindex(bnum,nT), payweight(bnum,nT), avglaborl(jnum), avglabord(jnum), avglaborh(jnum),&
        vhcollege(jc-1_ik,nabil,anum,bnum,chibnum), ghcollege(jc-1_ik,nabil,anum,bnum,chibnum), laborhcollege(jc-1_ik,nabil,anum,bnum,chibnum))

    ! calculate the index/weights in advance to save time
    call indexcal(nabil, anum, bnum, adnum, bdnum,  r_b, lambda,  bgrid, bdgrid, adgrid, a0grid, agrid, aindex, ainiindex, aweight, ainiweight, &
        bindex, bweight, payindex, payweight, nopayindex, nopayweight, ad0index, ad0weight)
        
    ! solve decision rules
    call decisionrule(tbar, nabil, eagridj, chibshock, anum, bnum, chir, psycost,  piel0, piec0,   netphi, bbar, &
        a0grid,  t, vh, vl, vhcollege, chidage, lambda, agrid, bgrid,  wh, wl, piec, piel, phid_abil, &
        payindex, payweight, nopayindex, nopayweight, gh, gl, vh, vl, &
        pchoicec, pchoicenc, laborh, laborl, edudecision, bfchoice, bhini, vhini,  vhcollege, laborhcollege, ghcollege)
    write(*,*) 'just finished decisionrule'


    allocate(muh(jnum,nT,nabil,e_num,adnum,bdnum), mul(jnum,nT,nabil,e_num,adnum,bdnum), mud(jnum,nT,nabil,e_num,adnum,bdnum))
    ! distribution
    call entiredist(pabil, phid_abil, chibprob, nabil, eagridj, anum, bnum, adnum, bdnum, wp, probcost,  pop, ad0index, ad0weight,              &
        adgrid, edudecision,  piel0, piec0,  bfchoice, bdgrid, laborhcollege, laborh, t,  r_b, lambda, pchoicec, pchoicenc,             &
        ainiindex, ainiweight, ghcollege, gh, gl, aindex, aweight, bindex, bweight, piec, piel,  laborl, muhcollege, muh,  mul,  mud,               &
        kagg, lagg, yagg,  debtout, navg, educated, noneducated, stddebtholder, partransfer_c, partransfer_nc, partransfer,             & 
        debtbyage, ddebt, fracborrower, avgdebt, avgdebt_c, avgdebt_d,totaldebtbyage, ddebtpercent, totalborrower,                      &
        ddebtbyage, ddebtmubyage, ddebtpercentmu, colrate, drateh, drated)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                   !
    !                targets and moments                !
    !                                                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*) '!                                                   !'
    write(*,*) '!      Targets Identifying 4 Parameters             !'
    write(*,*) '!                                                   !'
    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)
    write(*,*) '                                                   targets                                    data           model'
    write(*,*)
    write(*,'(1x, a, f8.4)')    ' A1.  K/Y                                                              :    3.00   ',    kagg/yagg
    write(*,*)
    write(*,'(1x, a, f8.4)')    ' A2.  average labor supply                                             :    0.30   ',    navg
    write(*,*)
    write(*,'(1x, a, f7.3)')    ' A3.  Frisch elasticity of male labor supply                           :    0.50   ',    (1/eta)
    write(*,*)
    write(*,'(1x, a, f8.4)')    ' A4.  Std loan borrowing limit(unsubs) $23,000/31123                   :    0.74   ',    -bbar/yagg
    write(*,*)
    write(*,'(1x, a, f8.4)')    ' A5.  Graduation rate (NLSY79 cohort/ 25-85 pop)                       :           ',    educated
    write(*,*)
    write(*,'(1x, a, f8.4)')    ' A6.  Graduation rate (NLSY79 cohort/ j=1:jc-1)                        :    0.1887 ',    colrate
    write(*,*)    
    write(*,'(1x, a, f8.4)')    ' A7.  students graduating with student debt (NLSY97)                   :    0.60   ',    stddebtholder
    write(*,*)
    write(*,'(1x, a, f8.4)')    ' A8.  Average student debt at graduation($5,229/31123)                 :    0.1680 ',    -avgdebt/yagg
    write(*,*)
    write(*,'(1x, a, f8.4)')    ' A9.  outstanding std debt as a fraction of gdp (51.33B/7339.58B)      :    0.0070	',    -debtout/yagg
    write(*,*)
    write(*,'(1x, a, f8.4)')    ' A10. fraciton of delinquent student debt ( out of total debt)         :    _____  ',    ddebtpercent
    write(*,*)
    write(*,'(1x, a, f8.4)')    ' A11.  Average net college tuition(w/o selection 5631/31123)           :    0.1809 ',    phicost/yagg
    write(*,*)
    write(*,'(1x, a, f8.4)')    ' Average parental transfer: $4,719.56 (0.1516)                         : ', partransfer/yagg
    write(*,'(1x, a, f8.4)')    ' Average parental transfer for college: $9,722.83 (0.3124)             : ', partransfer_c/yagg
    write(*,'(1x, a, f8.4)')    ' Average parental transfer for non-college: $3,230.14 (0.1038)         : ', partransfer_nc/yagg
    write(*,*)
    write(*,'(1x, a, f8.4)')    ' Average debt at graduation for all           : $5,229.176 (0.1680)    : ', -avgdebt/yagg
    write(*,'(1x, a, f8.4)')    ' Average debt at graduation for college       : $9,493.583 (0.3050)    : ', -avgdebt_c/yagg
    write(*,'(1x, a, f8.4)')    ' Average debt at graduation for dropout       : $2,618.589 (0.0841)    : ', -avgdebt_d/yagg    
    write(*,*)
    write(*,'(1x, a, f8.4)')    ' delinqunecy rate for all workers (SCF1998)      : 0.1863              : ',  ddebtpercentmu
    write(*,'(1x, a, f8.4)')    ' delinqunecy rate for college graduate workers   : 0.1388              : ',  drateh
    write(*,'(1x, a, f8.4)')    ' delinqunecy rate for college dropouts workers   : 0.2523              : ',  drated        
    write(*,*)
    write(*,'(1x, a, f8.4)')    ' total amount of defaulted debt as a fraction of gdp                   : ', sum(ddebt)/yagg
    write(*,'(1x, a, f8.4)')    ' Total fraction of borrowers from whole popluation                     : ', fracborrower    
    write(*,*)
    write(*,*)
    !write(*,'(1x,(a,f8.2),a)')    'debt at age under 30 (Data: 0.36, Model:', debtbyage(1)/sum(debtbyage(1:5))  , ')'
    !write(*,'(1x,(a,f8.2),a)')    'debt at age 30-39    (Data: 0.33, Model:', debtbyage(2)/sum(debtbyage(1:5))  , ')'
    !write(*,'(1x,(a,f8.2),a)')    'debt at age 40-49    (Data: 0.16, Model:', debtbyage(3)/sum(debtbyage(1:5))  , ')'
    !write(*,'(1x,(a,f8.2),a)')    'debt at age 50-59    (Data: 0.11, Model:', debtbyage(4)/sum(debtbyage(1:5))  , ')'
    !write(*,'(1x,(a,f8.2),a)')    'debt at age 60+      (Data: 0.04, Model:', debtbyage(5)/sum(debtbyage(1:5))  , ')'
    !write(*,*)
    !write(*,*)
    write(*,'(1x, a, f12.6)') 'population hitting blow', sum(muhcollege(:,:,:,1,:))+sum(muh(:,:,:,:,:,1))+sum(mul(:,:,:,:,:,1))+sum(mud(:,:,:,:,:,1))
    write(*,'(1x, a, f12.6)') 'population hitting ahigh', sum(muhcollege(:,:,adnum,:,:))+sum(muh(:,:,:,:,adnum,:))+sum(mul(:,:,:,:,adnum,:))+sum(mud(:,:,:,:,adnum,:))
    write(*,*) 
    



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !       SIMULATION CODE             !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ( isim .eq. 1_ik) then
        allocate(asim(jnum,nsim), bsim(jnum,nsim), idxprod(jnum,nsim), edusim(jnum,nsim), prodsim(jnum,nsim), paysim(jnum,nsim), &
            idxability(nsim), abilsim(jnum, nsim), model(ntarg), idxchib(nsim))
        call simulate(chibshock, chibprob, nabil, eagridj, pabil, abilgrid, phid_abil, anum, bnum,  yagg,  probcost, pi, wh, wl, psycost,  &
             nsim,  egridc, piec, piel, piec0, piel0, piea0, a0grid, edudecision, bfchoice, gl, gh, ghcollege, laborl, laborh, laborhcollege, pchoicec, pchoicenc, &
            agrid, bgrid,  r_b, lambda, netphi, asim, bsim, idxprod, edusim, prodsim, paysim,idxability, abilsim, iss,  &
            idxchib)

    end if

    if (iarray .eq. 1_ik) then
        call hdf5_openf(datafile2, fileid)

        np = 18;
        allocate(para(np))
        para(1) = wh
        para(2) = wl
        para(3) = jnum
        para(4) = e_num
        para(5) = anum
        para(6) = bnum
        para(7) = kagg
        para(8) = lagg
        para(9) = adnum
        para(10) = bdnum
        para(11) = nT
        para(12) = psi
        para(13) = eta
        para(14) = epnum
        para(15) = etnum
        para(16) = nabil
        para(17) = ncost
        para(18) = bbar
        
        call hdf5_write(para, fileid, 'para')
        call hdf5_write(egridc, fileid, 'egridc')
        call hdf5_write(egridl, fileid, 'egridl')
        call hdf5_write(lambda, fileid, 'lambda')
        !call hdf5_write(gl, fileid, 'gl')
        !call hdf5_write(laborl, fileid, 'laborl')
        call hdf5_write(vl, fileid, 'vl')
        !call hdf5_write(gh, fileid, 'gh')
        call hdf5_write(vh, fileid, 'vh')
        call hdf5_write(vhini, fileid, 'vhini')
        call hdf5_write(bhini, fileid, 'bhini')
        call hdf5_write(edudecision, fileid, 'edudecision')
        call hdf5_write(bfchoice, fileid, 'bfchoice')
        !call hdf5_write(muh, fileid, 'muh')
        !call hdf5_write(mul, fileid, 'mul')
        call hdf5_write(bgrid, fileid, 'bgrid')
        call hdf5_write(a0grid, fileid, 'a0grid')
        call hdf5_write(pia0, fileid, 'pia0')
        !call hdf5_write(pchoicec, fileid, 'pchoicec')
        !call hdf5_write(pchoicenc, fileid, 'pchoicenc')
        call hdf5_write(laborh, fileid, 'laborh')
        !call hdf5_write(laborl, fileid, 'laborl')
        call hdf5_write(epgridc, fileid, 'epgridc')
        call hdf5_write(abilgrid, fileid, 'abilgrid')
        call hdf5_write(eagridj, fileid, 'eagridj')
        call hdf5_write(pabil, fileid, 'pabil')
        call hdf5_write(chibshock, fileid, 'chibshock')
        call hdf5_write(chibprob, fileid, 'chibprob')
        call hdf5_write(psycost, fileid, 'psycost')
        call hdf5_write(probcost, fileid, 'probcost')

        
        call hdf5_closef(fileid)

    end if

    call date_and_time(datestring, timestring)
    write (*,*) ' time: ', timestring(1:2), ':', timestring(3:4), '  date: ', datestring(7:8), '/', &
        datestring(5:6), '/', datestring(1:2), datestring(3:4)
    write(*,*)

    if (itext .eq. 1_ik) then
        open(unit = 50, file = datafile1, status = 'unknown', action = 'write')
        write(50,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(50,*)'!                                                   !'
        write(50,*) '!      Targets Identifying 4 Parameters             !'
        write(50,*) '!                                                   !'
        write(50,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(50,*)
        write(50,'(1x, a, f8.4)')    ' A1.  K/Y                                                              :    3.00   ',    kagg/yagg
        write(50,*)
        write(50,'(1x, a, f8.4)')    ' A2.  average labor supply                                             :    0.30   ',    navg
        write(50,*)
        write(50,'(1x, a, f7.3)')    ' A3.  Frisch elasticity of male labor supply                           :    0.50   ',    (1/eta)
        write(50,*)
        write(50,'(1x, a, f8.4)')    ' A4.  Std loan borrowing limit(unsubs) $23,000/31123                   :    0.74   ',    -bbar/yagg
        write(50,*)
        write(50,'(1x, a, f8.4)')    ' A5.  Graduation rate (NLSY79 cohort/ 25-85 pop)                       :           ',    educated
        write(50,*)
        write(50,'(1x, a, f8.4)')    ' A6.  Graduation rate (NLSY79 cohort/ j=1:jc-1)                        :    0.1887 ',    colrate
        write(50,*)    
        write(50,'(1x, a, f8.4)')    ' A7.  students graduating with student debt (NLSY97)                   :    0.60   ',    stddebtholder
        write(50,*)
        write(50,'(1x, a, f8.4)')    ' A8.  Average student debt at graduation($5,229/31123)                 :    0.1680 ',    -avgdebt/yagg
        write(50,*)
        write(50,'(1x, a, f8.4)')    ' A9.  outstanding std debt as a fraction of gdp (51.33B/7339.58B)      :    0.0070	',    -debtout/yagg
        write(50,*)
        write(50,'(1x, a, f8.4)')    ' A10. fraciton of delinquent student debt ( out of total debt)         :    _____  ',    ddebtpercent
        write(50,*)
        write(50,'(1x, a, f8.4)')    ' A11.  Average net college tuition(w/o selection 5631/31123)           :    0.1809 ',    phicost/yagg
        write(50,*)
        write(50,'(1x, a, f8.4)')    ' Average parental transfer: $4,719.56 (0.1516)                         : ', partransfer/yagg
        write(50,'(1x, a, f8.4)')    ' Average parental transfer for college: $9,722.83 (0.3124)             : ', partransfer_c/yagg
        write(50,'(1x, a, f8.4)')    ' Average parental transfer for non-college: $3,230.14 (0.1038)         : ', partransfer_nc/yagg
        write(50,*)
        write(50,'(1x, a, f8.4)')    ' Average debt at graduation for all           : $5,229.176 (0.1680)    : ', -avgdebt/yagg
        write(50,'(1x, a, f8.4)')    ' Average debt at graduation for college       : $9,493.583 (0.3050)    : ', -avgdebt_c/yagg
        write(50,'(1x, a, f8.4)')    ' Average debt at graduation for dropout       : $2,618.589 (0.0841)    : ', -avgdebt_d/yagg    
        write(50,*)
        write(50,'(1x, a, f8.4)')    ' delinqunecy rate for all workers (SCF1998)      : 0.1863              : ',  ddebtpercentmu
        write(50,'(1x, a, f8.4)')    ' delinqunecy rate for college graduate workers   : 0.1388              : ',  drateh
        write(50,'(1x, a, f8.4)')    ' delinqunecy rate for college dropouts workers   : 0.2523              : ',  drated      
        write(50,*)
        write(50,'(1x, a, f8.4)')    ' total amount of defaulted debt as a fraction of gdp                   : ', sum(ddebt)/yagg
        write(50,'(1x, a, f8.4)')    ' Total fraction of borrowers from whole popluation                     : ', fracborrower    
        write(50,*)
        close(50)
    end if



    deallocate(vl,gl,  vh, gh,  laborl, vhini, bhini, bindex, edudecision, bfchoice, pchoicec, pchoicenc, &
        aweight, ainiweight,aindex, ainiindex, ad0index, ad0weight, bweight, nopayindex, nopayweight, payindex, &
        payweight, avglaborl, avglabord, avglaborh)


    end subroutine ss_Solve

    end module steadystate_solve