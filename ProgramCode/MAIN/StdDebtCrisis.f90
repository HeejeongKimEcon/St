    program StdDebtCrisis

    use Lib_kindset
    use Lib_conshock
    use Lib_inequality
    use Lib_sort
    use Lib_grid
    use StdDebtCrisis_mod
    use omp_lib
    use steadystate_solve
    use transition_solve
    use get_params

    !To optimize programs to run faster, set the following
    !Configuation Properties\Fortran\Run-time\set Check Stack Frame to no
    !Fortran\General\Set Debug information format to none
    !Fortran\ Optimization to Maximize Speed
    !Use Intel Math Kernel Library in Properties\Fortran\Libraries
    !In properties\linker\manifest\generate manifest file 'NO'

    implicit none

    integer(ik) :: anum, bnum, adnum, bdnum, olgo, j, i, t, nsim, ia, ie,  ic, ib, ief, bmnum, nabil, iab,ncut,k,k2
    integer(ik), allocatable:: idxprod(:,:), paysim(:,:), idxability(:), idxchib(:), ptcut(:), abilcut(:)

    real(rk) ::   r_b, wp, wl, wh,  lb, ub, tbar, pi, chir, chid1, chid2, chid3, chid4, chid5,  &
        psi1, psi2, eval, phi0, phi1, phi2, aval, psi3, phic0, phic1, phic2, phic3, phic4, psi4, chib, &
        grah, gral, abillow, abilhigh, chibl

    real(rk):: starttime, endtime

    real(rk), allocatable :: lambda(:), agrid(:), adgrid(:), bgrid(:), bdgrid(:), bgridtemp(:), bdgridtemp(:),   &
        egridc(:), epgridc(:), etgridc(:), egridl(:), piec(:,:),  piel(:,:), a0grid(:,:), piel0(:), piec0(:),  pia0(:),  &
        piea0(:,:), muh(:,:,:,:,:,:), expel(:), expeh(:), netphi(:,:), mul(:,:,:,:,:,:), mud(:,:,:,:,:,:), psycost(:,:,:), probcost(:,:,:), pia(:,:), &
        ashock(:), pop(:), psyshock(:), pic(:,:), pic0(:), chidage(:), asim(:,:), bsim(:,:),  edusim(:,:), prodsim(:,:), &
        csim(:,:), laborsim(:,:),  piea0_t(:,:,:), abilgrid(:), pabil(:), &
        eagridj(:,:,:,:), eagridj_t(:,:,:,:,:), abilsim(:,:), chibshock(:), chibprob(:), chibprobt(:,:), &
        netphi_pt(:), netphi_abil(:), ptcutamt(:), abilcutamt(:), countpt(:), countabil(:), muhcollege(:,:,:,:,:), phid_abil(:)


    real(rk), allocatable:: piec_t(:,:,:), piel_t(:,:,:), piel0_t(:,:), piec0_t(:,:),phi_t(:,:,:),  wh_t(:), &
        psycost_t(:,:,:,:), probcost_t(:,:,:,:), bbar_t(:)
    character(30):: datafile, paramfile, datestring, timestring


    write(*,*) ' Student Debt Crisis '
    write(*,*)

    call date_and_time(datestring, timestring)

    write (*,*) ' time: ', timestring(1:2), ':', timestring(3:4), '  date: ', datestring(7:8), '/', &
        datestring(5:6), '/', datestring(1:2), datestring(3:4)
    write(*,*)

    call cpu_time ( starttime )

    datafile = 'KK_stddebt_'
    olgo = len_trim(datafile)
    datafile(olgo + 1:olgo+4) = timestring(1:4)

    paramfile = ' '
    paramfile(1:olgo+4) = datafile(1:olgo+4)
    paramfile(olgo+5:olgo+13) = '_para.txt'

    if (idefault .eq. 1_ik) then
        write(*,*) ' Default is forbidden at all (i.e. idefault = 1_ik)'
    end if

    !!!!!!!!!!!!!!!!!!!!!
    !!    Parameters   !!
    !!!!!!!!!!!!!!!!!!!!!

    ! Calibrate to 1984 economy
    r_b   = r+0.031_rk                          ! interest rate on student debt (91-day t-bill + 3.1% in 1997)
    wp    = 1.107058003         !1.240553629_rk! 1.415522407 !1.360336471!                           ! Wage premium for 1984 (NLSY)
    pi    = 0.25!1.0!95!1.8!7! 2.46504960294672    !2.88!1.95 ! adjusting egridc at age j = 1 to 4.
    chir  = 0.0                                 ! disutility of holding negative debt after retirement.
    psi1  = -3.20!0.03!-3.5!-2.9!-2.9!-2.60813402343018             ! Parental transfer paramete
    psi2  = 0.22!0.005!0.015!0.007                                ! a0grid(ie,ia) = psi1 + psi2log(eval) + psi3*ashock(ia)
    psi3  = 0.23!0.043!8!6!0.26!0.15                             !
    psi4  = 0.0                                 !
    phi0  = -1.21!-1.57!-2.7               ! Education cost parameter
    phi1  = 0.60!0.12              ! netphi(ie, ia) =  (phi0)+(phi1)*(aval)+ (phi2)*(eval)
    phi2  = 0.70!
    if (inormal .eq. 1_ik) then
        phic0 = 5.0!2.0!              ! Psychic cost parameter
        phic1 = 24.8!3.0!                  ! phic3+phic0*psyshock(ic)-(phic1)*(aval)-(phic2)*(eval) -(phic4)*(aval)*(eval)
        phic2 = 35.5!11.5!16.1!1.0!
        phic3 = 18.0_rk!0.4!16.0!19
        phic4 = 0.0
    else
        phic0 = 10.0_rk
        phic1 = 0.0
        phic2 = 0.0!1.85
        phic3 = 0.0_rk!14.7_rk
    end if


    wl = 1.0_rk                          ! wage for non-college worker (normalized)
    wh = wp*wl                           ! wage for college grad
    tbar = 0.25!                         ! study time (abbott et al. 2019)

    chib =0.245 !0.249_rk               ! disutility of borrowing student debt during college education(must not exceed 0.25). 
    chibl = 0.240_rk !0.221
    allocate(chibshock(chibnum), chibprob(chibnum), chibprobt(chibnum,chibnum))
    call tauchen(chib, dsqrt(1.0_rk) , 0.0_rk,  3.0_rk, chibnum, chibshock, chibprobt)
    chibprob = chibprobt(1,:)
    chibshock = exp(chibshock)



    call linspace(chibl, chib, chibnum, chibshock)!(0.234_rk
    !chibshock(1:chibnum) = chib
    chibprob(1:chibnum) = 1.0_rk/chibnum

    write(*,*) 'chibshock = ', chibshock

    chid1 = 0.003_rk!0.07!0.05_rk                     ! disutility from being delinquent
    chid2 = chid1!0.012!0.00
    chid3 = chid1!0.012!0.00
    chid4 = chid1!0.012!0.00
    chid5 = chid1!0.012!0.00
    grah = 0.81765398_rk!0.9724652_rk
    gral = 0.64816225_rk!0.5658639_rk

    allocate(chidage(jnum))
    do j = 1,jnum
        if (j.lt.12_ik) then
            chidage(j) = chid1
        elseif (j.ge.12 .and. j.lt.22) then
            chidage(j) = chid2
        elseif (j.ge.22 .and. j.lt.32) then
            chidage(j) = chid3
        elseif (j.ge.32 .and. j.lt.42) then
            chidage(j) = chid4
        else
            chidage(j) = chid5
        end if
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !calculate the probability mass for each age
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(pop(jnum))
    pop(1:jnum) = 1.0_rk/real(jnum)


    !Payment schedule
    allocate(lambda(nT))
    do i = 1, nT
        lambda(i) = (1.0_rk)/(nT + 1_rk - i)
    end do

    !Experience premium
    allocate(expel(jrnum), expeh(jrnum))
    do j = 1, jrnum
        !expeh(j) = 0.0517288*(j) - 0.0009269*((j)**2)     ! College
        expeh(j) = 0.0189489*(j) - 0.0006032*((j)**2)     ! College        
    end do
    expeh = dexp(expeh)

    do j = 1, jrnum
        !expel(j) = 0.0352376*(j) -0.0006287*((j)**2)     ! No-College
        expel(j) = 0.0189489*(j) - 0.0006032*((j)**2)     ! No-College         
    end do
    expel = dexp(expel)

    !Asset grid & debt grid(20,30,20,20) !bnum should be at least 20.
    anum = 20_ik; adnum = 20_ik; bnum = 20_ik; bdnum = 20_ik; bmnum = 10_ik
    allocate(agrid(anum), bgrid(bnum), bgridtemp(bnum), adgrid(adnum), bdgrid(bdnum), bdgridtemp(bnum))
    lb = al - 1.0*al +1.0
    ub = ah - 1.0*al +1.0
    call logspace(lb, ub, anum, agrid)
    call logspace(lb, ub, adnum, adgrid)
    agrid(1:anum) = agrid(1:anum) + (al-1.0)
    adgrid(1:adnum) = adgrid(1:adnum) + (al-1.0)

    !call linspace(blow, bbar, bmnum, bgridtemp)
    !bgrid(1:bmnum) = bgridtemp(1:bmnum)
    !
    !call linspace(blow, bbar, bmnum, bdgridtemp)
    !bdgrid(1:bmnum) = bdgridtemp(1:bmnum)
    !
    !deallocate(bgridtemp, bdgridtemp)
    !allocate(bgridtemp(11), bdgridtemp(11))
    !
    !call linspace(bbar, bhigh, 11, bgridtemp)
    !bgrid(bmnum+1:bnum) = bgridtemp(2:11)
    !
    !call linspace(bbar, bhigh,11, bdgridtemp)
    !bdgrid(bmnum+1:bdnum) = bdgridtemp(2:11)

    lb = bhigh - 1.0*bhigh +1.0
    ub = -blow - 1.0*bhigh +1.0
    call logspace(lb, ub, bnum, bgridtemp)
    call logspace(lb, ub, bdnum, bdgridtemp)

    bgridtemp(1:bnum) = bgridtemp(1:bnum) + (bhigh-1.0)
    bdgridtemp(1:bdnum) = bdgridtemp(1:bdnum) + (bhigh-1.0)

    do i = 1, bnum
        bgrid(i) = -bgridtemp(bnum+1-i)
    end do
    do i = 1, bdnum
        bdgrid(i) = - bdgridtemp(bdnum+1-i)
    end do


    !ability grid
    OPEN(26, FILE='ability_10.txt', STATUS='OLD', ACTION='READ')
    read(26,*) nabil
    read(26,*) abillow
    read(26,*) abilhigh

    allocate(abilgrid(nabil), pabil(nabil), phid_abil(nabil))

    do ia= 1, nabil
        read(26,*) abilgrid(ia)
    end do
    do ia= 1, nabil
        read(26,*) pabil(ia)
    end do
    close(26)

    !dropout over initial ability
    OPEN(27, FILE='dropout_abil.txt', STATUS='OLD', ACTION='READ')

    do ia= 1, nabil
        read(27,*) phid_abil(ia)
    end do
    close(27)
    !write(*,*) phid_abil
    !phid_abil(:) = 0.3_rk
    !Labor productivity shock process for college-educated & non-college-educated households
    !!!!!Earning shocks (1984)

    allocate(egridc(e_num),epgridc(epnum),etgridc(etnum), egridl(e_num), piec(e_num,e_num), piel(e_num,e_num),&
        piel0(e_num), piec0(e_num), eagridj(jnum, nabil, e_num, 2))

    call steadyshock(pi, nabil, grah, gral, abilgrid, expeh, expel, egridc, epgridc, etgridc, egridl, piec, piel, piel0, piec0, eagridj)

    !ebar= pi*dot_product(egridc, piec0)*expeh(1)


    !Distribution for parental transfer
    allocate(ashock(anum), a0grid(nabil, anum),  pia0(anum), pia(anum, anum), piea0(nabil, anum))

    call tauchen (meana0, stda, 0.0_rk, multiple, anum, ashock, pia)
    call ergodicdist(anum, pia, precerg, pia0)
    do ie= 1, nabil
        !eval = exp(abilgrid(ie))
        eval = abilgrid(ie)
        do ia = 1,anum
            !a0grid(ie,ia) = psi1 +psi3*exp(ashock(ia))+ psi2*(eval)
            a0grid(ie,ia) = psi1 +psi3*ashock(ia)+ psi2*(eval)
            piea0(ie,ia) = pia0(ia)*pabil(ie)
        end do
    end do
    a0grid = exp(a0grid)


    !education cost
    allocate(netphi(nabil,anum))
    netphi = phi0;
    do ie = 1,nabil
        !eval = exp(abilgrid(ie))
        eval =abilgrid(ie)
        do ia = 1,anum
            aval = log(a0grid(ie, ia))
            netphi(ie, ia) =  (phi0)+(phi1)*(aval)+ (phi2)*(eval)
            !write(*,*) ie, ia, aval, netphi(ie,ia)
        end do
    end do
    netphi = exp(netphi)


    allocate(netphi_abil(nabil), netphi_pt(anum), ptcut(anum), abilcut(nabil), ptcutamt(4), abilcutamt(3), countpt(4), countabil(3))
    ptcutamt = 0.0_rk; abilcutamt = 0.0_rk; countpt = 0.0_rk; countabil = 0.0_rk

    do iab = 1,nabil
        netphi_abil(iab) = dot_product(netphi(iab,1:anum),piea0(iab,1:anum))
        netphi_abil(iab) = netphi_abil(iab)/sum(piea0(iab,1:anum))
    end do

    do ia = 1,anum
        netphi_pt(ia) = dot_product(netphi(1:nabil,ia),piea0(1:nabil,ia))
        netphi_pt(ia) = netphi_pt(ia)/sum(piea0(1:nabil,ia))
    end do

    ncut = int(anum/4);
    do k =1,4
        if (k.eq. 4_ik) then
            k2 = anum
        else
            k2 = k*ncut
        end if
        do i = (k-1)*ncut+1, k2
            ptcut(i)= k
            ptcutamt(k) = ptcutamt(k) + netphi_pt(i)
            countpt(k) = countpt(k) + 1.0_rk
        end do

    end do

    ncut = int(nabil/3);
    do k =1,3
        if (k.eq. 3_ik) then
            k2 = nabil
        else
            k2 = k*ncut
        end if
        do i = (k-1)*ncut+1, k2
            abilcut(i)= k
            abilcutamt(k) = abilcutamt(k) + netphi_abil(i)
            countabil(k) = countabil(k) + 1.0_rk
        end do
    end do
     
    ptcutamt = ptcutamt/countpt
    abilcutamt = abilcutamt/countabil

    
    write(*,*)
    write(*,*)
    write(*,*) 'netphi over ability '
    write(*,*) (abilcutamt/0.4084)*gdp_data
    write(*,*) 'netphi abil ratio = ', abilcutamt(3)/abilcutamt(1)
    write(*,*)
    write(*,*)
    write(*,*) 'netphi over parental transfer '
    write(*,*) (ptcutamt/0.4084)*gdp_data
    write(*,*) 'netphi PT ratio = ', ptcutamt(4)/ptcutamt(1)

      
    
    !psychic education cost
    allocate(psycost(nabil, anum, ncost), probcost(nabil, anum, ncost), psyshock(ncost), pic(ncost, ncost), pic0(ncost))

    if ( ncost .gt. 1_ik) then
        if ( inormal .eq. 1_ik) then
            call tauchen (0.0_rk, 1.0_rk, 0.0_rk,  multiple, ncost, psyshock, pic)  ! standard normal
            call ergodicdist(ncost, pic, precerg, pic0)
            !            psyshock = exp(psyshock)
        else
            call linspace(0.0_rk, phic0, ncost, psyshock)
            pic0(1:ncost) = 1.0_rk/ncost
        end if
    else
        psyshock(ncost) = 0.0_rk
        pic0 = 1.0_rk
        pic = 1.0_rk
    end if


    do ie = 1, nabil
        eval = exp(abilgrid(ie))
        do ia = 1,anum
            aval = a0grid(ie,ia)
            do ic = 1,ncost
                if ( ncost .gt. 1_ik) then
                    if (inormal.eq.1_ik) then
                        psycost(ie,ia,ic)= phic3+phic0*psyshock(ic)-(phic1)*(aval)-(phic2)*(eval) -(phic4)*(aval)*(eval)
                    else
                        psycost(ie,ia,ic)= phic3+psyshock(ic)-(phic1)*(aval)-(phic2)*(eval) -(phic4)*(aval)*(eval)
                    end if
                else
                    psycost(ie,ia,ic)= phic3+psyshock(ic)-(phic1)*(aval)-(phic2)*(eval) -(phic4)*(aval)*(eval)
                end if
                probcost(ie,ia,ic) =pic0(ic)*pabil(ie)*pia0(ia)
            end do
        end do
    end do


    

    !psycost = psycost - 33.3;

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!     solve steady state here
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    t=0_ik
    if (isim .eq. 1_ik) then
        nsim = 10000   ! number of simulation
    else
        nsim = 1
    end if
    allocate(asim(jnum, nsim), bsim(jnum, nsim), idxprod(jnum, nsim), edusim(jnum, nsim), prodsim(jnum, nsim), paysim(jnum, nsim), &
        csim(jnum, nsim), laborsim(jnum, nsim), idxability(nsim), abilsim(jnum,nsim), idxchib(nsim))
    allocate(muh(jnum,nT,nabil,e_num,adnum,bdnum), mul(jnum,nT,nabil,e_num,adnum,bdnum), mud(jnum,nT,nabil,e_num,adnum,bdnum),  &
        muhcollege(jc-1, nabil, adnum, bdnum, chibnum))

    call ss_solve(tbar, chibprob, nabil, chibshock, anum, bnum, adnum, bdnum, probcost, eagridj,pabil, abilgrid,  &
        pop,  chir,psycost, datafile, olgo,netphi, pi, phid_abil,   chidage,  t, nsim, r_b,&
        wp, wl, wh, lambda, agrid, adgrid, bgrid, bdgrid, egridc,  egridl,piec, piel, a0grid, &
        piel0, piec0, pia0, piea0,  mul, mud, muh, muhcollege, asim, bsim, idxprod, edusim, prodsim, paysim,  &
        idxability, abilsim,  epgridc, idxchib)


    if (itext .eq. 1_ik) then
        open(unit = 50, file = paramfile, status = 'unknown', action = 'write')
        write(50,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(50,*)'!               Parameter Values                    !'
        write(50,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(50,'(1x, a, f10.4)') 'r = ', r
        write(50,'(1x, a, f10.4)') 'wl = ', wl
        write(50,'(1x, a, f10.4)') 'wh = ', wh
        write(50,'(1x, a, i2.0)')    'jnum   :  ', jnum
        write(50,'(1x, a, i2.0)')    'jrnum   :  ', jrnum
        write(50,'(1x, a, i1.0)')    'jc   :  ', jc
        write(50,'(1x, a, i2.0)')    'jpaynum   :  ', jpaynum
        write(50,'(1x, a, f8.3)')    'beta   :  ', beta
        write(50,'(1x, a, f8.3)')    'r_b   :  ', r_b
        write(50,'(1x, a, i2.0)')    'nT   :  ', nT
        write(50,'(1x, a, f8.3)')    'psi   :  ', psi
        write(50,'(1x, a, f8.3)')    'psi_c   :  ', psi_c
        write(50,'(1x, a, f8.3)')    'pi   :  ', pi
        write(50,'(1x, a, f8.4)')    'chir   :  ', chir
        write(50,'(1x, a, f8.4)')    'psi1(parental transfer)   :  ', psi1
        write(50,'(1x, a, f8.4)')                       'psi2   :  ', psi2
        write(50,'(1x, a, f8.4)')                       'psi3   :  ', psi3
        write(50,'(1x, a, f8.4)')    'phi0(edu_cost)   :  ', phi0
        write(50,'(1x, a, f8.4)')              'phi1   :  ', phi1
        write(50,'(1x, a, f8.4)')              'phi2   :  ', phi2
        write(50,'(1x, a, f8.4)')    'phic0(psychic_edu_cost)   :  ', phic0
        write(50,'(1x, a, f8.4)')                      'phic1   :  ', phic1
        write(50,'(1x, a, f8.4)')                      'phic2   :  ', phic2
        write(50,'(1x, a, f8.4)')                      'phic3   :  ', phic3
        write(50,'(1x, a, f8.4)')                      'phic4   :  ', phic4
        write(50,'(1x, a, f8.4)')      'bbar   :  ', bbar
        write(50,'(1x, a, f8.4)')      'chibl   :  ', chibl
        write(50,'(1x, a, f8.4)')      'chib   :  ', chib        
        close(50)
    end if

    !!!exporting parameters to use them in estimation

    open(unit = 33, file = 'ss_par.txt', status = 'unknown', action = 'write')
    write(33,*) r_b
    write(33,*) wp
    write(33,*) wl
    write(33,*) wh
    write(33,*) tbar
    write(33,*) pi
    write(33,*) chir
    write(33,*) psi1
    write(33,*) psi2
    write(33,*) psi3
    write(33,*) phi0
    write(33,*) phi1
    write(33,*) phi2
    write(33,*) phic0
    write(33,*) phic1
    write(33,*) phic2
    write(33,*) phic3
    write(33,*) phic4
    write(33,*) chid1
    write(33,*) chid2
    write(33,*) chid3
    write(33,*) chid4
    write(33,*) chid5
    write(33,*) anum
    write(33,*) bnum
    write(33,*) adnum
    write(33,*) bdnum
    write(33,*) chib
    write(33,*) nabil

    do j = 1,jnum
        write(33,*) chidage(j)
    end do
    do j = 1,jnum
        write(33,*) pop(j)
    end do
    do i = 1,nT
        write(33,*) lambda(i)
    end do
    do j = 1,jrnum
        write(33,*) expel(j)
    end do
    do j = 1,jrnum
        write(33,*) expeh(j)
    end do
    do ia = 1,anum
        write(33,*) agrid(ia)
    end do
    do ib = 1,bnum
        write(33,*) bgrid(ib)
    end do
    do ia = 1,adnum
        write(33,*) adgrid(ia)
    end do
    do ib = 1,bdnum
        write(33,*) bdgrid(ib)
    end do
    do ia= 1, anum
        write(33,*) ashock(ia)
    end do
    do ia= 1, anum
        write(33,*) pia0(ia)
    end do
    do ie = 1, e_num
        write(33,*) egridc(ie)
    end do
    do ie = 1, e_num
        write(33,*) egridl(ie)
    end do
    do ie = 1, e_num
        do ief = 1, e_num
            write(33,*) piec(ie, ief)
        end do
    end do
    do ie = 1, e_num
        do ief = 1, e_num
            write(33,*) piel(ie, ief)
        end do
    end do
    do ie = 1, e_num
        write(33,*) piel0(ie)
    end do
    do ie = 1, e_num
        write(33,*) piec0(ie)
    end do
    do j = 1, jnum
        do ia = 1, nabil
            do ie = 1, e_num
                do ief = 1, 2
                    write(33,*) eagridj(j,ia,ie,ief)
                end do
            end do
        end do
    end do
    do ie = 1, nabil
        write(33,*) abilgrid(ie)
    end do
    do ie = 1, nabil
        write(33,*) pabil(ie)
    end do

    do ie = 1, chibnum
        write(33,*) chibshock(ie)    
    end do
    do ie = 1, chibnum
        write(33,*) chibprob(ie)
    end do

    close(33)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!     solve transitional dynamics
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (itran.eq.1_ik) then
        !Transitional dynamics period
        allocate(piec_t(e_num,e_num,tnum), piel_t(e_num,e_num,tnum), piel0_t(e_num,tnum), piec0_t(e_num,tnum),phi_t(nabil, anum,tnum),  wh_t(tnum), &
            psycost_t(nabil,anum,ncost,tnum), probcost_t(nabil, anum, ncost, tnum), bbar_t(tnum), piea0_t(nabil, anum,tnum),                        &
            eagridj_t(jnum, nabil, e_num, 2, tnum))
        call transtioninput( nabil, pabil, anum,  wh, netphi, probcost, psycost,  piec, piel, epgridc, etgridc, eagridj, &
            psycost_t, probcost_t, wh_t, phi_t, bbar_t, piec_t, piel_t, piel0_t, piec0_t, pia0, pic0, piea0_t, eagridj_t)

        !initial distribution estimation
        call transition(tbar, chibprob, eagridj_t, abilgrid, pabil, phid_abil, chibshock, nabil,anum, bnum, adnum, bdnum, wp,  pop, chir,     &
            datafile, olgo,r_b, chidage,  pi,  psycost_t, probcost_t, wh_t, phi_t, bbar_t, piec_t, piel_t, piel0_t, piec0_t, asim, bsim, idxprod,   &
            edusim, prodsim,paysim,  idxability, abilsim, lambda, pia0, a0grid, agrid, bgrid, adgrid, bdgrid, egridl, egridc,                       &
            muhcollege, muh, mul, mud, piea0_t,nsim, idxchib)

    end if

    ! Report time
    call cpu_time ( endtime )
    write(*,*)
    write(*,'(1x,a,f8.2,a)') ' elapsed time: ', endtime - starttime, ' seconds '


    end program stddebtcrisis