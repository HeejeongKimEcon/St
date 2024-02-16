    module estimation_mod


    use Lib_kindset
    use Lib_conshock
    use Lib_inequality
    use Lib_sort
    use Lib_grid
    use omp_lib
    use StdDebtCrisis_mod
    use get_params
    implicit none


    public:: educost, ss_solve_est, sortme


    contains


    subroutine  educost(nabil, abilgrid, pabil, anum,  nparam,  x, egridc, pia0, piec0, ashock,  a0grid, psycost, probcost, netphi, piea0 )
    integer(ik)::  nabil, anum,ie, ia, nparam,   ic, j
    real(rk)::     psi1, psi2, psi3, phi0, phi1, phi2, egridc(e_num), a0grid(nabil, anum),  aval, eval, x(nparam), psyshock(ncost), &
        psycost(nabil, anum,ncost), costbound, netphi(nabil, anum), probcost(nabil, anum, ncost), piec0(e_num), pia0(anum), chidage(jnum), &
        ashock(anum), piea0(nabil, anum), pic(ncost,ncost), pic0(ncost), phic0, phic1, phic2, phic3, phic4, chid1, chid2, chid3, chid4, chid5,&
        abilgrid(nabil), pabil(nabil)
    intent(in):: nabil, abilgrid, pabil, anum, egridc, ashock, x, nparam,  pia0, piec0
    intent(out):: psycost, probcost, netphi, a0grid, piea0

    psi1 = x(1)
    psi2 = x(2)
    psi3 = x(3)
    phi0 =(x(4))
    phi1 =(x(5))
    phi2 =(x(6))
    phic0 = (x(7))
    phic1 = (x(8))
    phic2 = (x(9))

    if (inormal .eq. 1_ik) then
        phic3 = (x(10))
        phic4 = 0.0_rk
    else !when uniform distribution. 
        phic3 = (x(10)) !we use the value here. 
        phic4 = 0.0_rk
    end if

    do ie= 1, nabil
        eval = exp(abilgrid(ie))
        do ia = 1,anum
            a0grid(ie,ia) = psi1 +psi3*exp(ashock(ia)) + psi2*(eval)
            piea0(ie,ia) = pia0(ia)*pabil(ie)
        end do
    end do
    !a0grid = exp(a0grid)

    !monetary education cost
    netphi = phi0;
    do ie = 1,nabil
        eval = exp(abilgrid(ie))
        do ia = 1,anum
            aval = a0grid(ie, ia)
            netphi(ie, ia) =  (phi0)+(phi1)*(aval)+ (phi2)*(eval)
        end do
    end do

    if ( ncost .gt. 1_ik) then
        if ( inormal .eq. 1_ik) then
            call tauchen (0.0_rk, 1.0_rk, 0.0_rk,  multiple, ncost, psyshock, pic)  ! standard normal
            call ergodicdist(ncost, pic, precerg, pic0)
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

    end subroutine educost


    subroutine ss_solve_est(nabil, abilgrid, pabil, eagridj, chibshock, anum, bnum, adnum, bdnum,ebar,  pop,  weights, targ, ntarg, chir,psycost, probcost, netphi, pi, &
        tbar, chidage,   t, nsim, r_b,  wp, wl, wh, lambda, agrid, adgrid, bgrid, bdgrid, egridc,  egridl,piec, piel, a0grid, &
        piel0, piec0, pia0, piea0,  error,chibprob)

    integer(ik) :: nabil,  anum, bnum, adnum, bdnum,  count, ntarg, iss

    integer(ik), allocatable:: bindex(:), pchoicec(:,:,:,:,:,:), pchoicenc(:,:,:,:,:,:),  aindex(:),ad0index(:,:), &
        ainiindex(:,:),payindex(:,:),nopayindex(:), idxprod(:,:), paysim(:,:), idxability(:), idxchib(:)

    integer(ik) ::ie, ia,  ip,   iter, t,  nsim, i, j

    real(rk) ::  r_b,  wp, wl, wh,  kaggh, laggh, debtout,  &
        kaggl, laggl, kagg, lagg, yagg, stddebtholder, lambdaval, rval, wlval, whval, educated, noneducated, totalstddebt, &
        avgdebt, chid,  tbar, navg,  phicost, kaggd,laggd, debtoutc, error, &
        debtoutnc, laggh_jc, pi,partransfer_c, partransfer_nc, partransfer, debtbyage(5), chir, totaldebtbyage(5), ddebtpercent,&
        colrate(3,4), totalborrower, colrate_abil(3), colrate_wealth(4), abilratio, wealthratio, ebar, ddebtbyage(5), ddebtmubyage(5), ddebtpercentmu, debtbyagesim(5)


    real(rk):: pop(jnum),  lambda(nT), agrid(anum), bgrid(bnum), adgrid(adnum), bdgrid(bdnum), egridc(e_num),  egridl(e_num), a0grid(nabil, anum), &
        pia0(anum), netphi(nabil,anum), psycost(nabil,anum,ncost), probcost(nabil,anum, ncost), piec(e_num,e_num),  piel(e_num,e_num), piel0(e_num), piec0(e_num), &
        prob(ncost), muh(jnum,nT,nabil,e_num,adnum,bdnum), mul(jnum,nT,nabil,e_num,adnum,bdnum), mud(jnum,nT,nabil,e_num,adnum,bdnum), fracborrower(5),&
        model(ntarg), targ(ntarg), weights(ntarg), piea0(nabil,anum), ddebt(jnum), chidage(jnum), afsim(jnum,nsim), bfsim(jnum,nsim), &
        edufsim(nsim), eagridj(jnum,nabil,e_num,2), pabil(nabil), abilgrid(nabil),chibprob(chibnum), chibshock(chibnum)


    real(rk), allocatable :: gl(:,:,:,:,:,:), vl(:,:,:,:,:,:),  bhini(:,:,:),  gh(:,:,:,:,:,:), vhini(:,:,:),  vh(:,:,:,:,:,:), edudecision(:,:,:,:), bfchoice(:,:,:,:), &
        avglabord(:), ad0weight(:,:), aweight(:), ainiweight(:,:), bweight(:),  nopayweight(:),payweight(:,:), laborl(:,:,:,:,:,:), &
        laborh(:,:,:,:,:,:),avglaborl(:),avglaborh(:), lclaborsim(:), lcbsim(:), lcasim(:), edulclaborsim(:,:), edulcbsim(:,:), edulcasim(:,:), &
        defaultsim(:,:), lcearnsim(:), edulcearnsim(:,:),  asim(:,:), bsim(:,:), &
        edusim(:,:), prodsim(:,:), consh(:,:,:,:,:,:), abilsim(:), vhcollege(:,:,:,:,:)

    intent(in):: nabil, abilgrid, pabil, anum, bnum, adnum, bdnum,  pop, weights, targ, ntarg,  chir,psycost, probcost, netphi, &
        pi, ebar, tbar, chidage,     t,  nsim,  r_b,  wp, wl, wh, lambda, agrid, adgrid, &
        bgrid, bdgrid, egridc,  egridl,piec, piel, a0grid, piel0, piec0, pia0, piea0, eagridj,chibprob, chibshock

    intent(out)::  error



    iter = 0_ik; phicost = 0.0_rk
    do ie  = 1, nabil
        do ia = 1, anum
            phicost = phicost + netphi(ie,ia)*pia0(ia)*pabil(ie)
        end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Decision Rules !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!

    allocate(vl(jnum,nT,nabil,e_num,anum,bnum),gl(jnum,nT,nabil,e_num,anum,bnum),  vh(jnum,nT,nabil,e_num,anum,bnum), gh(jnum,nT,nabil,e_num,anum,bnum),  & !college-educated retirees & workers
        laborl(jrnum,nT,nabil,e_num,anum,bnum), laborh(jrnum,nT,nabil,e_num,anum,bnum), vhini(nabil,anum,chibnum), bhini(nabil,anum,chibnum), bindex(bdnum), &
        edudecision(nabil,anum, ncost, chibnum), bfchoice(nabil,anum, ncost, chibnum), pchoicec(jpaynum,nT,nabil, e_num, anum, bnum), pchoicenc(jpaynum,nT,nabil, e_num, anum, bnum), &
        aweight(adnum), ainiweight(nabil, adnum),aindex(adnum), ainiindex(nabil, adnum), ad0index(nabil, anum),   ad0weight(nabil, anum),  bweight(bdnum), &
        nopayindex(bnum), nopayweight(bnum), payindex(bnum,nT), payweight(bnum,nT), consh(jnum,nT,nabil,e_num,anum,bnum), vhcollege(jc-1_ik,nabil,anum,bnum,chibnum))

    ! calculate the index/weights in advance to save time
    call indexcal(nabil, anum, bnum, adnum, bdnum, r_b, lambda,  bgrid, bdgrid,adgrid, a0grid, agrid, aindex, ainiindex, aweight, ainiweight, &
        bindex, bweight, payindex, payweight, nopayindex, nopayweight, ad0index, ad0weight)

        
        
    ! solve decision rules
    call decisionrule(tbar, nabil, eagridj, chibshock, anum, bnum, ebar, chir, psycost, pi, piel0, piec0, netphi, bbar, &
        a0grid,  t,  vh, vl, vhcollege, chidage, lambda, agrid, bgrid,  wh, wl,  piec, piel,  &
        payindex, payweight, nopayindex, nopayweight, gh, gl, vh, vl, &
        pchoicec, pchoicenc, laborh, laborl, edudecision, bfchoice, bhini, vhini, consh, vhcollege)


        
        
    ! distribution
    call entiredist(chibprob, nabil, eagridj, ebar, anum, bnum, adnum, bdnum,wp, pi, probcost,  pop, ad0index, ad0weight,  adgrid, edudecision, piel0, piec0, &
        bfchoice, bdgrid, laborh, t,  r_b, lambda, pchoicec, pchoicenc, ainiindex, ainiweight, gh, gl,  aindex, aweight, &
        bindex, bweight, piec, piel,   laborl, muh,  mul,  mud, kagg, lagg, yagg, &
        debtout,navg, educated, noneducated, stddebtholder, partransfer_c, partransfer_nc, partransfer, &
        debtbyage, ddebt, fracborrower, avgdebt,totaldebtbyage, ddebtpercent, totalborrower, ddebtbyage, ddebtmubyage, ddebtpercentmu)

    totaldebtbyage = - totaldebtbyage/yagg

    allocate(asim(jnum,nsim), bsim(jnum,nsim), idxprod(jnum,nsim), edusim(jnum,nsim), prodsim(jnum,nsim), paysim(jnum,nsim), &
        idxability(nsim), abilsim(nsim) , idxchib(nsim))
    iss = 0_ik

    call simulate(chibshock, chibprob, nabil, eagridj, pabil, abilgrid, ebar, anum, bnum,  yagg,  probcost, pi, wh, wl, psycost,  &
        nsim,  egridc, egridl, piec, piel, piec0, piel0, piea0, a0grid, edudecision, bfchoice, gl, gh, laborl, laborh, pchoicec, pchoicenc, &
        agrid, bgrid,  r_b, lambda, netphi, asim, bsim, idxprod, edusim, prodsim, paysim, idxability, abilsim, model,ntarg, iss, idxchib)


    error = 0.0_rk
    do i = 1, ntarg
        error = error + weights(i)*dabs(model(i) - targ(i))**2.0_rk
    end do

    write(*,'(1x,a)') ' data and model moments'
    
    write(*,'(1x,a, 2f8.4)') ' fraction of college educated       ', targ(8), model(8)!educated!model(8)
    write(*,'(1x,a, 2f8.4)') ' average net education cost         ', targ(9), model(9)!phicost/yagg!model(9)
    write(*,'(1x,a, 2f8.4)') ' avg parental transfers             ', targ(10), model(10)!partransfer/yagg!model(10)
    write(*,'(1x,a, 2f8.4)') ' avg cumulative std debt            ', targ(11), model(11)!-avgdebt/yagg!model(11)
    write(*,'(1x,a, 2f8.4)') ' graduating seniors with std debt   ', targ(12), model(12)!stddebtholder!model(12)
    write(*,'(1x,a, 2f8.4)') ' aggregate std debt                 ', targ(13), model(13)!-debtout/yagg!model(13)

    end subroutine ss_solve_est

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                           !
    !   Sort a set of n function values, f(n)   !
    !                                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Taken from downhillsimplex program

    subroutine sortme(n, f, in0)

    integer(ik):: n, i, i0, i1, in0(n), fail
    real(rk):: f(n), f0, f1

    do i = 1, n
        in0(i) = i
    end do

    fail = 1

    do

        if (fail.eq.0) then
            return
        end if

        fail = 0

        do i = 1, n-1

            i0 = in0(i)
            i1 = in0(i+1)
            f0 = f(i0)
            f1 = f(i1)

            if (f1.lt.f0) then
                in0(i) = i1
                in0(i+1) = i0
                fail = 1
            end if

        end do

    end do

    end subroutine sortme





    end module estimation_mod