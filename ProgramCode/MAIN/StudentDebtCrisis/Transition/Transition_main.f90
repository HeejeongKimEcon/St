
    program Transition_main

    use Lib_kindset
    use Lib_conshock
    use Lib_inequality
    use Lib_sort
    use Lib_grid
    use StdDebtCrisis_mod
    use omp_lib
    use transition_solve
    use get_params
    !To optimize programs to run faster, set the following
    !Configuation Properties\Fortran\Run-time\set Check Stack Frame to no
    !Fortran\General\Set Debug information format to none
    !Fortran\ Optimization to Maximize Speed
    !Use Intel Math Kernel Library in Properties\Fortran\Libraries
    !In properties\linker\manifest\generate manifest file 'NO'

    implicit none

    integer(ik) :: anum, bnum, adnum, bdnum, jc,  olgo, j, i, t, nsim, ia, ie,  ic, ib, ief,  bmid
        
    real(rk) ::   r_b, wp, wl, wh,  lb, ub, tbar, addparental, pi, chir, ebar, chid1, chid2, chid3, chid4, chid5,  &
        psi1, psi2, eval, phi0, phi1, phi2, aval, psi3, phic0, phic1, phic2, phic3, phic4, psi4
        
    real(rk):: starttime, endtime

    real(rk), allocatable :: lambda(:), agrid(:), adgrid(:), bgrid(:), bdgrid(:),  &
        egridc(:), epgridc(:), etgridc(:), egridl(:),epgridl(:),etgridl(:), piec(:,:), piepc(:,:), piea0(:,:), &
        pietc(:,:), piel(:,:),  piepl(:,:), pietl(:,:), a0grid(:,:), piel0(:), piec0(:),  pia0(:), muh(:,:,:,:,:), &
        expel(:), expeh(:), netphi(:,:), mul(:,:,:,:,:), mud(:,:,:,:,:), psycost(:,:,:), probcost(:,:,:), pia(:,:), &
        ashock(:), pop(:), psyshock(:), pic(:,:), pic0(:), chidage(:)
    
    real(rk), allocatable:: piec_t(:,:,:), piel_t(:,:,:), piel0_t(:,:), piec0_t(:,:),phi_t(:,:,:),  wh_t(:), &
            psycost_t(:,:,:,:), probcost_t(:,:,:,:), bbar_t(:)
    character(30):: datafile, datestring, timestring


    
    !!!!!Read parameters from steady state text file!!!!!
    OPEN(33, FILE='ss_par.txt', STATUS='OLD', ACTION='READ')  
    read(33,*) jc
    read(33,*) r_b
    read(33,*) wp
    read(33,*) wl
    read(33,*) wh
    read(33,*) tbar
    read(33,*) addparental
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
    allocate(pop(jnum), lambda(nT), expel(jnum), expeh(jnum), chidage(jnum))

    do j = 1,jnum
        read(33,*) chidage(j)
    end do    
    do j = 1,jnum
        read(33,*) pop(j)
    end do
    do i = 1,nT
        read(33,*) lambda(i)
    end do
    do j = 1,jnum
        read(33,*) expel(j)
    end do    
    do j = 1,jnum
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
        ashock(anum), a0grid(e_num, anum), pia0(anum),  piea0(e_num,anum))
    do ia= 1, anum
        read(33,*) ashock(ia)
    end do    
    do ia= 1,anum
        read(33,*) pia0(ia)
    end do  
    do ie = 1, e_num
        do ia = 1, anum
            read(33,*) a0grid(ie, ia)
        end do
    end do
     do ie = 1, e_num
        do ia = 1, anum
            read(33,*) piea0(ie, ia)
        end do
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
    
    close(33)
    
  !Distribution for parental transfer
  !  allocate(ashock(anum), a0grid(e_num, anum),  pia0(anum), pia(anum, anum), piea0(e_num, anum))
    allocate(pia(anum, anum))
    psi1 = -3.47834396362305!-3.50379566793351
    psi2 = 0.174190846085548! 0.174965991216521
    psi3 = 0.355025386810303! 0.302480217850425 
    psi4 = 0.0
    call tauchen (meana0, stda, 0.0_rk, 3.0_rk, anum, ashock, pia)
    call  ergodicdist(anum, pia, precerg, pia0)
    do ie= 1, e_num
        eval = egridc(ie)
        do ia = 1,anum
            a0grid(ie,ia) = psi1 +psi3*ashock(ia)+ psi2*log(eval/varp1) 
            piea0(ie,ia) = pia0(ia)*piec0(ie)
        end do
    end do
    a0grid = exp(a0grid)


    !education cost
    phi0 = 1.826754003018141E-002 !0.766988366952555E-002 
    phi1 = 1.519245207309723E-002 ! 1.264728636508787E-002 
    phi2 = 1.766294144093990E-002 ! 2.511777583806133E-002 
    allocate(netphi(e_num,anum))
    netphi = phi0;
    do ie = 1,e_num
        eval = egridc(ie)
        do ia = 1,anum
            aval = a0grid(ie, ia)
            netphi(ie, ia) =  (phi0)+(phi1)*(aval)+ (phi2)*(eval)
        end do
    end do
    write(*,*) ' while borrowing limit per year is ', bbar*(gdp_data/(4.0_rk))

    
    !psychic education cost
    allocate(psycost(e_num, anum, ncost), probcost(e_num, anum, ncost), psyshock(ncost), pic(ncost, ncost), pic0(ncost))
    if (inormal .eq. 1_ik) then 
        phic0 = 20.0254516601562! 20.0!1.01244996462876 !1.00
        phic1 = 1.02545166015625! 1.0!2.83789357975926 !1.00
        phic2 = 5.02545166015625! 5.0!1.02797565664496 !1.00
        phic3 = 16.5254516601562! 10.0!15.0
        phic4 = 0.0
    else 
        phic0 = 50.0_rk
        phic1 = 1.0
        phic2 = 11.0
        phic3 = 0.0
    end if

    if ( inormal .eq. 1_ik) then
        call tauchen (0.0_rk, 1.0_rk, 0.0_rk, multiple, ncost, psyshock, pic)  ! standard normal 
        call ergodicdist(ncost, pic, precerg, pic0)
    else
        call linspace(-phic0, phic0, ncost, psyshock)
        pic0(1:ncost) = 1.0_rk/ncost
    end if

      
    do ie = 1, e_num
        eval = egridc(ie)
        do ia = 1,anum
            aval = a0grid(ie,ia)
            do ic = 1,ncost
                psycost(ie,ia,ic)= phic3+phic0*psyshock(ic)-(phic1)*(aval)-(phic2)*(eval) -(phic4)*(aval)*(eval)
                probcost(ie,ia,ic) =pic0(ic)*piec0(ie)*pia0(ia)
            end do
        end do
    end do
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!     solve transitional dynamics
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (itran.eq.1_ik) then
        
        !Transitional dynamics period
        allocate(piec_t(e_num,e_num,tnum), piel_t(e_num,e_num,tnum), piel0_t(e_num,tnum), piec0_t(e_num,tnum),phi_t(e_num, anum,tnum),  wh_t(tnum), &
            psycost_t(e_num,anum,ncost,tnum), probcost_t(e_num,anum,ncost,tnum), bbar_t(tnum))
        call transtioninput( anum,  wh, netphi, probcost, psycost,  piec, piel, epgridc, etgridc,  &
            psycost_t, probcost_t, wh_t, phi_t, bbar_t, piec_t, piel_t, piel0_t, piec0_t)
        
        !initial distribution estimation
        call transition(anum, bnum, adnum, bdnum, ebar, wp,  pop, chir, datafile, olgo,  jc, r_b, chidage,tbar,  pi, addparental,  &
            psycost_t, probcost_t, wh_t, phi_t, bbar_t, piec_t, piel_t, piel0_t, piec0_t, &
            lambda, expel, expeh,  pia0, a0grid, agrid, bgrid, adgrid, bdgrid, egridc, egridl,muh, mul, mud )

    end if

    ! Report time
    call cpu_time ( endtime )
    write(*,*)
    write(*,'(1x,a,f8.2,a)') ' elapsed time: ', endtime - starttime, ' seconds '



    end program Transition_main

