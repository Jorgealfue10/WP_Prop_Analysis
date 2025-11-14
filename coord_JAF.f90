program psitrJAF 
    
    ! kinds
    use decimal,     only: long, sip, dop

    ! datos/arrays de la función de onda y rejilla
    use psidef
    use griddatmod
    use dvrdatmod

    ! interfaces de lectura / dvr
    use rdfilesmod
    use rddvrmod

    ! varios de Quantics (si ya los tenías, mantenlos)
    use datenmod
    use compdatmod
    use global
    use globmemmod
    use aglobal
    use memdimmod
    use keyunitsmod
    use compdatmod
    use rdopermod
    use operdef
    use psidef
    use iofile
    use iopsidef
    use datenmod
    use compdatmod 
    use compdat1mod
    use iodvrdef
    use iorst, only: rdrstpsi
    use channels, only: irst
    use logdat
    use gwplib, only: param_from_psi
    use channels
    use logdat
    use maxv
    use psioutmod

    !==============================!
    !    Variables declarations    !
    !==============================!

    implicit none
    ! --- indexes and counters ---
    integer(long) :: i, j, k, ilbl
    integer(long) :: ip, ig, ith, idx
    integer(long) :: ierr
    integer(long) :: chkdvr, chkgrd, chkpsi, chkdat
    integer(long) :: check, check_dvr
    integer(long) :: nrp, nrg, nth
    integer(long) :: tot_dim
    integer(long) :: ndof1
    integer(long), allocatable :: fdvr(:)

    ! --- arguments ---
    integer(long) :: ntime

    ! --- real arrays ---
    real(dop), allocatable, target :: newgrid(:,:)      ! grid completo transformado
    real(dop), allocatable, target :: auxort(:)

    ! --- pointers ---
    real(dop), pointer :: rp_grid(:), rg_grid(:), th_grid(:)
    real(dop), pointer :: rpp_grid(:), rgp_grid(:), thp_grid(:)

    ! --- real scalar variables ---
    real(dop) :: rp, rg, theta                  ! coordenadas iniciales
    real(dop) :: rpp, rgp, thetap               ! coordenadas transformadas
    real(dop) :: r12, r13, r23                  ! coordenadas internas
    real(dop) :: m1, m2, m3                    ! masas atómicas
    real(dop) :: time

    ! --- namefiles ---
    character(len=200) :: filename

    ! --- DVR-related variables ---
    real(dop), allocatable :: ort(:)            ! grid points
    real(dop), allocatable :: trafo(:,:)        ! transformation matrix
    real(dop), allocatable :: dvrmat(:,:)       ! DVR matrix
    real(dop), allocatable :: fftp(:), kinsph(:)
    complex(dop), allocatable :: hin(:), rueck(:)
    complex(dop), allocatable :: exphin(:), exprueck(:)
    integer(long), allocatable :: jsph(:), msph(:), fftfak(:)

    ! --- arrays para lectura de psi ---
    integer(long), allocatable :: jindx(:), trajst(:), gwpdep(:)
    complex(dop), allocatable :: psi(:)
    complex(sip), allocatable :: spsi(:), spsigrd(:)
    complex(dop), allocatable :: adgwp(:)
    complex(dop), allocatable :: psigrd(:), workc(:)
    logical :: lrst,lerr

    !==============!
    !  Tiempo WF   !
    !==============!
    write(*,fmt="(a)",advance="no") "WF Time: "
    read(*,*) ntime

    !==============================!
    !   Variables initialization   !
    !==============================!

    chkdvr = 1
    chkgrd = 1
    chkpsi = 1
    chkdat = 1
    lrddvr = .true.
    lerr = .true.
    lrst=.false.

    pi = dacos(-1.d0)

    call default
    call adefault

    allocate(fdvr(maxdim))

    !==============================!
    ! 1) Leer memdim del OPER      !
    !==============================!
    open(ioper,file="./oper",form='unformatted',status='old')
    call rdmemdim(ioper)
    close(ioper)

    !==============================!
    ! 2) Reservar memoria módulos  !
    !==============================!
    allocmemory=0
    call alloc_dvrdat
    call alloc_grddat    
    call alloc_psidef

    ! (no añadimos lógica; mantengo tus módulos tal cual)

        !==============================!
    ! Flags para leer DVR arrays   !
    !==============================!
    dvrdata(1) = .true.
    dvrdata(2) = .true.
    if (lfft) then 
        dvrdata(6)  = .true.
        dvrdata(7)  = .true.
        dvrdata(11) = .true.
    endif

    !==============================!
    ! 3) dvrinfo (lee defs/labels) !
    !==============================!
    filename = "./dvr"
    ! ilbl = index(filename,' ')-1
    open(idvr,file=filename,form='unformatted',status='old')
    write(6,'(a)') ' Reading DVR data from '//filename
    chkdvr=1
    call dvrinfo(lerr,chkdvr)
    rewind(idvr)
    chkdvr=0
    call rddvrdef(idvr,chkdvr,ndof1,fdvr)
    close(idvr)

    nrp = zort(2)-zort(1)
    nrg = zort(3)-zort(2)
    nth = zort(4)-zort(3)
    print*,'nrp,nrg,nth=',nrp,nrg,nth

    !==============================!
    ! 4) Abrir psi y leer headers  !
    !==============================!
    open(ipsi,file="./psi",form='unformatted',status='old')
    read(ipsi) filever(ipsi)
    chkdvr=0
    chkgrd=1
    chkpsi=1
    chkdat=1
    call rdpsiinfo(ipsi,chkdvr,chkgrd,chkpsi,chkdat)
    call rdpsidef(ipsi,check)
    close(ipsi)

    print*,ortdim,fftdim,expdim,sphdim,dgldim,dvrdim
    ! !==============================!
    ! ! 7) Reservas locales          !
    ! !==============================!
    allocate(rp_grid(nrp),rg_grid(nrg),th_grid(nth))
    allocate(rpp_grid(subdim(1)),rgp_grid(subdim(1)),thp_grid(subdim(1)))
    allocate(newgrid(subdim(1),nmode),auxort(ortdim))
    allocate(ort(ortdim))
    allocate(trafo(ortdim, ortdim))
    allocate(dvrmat(dvrdim, dvrord))
    
    allocate(fftp(fftdim))
    allocate(hin(fftdim))
    allocate(rueck(fftdim))
    allocate(fftfak(fftdim))
    allocate(exphin(expdim))
    allocate(exprueck(expdim))
    allocate(jsph(sphdim))
    allocate(msph(sphdim))
    allocate(kinsph(sphdim))
    allocate(psi(dgldim), spsi(dgldim))
    allocate(jindx(dgldim))
    allocate(trajst(dgldim))
    allocate(adgwp(dgldim))
    allocate(gwpdep(dgldim))
    allocate(psigrd(ortdim), spsigrd(ortdim), workc(ortdim))

    !==============================!
    ! 5) Leer arrays del DVR       !
    !==============================!
    open(idvr,file="./dvr",form='unformatted',status='old')
    chkdvr=2
    call rddvr(ort,trafo,dvrmat,fftp,hin,rueck,fftfak,exphin,&
            exprueck,jsph,msph,kinsph,chkdvr)
    close(idvr)

    !==============================!
    ! 8) Leer datos de la psi      !
    !==============================!
    open(ipsi,file="./psi",form='unformatted',status='old')
    read(ipsi) filever(ipsi)

    do i = 1, ntime
        call rdpsi(ipsi,psi,spsi,jindx,agmat,trajst,adgwp,gwpdep)
    enddo

    workcdim = max(workcdim, dgldim + 2*maxspf + 4)
    if (allocated(workc)) deallocate(workc)
    ! print*,workcdim
    allocate(workc(workcdim))
    ! call rdpsigrid(ipsi,psigrd,spsigrd,jindx,workc,lrst)
    close(ipsi)

    !==============================!
    ! 9) Preparar grid transform   !
    !==============================!
    auxort = ort

    print*,"subdim= ",subdim
    print*,"ortdim= ",ortdim
    print*,"zort= ",zort
    print*,"shape zort= ",shape(ort)
    rp_grid => auxort(zort(1):(zort(2)-1))
    rg_grid => auxort(zort(2):(zort(3)-1))
    th_grid => auxort(zort(3):(zort(4)-1))

    rpp_grid => newgrid(:,1)
    rgp_grid => newgrid(:,2)
    thp_grid => newgrid(:,3)

    m1 = 30.9737620d0
    m2 = 2.0141017778d0
    m3 = 2.0141017778d0

    open(unit=100,file="grid_test.dat",status="replace",action="write")
    do ith = 1, nth
        do ig = 1, nrg
            do ip = 1, nrp
                idx = ip + (ig-1)*nrp + (ith-1)*nrp*nrg
                rp = rp_grid(ip)
                rg = rg_grid(ig)
                theta = th_grid(ith)
                call transform_coords(rp,rg,theta,rpp,rgp,thetap,m1,m2,m3)
                rpp_grid(idx) = rpp
                rgp_grid(idx) = rgp
                thp_grid(idx) = thetap

                theta = theta*180.0d0/pi
                write(100,fmt="(6(f12.8,1x),2x,2f18.12)") rpp, rgp, thetap, rp, rg, theta, psi(idx)
            enddo
        enddo
        write(100,*) " "
    enddo
    close(100)

    call psiout(psi,dicht4,time,jindx,agmat,trajst,gwpdep)

    deallocate(auxort, newgrid, ort)
    deallocate(trafo, dvrmat, fftp, hin, rueck, fftfak, exphin)
    deallocate(exprueck, jsph, msph, kinsph, jindx, trajst)
    deallocate(adgwp, gwpdep, psigrd, spsigrd, workc)

    call alloc_dvrdat
    call alloc_grddat    
    call alloc_psidef

contains

    !==============================!
    ! Transforms the coordinates   !
    ! from Jacobi reactants to     !
    ! Jacobi products              !
    !==============================!

    subroutine Jac_to_ic(rg, rp, theta, r12, r13, r23, m1, m2, m3)
        implicit none
        real(dop), intent(in)  :: rg, rp, theta, m1, m2, m3
        real(dop), intent(out) :: r12, r13, r23
        real(dop) :: thrad, mred2, mred3, pi 

        pi = dacos(-1.d0)
        ! thrad = theta*pi/180.d0
        thrad = theta
        mred2 = m2/(m2+m3)
        mred3 = m3/(m2+m3)

        r12 = sqrt(rg*rg + (mred3*rp)**2 + 2.d0*rg*mred3*rp*cos(thrad))
        r13 = sqrt(rg*rg + (mred2*rp)**2 - 2.d0*rg*mred2*rp*cos(thrad))
        r23 = rp
    end subroutine Jac_to_ic

    subroutine ic_to_Jac(rg, rp, theta, r12, r13, r23, m1, m2, m3)
        implicit none
        real(dop), intent(in)  :: r12, r13, r23, m1, m2, m3
        real(dop), intent(out) :: rg, rp, theta
        real(dop) :: cbeta, ctheta, mred2, d1, pi

        pi = dacos(-1.d0)
        mred2 = m2/(m1+m2)
        d1    = mred2*r13
        cbeta = (r12*r12 + r13*r13 - r23*r23)/(2.d0*r12*r13)

        ! print*,"IC to JAC ",r12,r13,r23,cbeta,d1,mred2

        rp    = r13
        rg    = sqrt(r12*r12 + d1*d1 - 2.d0*d1*r12*cbeta)
        ctheta = (r12*r12 - d1*d1 - rg*rg)/(-2.d0*d1*rg)

        theta  = acos(ctheta)
        theta  = theta*180.d0/pi
    end subroutine ic_to_Jac

    subroutine transform_coords(rgR, rpR, thetaR, rgP, rpP, thetaP, m1, m2, m3)
        implicit none
        real(dop), intent(in)  :: rgR, rpR, thetaR, m1, m2, m3
        real(dop), intent(out) :: rgP, rpP, thetaP
        real(dop) :: r12, r13, r23

        call Jac_to_ic(rgR, rpR, thetaR, r12, r13, r23, m1, m2, m3)
        call ic_to_Jac(rgP, rpP, thetaP, r12, r13, r23, m1, m2, m3)
    end subroutine transform_coords

end program psitrJAF
