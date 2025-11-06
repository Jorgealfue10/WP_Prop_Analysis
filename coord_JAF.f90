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

    !==============================!
    !    Variables declarations    !
    !==============================!

    implicit none
    ! --- indexes and counters ---
    integer(long) :: i, j, k
    integer(long) :: ip, ig, ith, idx
    integer(long) :: ierr
    integer(long) :: chkdvr, chkgrd, chkpsi, chkdat
    integer(long) :: check, check_dvr
    integer(long) :: tot_dim
    integer(long) :: ndof1
    integer(long), allocatable :: fdvr(:)

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
    real(dop) :: m1, m2, m3                     ! masas atómicas

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
    logical :: lrst

    !==============================!
    !   Variables initialization   !
    !==============================!

    chkdvr = 1
    chkgrd = 1
    chkpsi = 1
    chkdat = 1
    lrddvr = .true.

    !==============================!
    !     Reading psi metadata     !
    !==============================!

    dname="./"
    filename=dname//'psi'
    open(ipsi,file="./psi",form='unformatted',status='old')
    rewind(ipsi)
    read(ipsi) filever(ipsi)
    call rdpsiinfo(ipsi,chkdvr,chkgrd,chkpsi,chkdat)

    tot_dim = griddim * nstate
    
    !==============================!
    !       Allocating vars        !
    !==============================!

    allocate(rp_grid(subdim(1)),rg_grid(subdim(2)),th_grid(subdim(3)))
    allocate(rpp_grid(griddim),rgp_grid(griddim),thp_grid(griddim))
    allocate(newgrid(griddim,nmode),auxort(griddim))
    allocate(ort(ortdim))
    allocate(trafo(ortdim, ortdim))
    allocate(dvrmat(dvrdim, dvrdim))
    allocate(fftp(fftdim))
    allocate(hin(fftdim))
    allocate(rueck(fftdim))
    allocate(fftfak(fftdim))
    allocate(exphin(expdim))
    allocate(exprueck(expdim))
    allocate(jsph(sphdim))
    allocate(msph(sphdim))
    allocate(kinsph(sphdim))
    allocate(psi(tot_dim), spsi(tot_dim))
    allocate(jindx(tot_dim))
    allocate(trajst(tot_dim))
    allocate(adgwp(tot_dim))
    allocate(gwpdep(tot_dim))
    allocate(psigrd(tot_dim), spsigrd(tot_dim), workc(tot_dim))
    allocate(fdvr(maxdim))

    
    !==============================!
    !       Reading psi data       !
    !==============================!

    call rdpsi(ipsi,psi,spsi,jindx,agmat,trajst,adgwp,gwpdep)
    call rdpsigrid(ipsi,psigrd,spsigrd,jindx,workc,lrst)

    check=1
    call rdpsidef(ipsi,check)
    close(ipsi)

    !==============================!
    !       Reading DVR data       !
    !==============================!
    if (lrddvr) then
        chkdvr=2
        filename=dname//'dvr'
        open(idvr,file=filename,form='unformatted',status='old')
        
        call dvrinfo(lrddvr,chkdvr)
        call rddvr(ort,trafo,dvrmat,fftp,hin,rueck,fftfak,exphin,&
                exprueck,jsph,msph,kinsph,chkdvr)

        check_dvr=1
        call rddvrdef(idvr,check_dvr,ndof1,fdvr)
        close(idvr)
    endif
    auxort = ort

    !==============================!
    !      Transforming grid       !
    !==============================!

    rp_grid => auxort(zort(1):(zort(1)+subdim(1)-1))
    rg_grid => auxort(zort(2):(zort(2)+subdim(2)-1))
    th_grid => auxort(zort(3):(zort(3)+subdim(3)-1))

    rpp_grid => newgrid(:,1)
    rgp_grid => newgrid(:,2)
    thp_grid => newgrid(:,3)

    do ith = 1, subdim(3)
        do ig = 1, subdim(2)
            do ip = 1, subdim(1)
                idx = ip + (ig-1)*subdim(1) + (ith-1)*subdim(1)*subdim(2)

                rp = rp_grid(ip)
                rg = rg_grid(ig)
                theta = th_grid(ith)

                call transform_coords(rp,rg,theta,rpp,rgp,thetap,m1,m2,m3)

                rpp_grid(idx) = rpp
                rgp_grid(idx) = rgp
                thp_grid(idx) = thetap
            enddo
        enddo
    enddo


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
        real(dop) :: thrad, mred2, mred3

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
        real(dop) :: cbeta, ctheta, mred2, d1

        mred2 = m2/(m1+m2)
        d1    = mred2*r13
        cbeta = (r12*r12 + r13*r13 - r23*r23)/(2.d0*r12*r13)

        rp    = r13
        rg    = sqrt(r12*r12 + d1*d1 - 2.d0*d1*r12*cbeta)
        ctheta = (r12*r12 - d1*d1 - rg*rg)/(-2.d0*d1*rg)
        theta  = acos(max(-1.d0,min(1.d0, ctheta)))

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