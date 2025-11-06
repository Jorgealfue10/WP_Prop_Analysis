program psi_transform 
    
    use decimal, only: long, dop, sip
    use global
    use globmemmod
    use aglobal
    use memdimmod
    use keyunitsmod
    use dvrdatmod
    use griddatmod
    use psidef
    use datenmod
    use compdatmod
    use operdef
    use rddvrmod
    use rdfilesmod, only: rdpsi
    use iofile, only: rdpsiinfo

    !==============================!
    !    Variables declarations    !
    !==============================!

    implicit none
    ! --- indexes and counters ---
    integer(long) :: i, j, k
    integer(long) :: ip, ig, ith, idx
    integer(long) :: ipsi, ierr
    integer(long) :: chkdvr, chkgrd, chkpsi, chkdat
    integer(long) :: check, check_dvr
    integer(long) :: tot_dim

    ! --- dynamic arrays (integers) ---
    integer(long), allocatable :: zetf(:,:)     ! índices por modo y estado
    integer(long), allocatable :: jindx(:)
    integer(long), allocatable :: workc(:)

    ! --- complex arrays ---
    complex(dop), allocatable :: psi(:), spsi(:)
    complex(dop), allocatable :: agmat(:)
    complex(dop), allocatable :: trajst(:), adgwp(:), gwpdep(:)
    complex(dop), allocatable :: psigrd(:), spsigrd(:)

    ! --- real arrays ---
    real(dop), allocatable :: newgrid(:,:)      ! grid completo transformado

    ! --- pointers ---
    real(dop), pointer :: rp_grid(:), rg_grid(:), th_grid(:)
    real(dop), pointer :: rpp_grid(:), rgp_grid(:), thp_grid(:)

    ! --- real scalar variables ---
    real(dop) :: rp, rg, theta                  ! coordenadas iniciales
    real(dop) :: rpp, rgp, thetap               ! coordenadas transformadas
    real(dop) :: r12, r13, r23                  ! coordenadas internas
    real(dop) :: m1, m2, m3                     ! masas atómicas

    ! --- namefiles ---
    character(len=200) :: filein, filename, dname

    ! --- logicals ---
    logical :: lrddvr

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
    open(ipsi,file=filename,form='unformatted',status='old',err=900)
    rewind(ipsi)
    read(ipsi,err=999) filever(ipsi)
    call rdpsiinfo(ipsi,chkdvr,chkgrd,chkpsi,chkdat)

    !==============================!
    !     Forming psi indexes      !
    !==============================!
    
    allocate(zetf(nmode,nstate))

    call get_indexes(zetf,tot_dim)

    !==============================!
    !       Allocating vars        !
    !==============================!

    allocate(psi(tot_dim),spsi(tot_dim),jindx(tot_dim),agmat(tot_dim))
    allocate(trajst(tot_dim),adgwp(tot_dim),gwpdep(tot_dim))
    allocate(workc(tot_dim),lrst(tot_dim),psigrd(tot_dim),spsigrd(tot_dim))
    allocate(newgrid(griddim,nmode))
    
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
        open(idvr,file=filename,form='unformatted',status='old',err=1000)
        
        call dvrinfo(lrddvr,chkdvr)
        call rddvr(ort,trafo,dvrmat,fftp,hin,rueck,fftfak,exphin,&
                exprueck,jsph,msph,kinsph,chkdvr)

        check_dvr=1
        call rddvrdef(idvr,check_dvr)
        close(idvr)
    endif

    !==============================!
    !      Transforming grid       !
    !==============================!

    rp_grid => ort(zort(1):zort(1)+subdim(1)-1)
    rg_grid => ort(zort(2):zort(2)+subdim(2)-1)
    th_grid => ort(zort(3):zort(3)+subdim(3)-1)

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

                call transform_coord(rp,rg,theta,rpp,rgp,thetap,m1,m2,m3)

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
        d1    = mred2*r12
        cbeta = (r12*r12 + r13*r13 - r23*r23)/(2.d0*r12*r13)

        rp    = r12
        rg    = sqrt(r13*r13 + d1*d1 - 2.d0*d1*r13*cbeta)
        ctheta = (r12*r12 - d1*d1 - rg*rg)/(2.d0*d1*rg)
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

    subroutine get_indexes(zetf,tot_dim):
        implicit none
        integer, allocatable, intent(out) :: zetf(:)
        integer(long), intent(out) :: tot_dim
        integer :: i,j
        integer :: val_ind

        allocate(zetf(nmode,nstate))

        tot_dim = griddim * nstate
        
        val_ind = tot_dim
        do i=1,nmode
            do j=1,nstate
                zetf(i,j) = val_ind
                val_ind = val_ind + griddim
            enddo
        enddo
        
    end subroutine get_indexes

end program psi_transform