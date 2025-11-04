program psi_transform 
    
    use decimal, only: long, dop
    use global
    use rdfilesmod
    use iofile
    use iopsifile

    implicit none

    !-----------------------------------------------------------------------
    ! Initialize variables and set defaults
    !-----------------------------------------------------------------------
    call default
    call adefault

    !-----------------------------------------------------------------------
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Declare files to be read
    ! set flag to .true. if file is to be read
    !
    !     lrddvr    DVR file
    !     lrdoper   OPER file
    !-----------------------------------------------------------------------
    lrddvr=.true.
    lrdoper=.true.

    !-----------------------------------------------------------------------
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Declare resources required (other than those needed to read a file)
    !
    ! ncomp  : no. of comparison data sets to be used.
    !          See auswutil.F for details of use of comparison sets
    ! lpsi   : wavefunction
    ! lpsi1  : second wavefunction
    ! lpsigrd : wavefunction to be used in grid representation (rather
    !           than spf)
    ! ldmat  : density matrices for psi
    ! ldmat1 : density matrices for psi1
    ! lort   : grid coordinates
    ! lgpop  : grid populations
    ! ltrafo : FBR/DVR trafo matrices
    !-----------------------------------------------------------------------
    ncomp=0
    lpsi=.true.
    lpsi1=.false.
    lpsigrd=.false.
    ldmat=.true.
    ldmat1=.false.
    lort=.true.
    lgpop=.false.
    ltrafo=.true.

    dnaem="./"
    filename=dname//'/psi'
    open(ipsi,file=filein,form='unformatted',status='old',err=900)
    !-----------------------------------------------------------------------
    ! get memory sizes 
    !-----------------------------------------------------------------------
    call rdmemdim(ipsi)
    ! if (lrdoper) then
    !     filename=operfile
    !     ilbl=index(filename,' ')-1
    !     open(ioper,file=filename,form='unformatted',status='old',err=900)
    !     call rdmemdim(ioper)
    !     close(ioper)
    ! endif

    !-----------------------------------------------------------------------
    ! Allocate memory 
    !-----------------------------------------------------------------------
    allocmemory=0
    call alloc_dvrdat
    call alloc_grddat
    ! call alloc_operdef
    call alloc_psidef
    ! call alloc_runpropmod
    ! allocate(lfast(maxrdf))
    ! lfastall=.false.
    ! lfast(:)=.false.

    !-----------------------------------------------------------------------
    ! Read system information from property file(s)
    ! chkflags = 0   read information and ignore
    !            1   read information and store
    !            2   read information and check
    ! First read file version number
    !-----------------------------------------------------------------------
    chkdvr=1
    chkgrd=1
    chkpsi=1
    chkdat=1
    chkprp=0
    filename=filein
    rewind(ipsi)
    read(ipsi,err=999) filever(ipsi)
    call rdpsiinfo(ipsi,chkdvr,chkgrd,chkpsi,chkdat)
    call rdpsi(ipsi,psi,spsi,jindx,agmat,trajst,adgwp,gwpdep)

    check=1
    call rdpsidef(ipsi,check)
    close(ipsi)

    if (lrddvr) then
        chkdvr=2
        filename=dname(1:dlaenge)//'/dvr'
        open(idvr,file=filename,form='unformatted',status='old',err=1000)
        
        call dvrinfo(lrddvr,chkdvr)
        call rddvr(ort,trafo,dvrmat,fftp,hin,rueck,fftfak,exphin,&
                exprueck,jsph,msph,kinsph,chkdvr)

        check_dvr=1
        call rddvrdef(idvr,check_dvr)
        close(idvr)
    endif



    contains

    ! ########################################
    ! Transforms the coordinates in the wavefunction
    ! from Jacobi reactants to Jacobi products
    ! ########################################

    subroutine Jac_to_ic(rg,rp,theta,r12,r13,r23,m1,m2,m3)
        implicit none
        real(long),intent(in) :: rg,rp,theta,m1,m2,m3
        real(long),intent(out) :: r12,r13,r23
        real(long) :: thrad,mred1,mred2,mred3
        
        thrad=theta
        
        mred2 = (m2/(m2+m3))
        mred3 = (m3/(m2+m3))

        r12=dsqrt(rg*rg + mred3*mred3*rp*rp + 2.d0*rg*mred3*rp*dcos(thrad))
        r13=dsqrt(rg*rg + mred2*mred2*rp*rp - 2.d0*rg*mred2*rp*dcos(thrad))
        r23=rp

    end subroutine get_ic

    subroutine ic_to_Jac(rg,rp,theta,r12,r13,r23,m1,m2,m3)
        implicit none
        real(long),intent(in) :: r12,r13,r23,m1,m2,m3
        real(long),intent(out) :: rg,rp,theta
        real(long) :: cbeta,ctheta,mred1,mred2,mred3,d1
        
        mred2 = (m2/(m1+m2))
        d1 = mred2*r12
        cbeta = (r12*r12+r13*r13-r23*r23)/(2.d0*r12*r13)

        rp = r12
        rg = dsqrt(r13*r13+d1*d1-2.d0*d1*r13*cbeta)
        ctheta = (r12*r12-d1*d1-rg*rg)/(2.d0*d1*rg)
        theta = acos(ctheta)

    end subroutine ic_to_Jac
        
    subroutine transform_coords(rgR,rpR,thetaR,rgP,rpP,thetaP,m1,m2,m3)
        implicit none
        real(long), intent(in) :: rgR,rpR,thetaR,m1,m2,m3
        real(long), intent(out) :: rgP,rpP,thetaP
        real(long) :: r12,r13,r23

        call Jac_to_ic(rgR,rpR,thetaR,r12,r13,r23,m1,m2,m3)
        call ic_to_Jac(rgP,rpP,thetaP,r12,r13,r23,m1,m2,m3)

    end subroutine transform_coords

    !-----------------------------------------------------------------------
    ! EXPINPUT reads the options and the argument.
    !-----------------------------------------------------------------------
    subroutine expinput(stat)

        use global
        use aglobal
        use griddatmod
        use psidef
        use datenmod
        use compdatmod
        use compdat1mod
        use operdef
        use runpropmod
        use hpsimod
        use channels
        use runpropmod

        implicit none

        integer(long), external      :: myiargc
        integer(long)                :: iarg,i,ierr
        character(len=80)            :: buf
        character(len=*), intent(out):: stat
        character*(c5) filename

    !-----------------------------------------------------------------------
    ! provide default input and output filenames 
    !-----------------------------------------------------------------------
        filein  = ' '
        fileout  = ' '
        operfile  = ' '
        operator = ' '
        unitlab(1)='a.u.'
        name='./'

    !-----------------------------------------------------------------------
    ! set default flags
    !
    ! stat : status of output file   
    !-----------------------------------------------------------------------
        stat   = 'new'
        lgnu   = .false.
        step   = 1
        nskip  = 0
        eunit  = 1.0
        lgnu   = .false.
        lagnu  = .false.
        lprint = .false.
        lnwexpect = .false.
        xmax   = 0.0_dop
        xmin   = 1.0_dop
        tfin   = 1.0e+50_dop

        !-----------------------------------------------------------------------
        ! read options 
        !-----------------------------------------------------------------------
        iarg = myiargc() 
        if (iarg .lt. 1) goto 200 

        i = 0
100     i=i+1
        call mygetarg(i,buf)
        if( buf(1:1) .eq. '-' ) then
            if( buf(1:4) .eq. '-op ') then
            i = i+1
            call mygetarg(i,buf)
            operfile = buf
            else if(buf(1:3) .eq. '-f ') then
            i = i+1
            call mygetarg(i,buf)
            filein = buf
            else if( buf(1:3) .eq. '-o ' ) then
            i = i+1
            call mygetarg(i,buf)
            fileout = buf
            else if(buf(1:3) .eq. '-i ') then
            i=i+1
            call mygetarg(i,buf)
            name = buf
            else if (buf(1:3) .eq. '-r ' ) then
            lnwexpect = .true.
            else if (buf(1:5) .eq. '-rst ' ) then
            irst = ipsi
            else if (buf(1:3) .eq. '-t ' ) then
            i = i+1
            call mygetarg(i,buf)
            read(buf,*) tfin

            else if( buf(1:4) .eq. '-DB ' ) then
            i = i+1
            call mygetarg(i,buf)
            ddpath = buf

            else if (buf(1:3) .eq. '-u ' ) then
            i = i+1
            call mygetarg(i,buf)
            call keylowc(buf)
            eunit = 1.0
            call keyunits(eunit,buf,ierr)
            if (ierr .ne. 0) then
                write(6,'(a)') 'Unknown Units with -u :'//buf
                stop
            endif
            unitlab(1)=buf(1:c1)
            eunit=1.0/eunit
            else if( buf(1:3) .eq. '-x ' ) then
            i = i+1
            call mygetarg(i,buf)
            read(buf,*,err=998) xmin
            i = i+1
            call mygetarg(i,buf)
            read(buf,*,err=998) xmax
            else if( buf(1:3) .eq. '-n ' ) then
            i = i+1
            call mygetarg(i,buf)
            read(buf,*,iostat=ierr) step
            if (ierr .ne. 0  .or.  step .le. 0) then
                write(6,'(2a/,a)')&
                    ' The  argument following the option -n',&
                    ' must be a positive integer. Not : ', buf
                stop
            end if
            else if( buf(1:5) .eq. '-skip' ) then
            i = i+1
            call mygetarg(i,buf)
            read(buf,*,iostat=ierr) nskip
            if (ierr .ne. 0  .or.  nskip .le. 0) then
                write(6,'(2a/,a)')&
                    ' The  argument following the option -skip',&
                    ' must be a positive integer. Not : ', buf
                stop
            end if
            else if( buf(1:6) .eq. '-nofs ' ) then
            fs = 1.0_dop
            else if( buf(1:3) .eq. '-G ' ) then
            lpgrid = .true.
            else if( buf(1:2) .eq. '-w' ) then
            lovwr=.true.
            else if( buf(1:2) .eq. '-g' ) then
            lgnu = .true.
            else if( buf(1:3) .eq. '-a ' ) then
            lgnu  = .true.
            lagnu = .true.
            stat  = 'unknown'
            else if( buf(1:3) .eq. '-p ' ) then
            lgnu   = .true.
            lprint = .true.
            else if( buf(1:2) .eq. '-h' .or. buf(1:2) .eq. '-?' ) then
            write(6,'(78a1)')  ('-',i=1,78)
            write(6,'(1x,2a,/1x,2a,/17x,a)')&
        'Purpose: ',&
        'Calculates the time evolution of an expectation value.',&
        'Usage  : ',&
        'expect [-op -f -i -o -nofs -g -a -r -rst -x -t -p -u -n -skip',&
                '-MC -Macc -G -w -ver -h -?] [operator].'
            write(6,'(/,a)') 'Options : '
            write(6,'(2a,/,9x,a)') &
                        ' -op FILE : The operator is read '&
                    ,'from file FILE rather than from ./oper'
            write(6,'(2a,/,9x,a)') &
                        ' -f FILE : The wavefunction is read '&
                    ,'from file FILE rather than from ./psi'&
            ,'  The string FILE may be a relative or a full path-name.'
            write(6,'(2a,/,11x,a/,11x,a/,11x,a)')&
                        ' -o FILE : The output is written to file ',&
                        'FILE.pl rather than to ./expect.pl',&
                'The string FILE may be a relative or a full path-name.',&
                'If FILE=no, i.e. "-o no", no output file will be opened,',&
                'and the output is directed to screen.'
            write(6,'(a)') ' -i DIR : data is stored in directory DIR'
            write(6,'(2a)') 
            write(6,'(2a)')&
                        ' -nofs   : No transformation to fs.',&
                        " Use when 'time-not-fs' was set in mctdh."
                write(6,'(2a)')&
                        ' -g      : GNUplot comand lines are written to ',&
                        'the output file(s).'
            write(6,'(2a)')&
                        ' -a      : GNUplot is called automatically.',&
                        ' (sets -g and -w).'
            write(6,'(a)')&
                        ' -DB PATH : Sets path to DD DB relative to current location.'
            write(6,'(a)')&
                        ' -p      : Print GNUplot to lpr (sets -g). '
            write(6,'(a)')&
                ' -MC     : Use Monte-Carlo integration (useful for G-MCTDH evaluation). '
            write(6,'(a)')&
                ' -MCacc R: Set accuracy of Monte-Carlo integration to R (default R=0.01). '
            write(6,'(2a)')&
                        ' -r      : Do not divide the expectation values',&
                        ' by norm^2.'
            write(6,'(2a)')&
                        ' -rst    : Use restart rather than psi file.'
            write(6,'(2a)')&
                        ' -x xi xf: Set the abscissa length.'
            write(6,'(2a)')&
                        ' -t tfin : The expectation values are computed',&
                        ' only up to time=tfin.'
            write(6,'(2a)')&
                        ' -u UNIT : The unit used is UNITs',&
                        ' (Default is a.u.).'
            write(6,'(2a)') ' -n step : Compute the property',&
                            ' only every step-th output.'
            write(6,'(2a)') ' -skip nskip : Skip the first nskip',&
                            ' wavefunctions.'
            write(6,'(2a)') &
                        ' -G      : Plot gridlines',&
                        ' (meaningfull only if -g or -a is set).'
            write(6,'(2a)')&
                        ' -w      : An existing property file is ',&
                        'overwritten.'
            write(6,'(a)')&
                        ' -ver    : Version information about the program.'
            write(6,'(a,/a)')&
                        ' -h      : Print this help text. ',&
                        ' -?      : Print this help text. '
            write(6,'(/2a/,2a)') ' The argument operator is the name',&
                            ' of the operator, as specified in an',&
                            ' HAMILTONIAN-SECTION_xxx, for which the',&
                            ' expectation value is required.'
            write(6,'(2a)') ' If no argument is given, ',&
                    'the user is prompted for the missing argument.'
            write(6,'(78a1)')  ('-',i=1,78)
            stop
            else if( buf(1:4) .eq. '-ver' ) then
            call wrversion(6)
            stop
            else
            write(6,'(2a)') ' unknown option : ', buf
            stop
            end if
            if (iarg .le. i) goto 200
            goto 100
        end if
        if (iarg .lt. i) goto 200

        !-----------------------------------------------------------------------
        ! get compulsory argument.
        !-----------------------------------------------------------------------
            if (iarg .eq. i) call mygetarg(i,operator)

        !-----------------------------------------------------------------------
        ! End of loop over options and arguments.
        !-----------------------------------------------------------------------
        200  continue

        !-----------------------------------------------------------------------
        ! set status variable to overwrite files
        !-----------------------------------------------------------------------
            if (lovwr) stat='unknown'

        !-----------------------------------------------------------------------
        ! convert filenames to absolute
        !-----------------------------------------------------------------------
        laenge=index(name,' ')-1
        call abspath(name,laenge)

        if( filein .eq. ' ') then
            if(ipsi.eq.irst) then
            filein=name(1:laenge)//'/restart'
            else
            filein=name(1:laenge)//'/psi'
            endif
        endif
        laein=index(filein,' ')-1
        if (filein(max(1,laein-6):laein) .eq. 'restart') irst = ipsi
        if (filein(laein:laein)   .eq. '/') then
            if(ipsi.eq.irst) then
            filein = filein(1:laein)//'restart'
            laein  = laein + 7
            else
            filein = filein(1:laein)//'psi'
            laein  = laein + 3
            endif
        end if

        call abspath(filein,laein)

        dname=name
        dlaenge=laenge

        if (operfile .eq. '  ') then
            operfile=name(1:laenge)//'/oper'
            opflaenge=index(operfile,' ')-1
        else
            opflaenge=index(operfile,' ')-1
            call abspath(operfile,opflaenge)
        endif

        if (fileout .eq. '  ') then
            fileout=name(1:laenge)//'/expect.pl'
            laeout=index(fileout,'  ')-1
        else if (fileout(1:3) .ne. 'no ') then
            laeout=index(fileout,'  ')-1
            call abspath(fileout,laeout)
        endif

        !-----------------------------------------------------------------------
        ! Open log file channel
        !-----------------------------------------------------------------------
        if(fileout .ne. 'no' ) then
            filename=fileout
            ilbl=index(filename,' ')-1
            if(filename(ilbl-2:ilbl) .eq. '.pl') then
            filename=filename(1:ilbl-3)
            endif
            filename=trim(filename)//'.log'
            open (unit=ilog,file=filename,form='formatted')
            logisopen=.true.
        else
            logisopen=.false.
        endif

        return

        998 continue
        write(6,*) 'Two real numbers expected after -x  option '
        stop
    end subroutine expinput

end program psi_transform