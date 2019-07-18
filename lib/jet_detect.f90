module jet_detect
  implicit none
  !
  ! Parameters
  integer(kind=4), parameter :: ni = 4, nr = 8     ! Numerical precision used throughout
  real(kind=nr), parameter :: nan = 0.0_nr/0.0_nr  ! Not-a-Number, used for missing values. Requires the gfortran compiler option -fno-range-check
  !
  ! Configuration variables
  real(kind=nr) :: jetint_thres = 5.5e-9_nr, &  ! K-threshold to isolate well-defined wind maxima
    &                 searchrad = 1.5_nr, &     ! Maximum distance between two points along the jet axis, in grid point indicies
    &                    minlen = 1.0e6         ! Minimum length of the jet lines in meters
  character, parameter :: cr = char(13_1)       ! Carriage Retirm is ASCII code 13
  logical :: grid_cyclic_ew = .false.           ! Is the grid periodic in east-west, i.e. in the x direction
contains
  !
  !@ Find jetaxes by zero-shear condition and a mask for well-defined wind maxima
  !@
  !@ The mask for well-defined maxima is defined by d/dn(U * dU/dn) < K. The
  !@ threshold K is to be defined in config module of dynfor.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ no : int
  !@     Maximum number of points along jet axes
  !@ nf : int
  !@     Maximum number of jet axis lines
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     u wind velocity component, typically filtered to T42 resolution.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     v wind velocity component, typically filtered to T42 resolution.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ 
  !@ np.ndarray with shape (nz,no,3) and dtype float64
  !@     List of points points belonging to jet axes for each time step. For each point, 
  !@     (0) the j-index, (1) the i-index and (2) the wind speed is saved.
  !@ np.ndarray with shape (nz,nf) and dtype float64
  !@     List of point indexes marking the beginning of jet axes within the point array.
  subroutine run_jet_detect(ja,jaoff,nx,ny,nz,no,nf,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), & 
                 &                dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: ja(nz,no,3_ni), jaoff(nz,nf)
    integer(kind=ni) :: nx,ny,nz, no, nf
    !f2py depend(nx,ny,nz) v
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nz) ja, jaoff
    !
    !real(kind=nr) :: us(nz,ny,nx), vs(nz,ny,nx)
    real(kind=nr) :: dsheardx(nz,ny,nx), dsheardy(nz,ny,nx), &
                 &   jetint(nz,ny,nx), shear(nz,ny,nx), ff(nz,ny,nx)
    integer(kind=ni) :: i,j,k, ip1,im1
    ! -----------------------------------------------------------------
    !
    write(*,*) 'preparing'
    !
    ff = sqrt( u(:,:,:)**2_ni + v(:,:,:)**2_ni )
    !
    ! Based on shear in natural coordinates, and the absolute wind speed u calculate d/dn(u * du/dn)
    call shear_nat(shear, nx,ny,nz, u,v, dx,dy)
    call ddx(dsheardx, nx,ny,nz, shear*ff, dx,dy)
    call ddy(dsheardy, nx,ny,nz, shear*ff, dx,dy)
    jetint(:,:,:) = (u(:,:,:)*dsheardy(:,:,:) - v(:,:,:)*dsheardx(:,:,:))/ff(:,:,:)
    !
    where(jetint > -jetint_thres)
       shear = nan
    end where
    !
    ! Save the wind speed along the jet axis in the third field along with the
    ! jet axis position. 
    !  (the original jetint is not needed anymore after the masking `where` block)
    jetint(:,:,:) = ff(:,:,:)
    !
    call line_locate(ja,jaoff, nx,ny,nz,no,nf, shear,jetint, dx,dy)
    !
    return
  end subroutine
  !
  !@ Calculate the local wind shear in natural coordinates
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated wind shear.
  subroutine shear_nat(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: ff(nz,ny,nx), ffx(nz,ny,nx), ffy(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    ff = sqrt(u(:,:,:)**2_ni + v(:,:,:)**2_ni)
    call grad(ffx, ffy, nx,ny,nz, ff, dx,dy)
    !
    res(:,:,:) = (u(:,:,:)*ffy(:,:,:) - v(:,:,:)*ffx(:,:,:))/ff(:,:,:)
  end subroutine
  !
  !@ Calculates partial x derivative: ddatdx = partial(dat)/partial(x)
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ longitude for non-cyclic grids.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     x-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`ddx_o4`, :meth:`ddx_on_q`
  subroutine ddx(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    forall(k = 1_ni:nz, j = 1_ni:ny, i = 2_ni:nx-1_ni)
       res(k,j,i) = (dat(k,j,i+1_ni)-dat(k,j,i-1_ni))/dx(j,i)
    end forall
    if (grid_cyclic_ew) then
       forall(k = 1_ni:nz, j = 1_ni:ny)
          res(k,j,1_ni) = (dat(k,j,  2_ni)-dat(k,j     ,nx))/dx(j,1_ni)
          res(k,j,nx  ) = (dat(k,j  ,1_ni)-dat(k,j,nx-1_ni))/dx(j,nx)
       end forall
    else 
       forall(k = 1_ni:nz, j = 1_ni:ny)
          res(k,j,1_ni) = nan
          res(k,j,nx  ) = nan
       end forall
    end if
  end subroutine
  !
  !@ Calculates partial y derivative: ddatdy = partial(dat)/partial(y)
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ latitude.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     y-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`ddy_o4`, :meth:`ddy_on_q`
  subroutine ddy(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    forall(k = 1_ni:nz, j = 2_ni:ny-1_ni, i = 1_ni:nx)
       res(k,j,i) = (dat(k,j+1_ni,i)-dat(k,j-1_ni,i))/dy(j,i)
    end forall
    forall(k = 1_ni:nz, i = 1_ni:nx)
       res(k,1_ni,i)=nan
       res(k,ny,i)=nan
    end forall
  end subroutine
  !
  !@ Calculates second partial derivative in x and y directions::
  !@
  !@     bx, by = grad(dat)
  !@
  !@ The routine uses 2nd-order centered differences. Returns NaN on first and last 
  !@ latitude and longitude for non-cyclic grids.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Data array
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     x-derivative of ``dat``.
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     y-derivative of ``dat``.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`grad_3d`
  subroutine grad(resx,resy,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: resx(nz,ny,nx), resy(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) resx, resy
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call ddx(resx,nx,ny,nz,dat,dx,dy)
    call ddy(resy,nx,ny,nz,dat,dx,dy)
    !
  end subroutine
  !
  !
  ! ###################################################################
  ! Internal helper functions, not intended to be called directly
  ! ###################################################################
  ! 
  subroutine line_locate(lines,lnoff,nx,ny,nz,no,nf,lnloc,lnint,dx,dy)
    real(kind=nr), intent(in) :: lnloc(nz,ny,nx), lnint(nz,ny,nx), &
                 &               dx(ny,nx), dy(ny,nx) 
    real(kind=nr), intent(out) :: lines(nz,no,3_ni), lnoff(nz,nf)
    integer(kind=ni) :: nx,ny,nz, no, nf
    !f2py depend(nx,ny,nz) lnint
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nz) lines, lnoff
    !
    integer(kind=ni), parameter :: nn = 30000_ni
    !
    real   (kind=nr), allocatable :: reci(:,:), recj(:,:), linelen(:)
    integer(kind=ni), allocatable :: lineptcnt(:)
    !
    real(kind=nr) :: zeroloc(nn,2_ni)
    integer(kind=ni) :: i,j,k, di, m, n, zerocnt, ptcnt, linecnt, off
    ! -----------------------------------------------------------------
    !
    do k = 1_ni,nz
       write(*,'(I5,A4,I5,A)', advance='no') k, 'of', nz, cr
       !
       ! find fronts
       call find_zeroloc(lnloc(k,:,:), nx,ny,nn, zeroloc,zerocnt)
       ! 
       allocate(recj(zerocnt,zerocnt), reci(zerocnt,zerocnt), lineptcnt(zerocnt), linelen(zerocnt) )
       reci(:,:) = nan
       recj(:,:) = nan
       linecnt    = 0_ni ! number of lines
       ptcnt      = 0_ni ! total numer of points
       lineptcnt(:) = 0_ni ! number of points per line
       call linejoin(zerocnt, nf*5_ni, nx,ny, zeroloc(:,2_ni), zeroloc(:,1_ni), &
               & recj, reci, linelen, lineptcnt, dx,dy) 
       !
       off = 0_ni
       do n = 1_ni,zerocnt
          if ( isnan(recj(n,1_ni)) ) then
             exit
          end if
          !
          ! filter fronts by length
          if (linelen(n) >= minlen) then
          !if (.true.) then
             linecnt = linecnt + 1_ni
             ptcnt = ptcnt + lineptcnt(n)
             !
             ! check if results larger than output array
             if (ptcnt > no) then
                write(*,*) 'Found more points than output array allows: ', no
                stop 1
             end if
             if (linecnt > nf) then
                write(*,*) 'Found more fronts than output array allows: ', nf
                stop 1
             end if
             !
             ! write into output arrays fr and froff
             do m = 1_ni,lineptcnt(n)
                lines(k,off+m,1_ni) = reci(n,m)
                lines(k,off+m,2_ni) = recj(n,m)
                lines(k,off+m,3_ni) = lnint(k,int(recj(n,m),ni),int(reci(n,m),ni))
             end do
             lnoff(k,linecnt) = off
             off = off + lineptcnt(n)
          end if
       end do
       ! Save the ending of the last front by saving the beginning of the
       ! first non-existant
       lnoff(k,linecnt+1_ni) = off
       !
       deallocate(reci, recj, lineptcnt, linelen)
       !
    end do ! loop over k
    !
    return
  end subroutine
  !
  ! Find locations where dat is zero by interpolating the 2-dim gridded data
  subroutine find_zeroloc(dat, nx,ny,nn, zeroloc, zerocnt)
    real(kind=nr), intent(in)  :: dat(ny,nx)
    real(kind=nr), intent(out) :: zeroloc(nn,2_ni)
    integer(kind=ni), intent(in) :: nx,ny, nn
    integer(kind=ni), intent(out) :: zerocnt 
    !
    integer(kind=ni) :: i,j, ip1
    ! -----------------------------------------------------------------
    !
    ! (1) scan along x
    zerocnt = 0_ni
    do i = 1_ni,nx-1_ni
       ip1 = i+1_ni
       do j = 1_ni,ny
          ! Zero line hits a grid point
          if (dat(j,i) == 0.0_nr) then
             if (dat(j,ip1) == 0.0_nr) then
                zerocnt = zerocnt + 1_ni
                zeroloc(zerocnt,2_ni) = j
                zeroloc(zerocnt,1_ni) = i + 0.5_nr
             else
                zerocnt = zerocnt + 1_ni
                zeroloc(zerocnt,1_ni) = i
                zeroloc(zerocnt,2_ni) = j
             end if
          ! interpolate to find line, i direction first
          else   
             if (dat(j,ip1) /= nan .and. dat(j,i) /= nan) then
                if ((dat(j,i) > 0.0_nr .and. dat(j,ip1) < 0.0_nr) .or. &
                  & (dat(j,i) < 0.0_nr .and. dat(j,ip1) > 0.0_nr)) then
                   zerocnt = zerocnt + 1_ni
                   zeroloc(zerocnt,2_ni) = j
                   zeroloc(zerocnt,1_ni) = i + dat(j,i)/(dat(j,i) - dat(j,ip1))
                end if    ! diff signs
             end if    ! Missin ip1
          end if    ! zero exactly at grid point
       end do   
    end do
    !
    ! take into account periodicity in x
    if ( grid_cyclic_ew ) then
       i = nx
       ip1 = 1_ni
       do j = 1_ni,ny
          ! Zero line hits a grid point
          if (dat(j,i) == 0.0_nr) then
             if (dat(j,ip1) == 0.0_nr) then
                zerocnt = zerocnt + 1_ni
                zeroloc(zerocnt,2_ni) = j
                zeroloc(zerocnt,1_ni) = i + 0.5_nr
             else
                zerocnt = zerocnt + 1_ni
                zeroloc(zerocnt,1_ni) = i
                zeroloc(zerocnt,2_ni) = j
             end if
          ! interpolate to find line, i direction first
          else   
             if (dat(j,ip1) /= nan .and. dat(j,i) /= nan) then
                if ((dat(j,i) > 0.0_nr .and. dat(j,ip1) < 0.0_nr) .or. &
                  & (dat(j,i) < 0.0_nr .and. dat(j,ip1) > 0.0_nr)) then
                   zerocnt = zerocnt + 1_ni
                   zeroloc(zerocnt,2_ni) = j
                   zeroloc(zerocnt,1_ni) = i + dat(j,i)/(dat(j,i) - dat(j,ip1))
                end if    ! diff signs
             end if    ! Missin ip1
          end if    ! zero exactly at grid point
       end do
    end if
    ! 
    ! (2) scan along y
    do i = 1_ni,nx
       do j = 1_ni,ny-1_ni
          ip1 = j + 1_ni
          ! Zero line hits a grid point
          if (dat(j,i) == 0.0_nr) then
             if (dat(ip1,i) == 0.0_nr) then
                zerocnt = zerocnt + 1_ni
                zeroloc(zerocnt,2_ni) = j + 0.5_nr          
                zeroloc(zerocnt,1_ni) = i
             else
                zerocnt = zerocnt + 1_ni
                zeroloc(zerocnt,1_ni) = i
                zeroloc(zerocnt,2_ni) = j
             end if
          ! interpolate to find line, j direction first
          else   
             if (dat(ip1,i) /= nan .and. dat(j,i) /= nan) then
                if ((dat(j,i) > 0.0_nr .and. dat(ip1,i) < 0.0_nr) .or. &
                  & (dat(j,i) < 0.0_nr .and. dat(ip1,i) > 0.0_nr)) then
                   zerocnt = zerocnt + 1_ni
                   zeroloc(zerocnt,2_ni) = j + dat(j,i)/(dat(j,i) - dat(ip1,i))
                   zeroloc(zerocnt,1_ni) = i
                end if
             end if    ! Missin ip1
          end if    ! zero exactly at grid point
       end do   
    end do
    !
    return
  end  subroutine find_zeroloc
  !
  ! Join a cloud of frontal points into frontal lines
  subroutine linejoin(cnt,nstruct_max,nx,ny,jidx,iidx,recj,reci,linelen,lineptcnt,dx,dy)
    real(kind=nr), intent(in)  :: jidx(cnt), iidx(cnt), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: recj(cnt,cnt), reci(cnt,cnt), linelen(cnt)
    integer(kind=ni), intent(in) :: cnt, nx,ny, nstruct_max
    integer(kind=ni), intent(out) :: lineptcnt(cnt)
    !
    real   (kind=nr) :: dists(cnt,cnt), longestdist
    integer(kind=ni) :: n,m,l, nn, inter(cnt,cnt), oidx(cnt), structcnt(nstruct_max), &
            &           donecnt,prevcnt, nstruct, startidx,endidx, &
            &           longestpath(2_ni), linecnt, pos,startpos,endpos
    logical :: used(cnt)
    ! -----------------------------------------------------------------
    !
    dists(:,:) = -1.0_nr ! distance metric between all pairs of points (sorted coherent structures);
    !                    ! -1.0_nr represents infinite distance ( = no connection)
    inter(:,:) = 0_ni    ! intermediate point defining the shortest path between each pair points
    used(:) = .false.    ! Used in the DFS? Is eqivalent to: Already assigned to a choherent structure?
    oidx(:) = 0_ni       ! Mapping of the new index (used in the dists matrix) to the old index used in jidx/iidx
    !
    !
    ! 1. Step: Find coherent structures via a Depth-First-Search
    nstruct = 0_ni
    donecnt = 0_ni
    do n = 1_ni,cnt
       if ( used(n) ) cycle
       !
       nstruct = nstruct + 1_ni
       if ( nstruct > nstruct_max ) then
          write(*,*) 'Found more than NSTRUCT_MAX=', nstruct_max, 'structures'
          stop 1
       end if
       prevcnt = donecnt
       call depth_first_search(cnt,donecnt, nx,ny, n,used,dists, oidx, jidx,iidx, dx,dy)
       structcnt(nstruct) = donecnt - prevcnt
    end do
    !
    ! For each of the structures
    linecnt = 0_ni
    startidx = 1_ni
    do nn = 1_ni,nstruct
       endidx = startidx - 1_ni + structcnt(nn)
       !
       ! 2. Step: Find shortest paths between all pairs of points in the structures 
       !          using the Floyd-Warshall algorithm 
       do n = startidx,endidx
          do m = startidx,endidx
             do l = startidx,endidx
                ! If l and m are connected via n
                if ( dists(n,m) >= 0.0_nr .and. dists(n,l) >= 0.0_nr ) then
                   ! Check if the connection via n is shorter than the
                   ! previously known connection (if any)
                   if ( dists(n,m) + dists(n,l) < dists(m,l) .or. dists(m,l) < 0.0_nr ) then
                      dists(m,l) = dists(n,m) + dists(n,l)
                      inter(m,l) = n
                   end if
                end if
             end do
          end do
       end do
       !
       ! 3. Step: Reconstruct path for the longest shortest paths within each structure 
       !          using the "inter(mediate)" information from step 2
       longestpath = maxloc(dists(startidx:endidx,startidx:endidx))
       longestdist = maxval(dists(startidx:endidx,startidx:endidx))
       startpos = longestpath(1_ni) + startidx - 1_ni
       endpos   = longestpath(2_ni) + startidx - 1_ni
       !
       if ( longestdist >= minlen ) then
          ! Initialise and record start position
          linecnt = linecnt + 1_ni
          lineptcnt(linecnt) = 1_ni
          linelen(linecnt) = longestdist
          pos = startpos
          recj(linecnt,1_ni) = jidx(oidx(pos))
          reci(linecnt,1_ni) = iidx(oidx(pos))
          !
          ! Recursively Construct and record the path to the end position
          call construct_path(cnt, startpos,endpos, inter, oidx, jidx,iidx, &
                  & lineptcnt(linecnt), recj(linecnt,:), reci(linecnt,:))
          !
          !write(*,*) linecnt, linelen(linecnt), lineptcnt(linecnt)
       end if
       !
       startidx = endidx + 1_ni
    end do
    !
    return
  end subroutine linejoin
  !
  !
  recursive subroutine depth_first_search(cnt,donecnt, nx,ny, n,used,dists, oidx, jidx,iidx, dx,dy)
    real(kind=nr), intent(in)  :: jidx(cnt), iidx(cnt), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(inout) :: dists(cnt,cnt)
    integer(kind=ni), intent(in) :: cnt, nx,ny, n
    integer(kind=ni), intent(inout) :: oidx(cnt), donecnt
    logical, intent(inout) :: used(cnt)
    !
    real   (kind=nr) :: di,dj, dist,distji
    integer(kind=ni) :: i,j, m, nn,nm
    ! -----------------------------------------------------------------
    !
    donecnt = donecnt + 1_ni
    nn = donecnt
    used(n) = .true.
    oidx(nn) = n
    dists(nn,nn) = 0.0_nr
    !
    i = iidx(n) ! Cast to integer
    j = jidx(n) ! Cast to integer
    !
    do m = 1_ni,cnt
       if ( used(m) ) cycle
       !
       di = iidx(m) - iidx(n)
       dj = jidx(m) - jidx(n)
       !
       ! Take into account periodicity in x
       if ( grid_cyclic_ew ) then
          if ( di > nx/2.0_nr ) then
             di = di - nx
          elseif ( di < -nx/2.0_nr ) then
             di = di + nx
          end if
       end if
       !
       dist = sqrt( (dx(j,i)*di/2.0_nr)**2.0_nr + &
                  & (dy(j,i)*dj/2.0_nr)**2.0_nr   )
       distji = sqrt(di**2_ni + dj**2_ni)
       !
       ! if we found a pair of close points
       if ( distji <= searchrad ) then
          nm = donecnt + 1_ni
          !write(*,*) n, m, nn, nm
          dists(nn,nm) = dist
          dists(nm,nn) = dist
          call depth_first_search(cnt,donecnt, nx,ny, m,used,dists, oidx, jidx,iidx, dx,dy)
       end if
    end do
  end subroutine
  !
  !
  recursive subroutine construct_path(cnt, startpos,endpos, inter, oidx,jidx,iidx, lineptcnt, recj, reci)
    real(kind=nr), intent(in)  :: jidx(cnt), iidx(cnt)
    real(kind=nr), intent(out) :: recj(cnt), reci(cnt)
    integer(kind=ni), intent(in) :: cnt, inter(cnt,cnt), oidx(cnt), startpos, endpos
    integer(kind=ni), intent(inout) :: lineptcnt
    !
    integer(kind=ni) :: interpos
    ! -----------------------------------------------------------------
    !
    interpos = inter(startpos,endpos)
    ! Nothing in between: directly connected
    if ( interpos == 0_ni ) then
       lineptcnt = lineptcnt + 1_ni
       recj(lineptcnt) = jidx(oidx(endpos))
       reci(lineptcnt) = iidx(oidx(endpos))
    ! No connection at all
    else if ( inter(startpos,endpos) < 0_ni ) then
       write(*,*) 'ERROR: Trying to join unconnected points'
       stop 1
    ! At least one point in between start and end
    else
       call construct_path(cnt, startpos,interpos, inter, oidx, jidx,iidx, lineptcnt, recj, reci)
       call construct_path(cnt, interpos,endpos, inter, oidx, jidx,iidx, lineptcnt, recj, reci)
    end if
  end subroutine
  !
end module
