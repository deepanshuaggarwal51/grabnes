module semicl

! Semiclassical Landau levels from zero-field Fermi contours.
!
! For each band sheet and each energy, the constant-energy contours in the
! moire Brillouin zone are traced with a periodic marching-squares algorithm,
! unwrapped across the periodic boundary, and their enclosed k-space area A is
! integrated.  Onsager quantisation
!
!       A(E_n) = 2*pi*e*B/hbar * (n + gamma)
!
! then gives the field at which Landau level n sits at energy E:
!
!       B[T] = SEMICL_C * A[Ang^-2] / (n + gamma),   SEMICL_C = hbar/(2*pi*e)*1e20
!
! Anchor: a contour enclosing the whole moire BZ at (n+gamma)=1 must return
! Phi_0/A_moire.  For the G/hBN CellSize 55 cell (13.53 nm) that is 26.0867 T,
! which SEMICL_C reproduces to 6 digits.  There is deliberately NO adjustable
! prefactor here -- if a fudge factor ever seems necessary, something upstream
! is wrong.
!
! Ported from P. M. Perez-Piskunow's python notebooks
! (PabloCodes/supermoire-master/contours), with three deliberate changes:
!   * areas are integrated in CARTESIAN k-space (the python integrated in the
!     2*pi-periodic fractional frame and absorbed the Jacobian into a hard
!     coded constant);
!   * the shoelace sign is KEPT, so electron- and hole-like orbits are
!     distinguishable instead of being abs()'d together;
!   * contours that fail to close are flagged as open orbits and excluded from
!     the fan rather than being assigned a meaningless area.
!
! Band ordering by wavefunction overlap (the python `band_order.py`) is NOT
! reproduced: with TAPW's valley projection exactly one energy-sorted sheet
! carries the contour throughout the primary miniband window, and elsewhere the
! union over sorted sheets is still the complete constant-energy surface.

   use mio
   use constants, only : hbar, qe, twopi

   implicit none

   PRIVATE

   ! hbar/(2*pi*e) * 1e20, so that B[T] = SEMICL_C * A[Ang^-2] / (n+gamma)
   real(dp), parameter, public :: SEMICL_C = 10475.7686_dp

   ! electron rest mass [kg], to report m* as a ratio
   real(dp), parameter :: SEMICL_ME = 9.10938291e-31_dp

   ! hbar^2/(2*pi) * 1e20/qe / m_e, so that  m*/m_e = SEMICL_M * dA/dE
   ! with A in Ang^-2 and E in eV.  Derived, not fitted: the same hbar and qe
   ! that give SEMICL_C = 10475.7686 give SEMICL_M = 1.2127626.
   real(dp), parameter, public :: SEMICL_M = &
        hbar*hbar/twopi * 1.0e20_dp / qe / SEMICL_ME

   ! Contour storage is FLAT: every traced point is a distinct boundary cell,
   ! so the total over all contours at one level is bounded by nx*ny.  A
   ! (points x contours) layout would need MAXCONTPTS*MAXCONT doubles, which is
   ! hundreds of MB for no reason.  Contour c occupies cx(cstart(c):cend(c)).
   integer, parameter :: MAXCONT = 4096        ! contours at one (band,energy)

   ! checkpoint format version -- bump if the layout changes
   integer, parameter :: SEMICL_CKPT_VERSION = 1

   ! Semicl.PsiHybrid: compute Psi by the hybrid route as well (see
   ! SemiclSuscFields).  Module-level so that the self-tests, which are the
   ! only place an exact Psi is known, can score it alongside the co-area and
   ! energy-difference routes.  Default .false. -> nothing new is computed or
   ! written and every output is byte-identical to before.
   logical :: psiHybOn = .false.

   ! Semicl.CoareaAnalytic: evaluate the co-area divergence from the ANALYTIC
   ! expansion (see semicl_coarea_div) instead of differencing F*gradE/|gradE|^2.
   ! Default .false. -> the differenced form, bit-for-bit as before.
   logical :: coareaAn = .false.

   ! Semicl.PsiTerms: also write the three terms of Psi = Phi_0 - Phi_1' + Phi_2''
   ! separately.  D is linear, so D(w1 - D(w2)) = D(w1) - D(D(w2)) and the split
   ! is exact; psi must equal t0 - t1 + t2, which the writer checks.  The point
   ! is to see the CANCELLATION: if the terms are large and nearly cancel, their
   ! individual accuracy is irrelevant and no discretisation of D can fix Psi.
   logical :: psiTermsOn = .false.

   public :: SemiclMarchingSquares
   public :: SemiclUnwrap
   public :: SemiclAreaCart
   public :: SemiclSelfTest
   public :: SemiclassicalLL

contains

!=======================================================================
! index wrap helpers (the k-grid is periodic: index nx+1 == index 1)
!=======================================================================
pure integer function iwrap(i, n)
   integer, intent(in) :: i, n
   iwrap = modulo(i-1, n) + 1
end function iwrap

!=======================================================================
! Marching squares on a periodic 2-D grid.
!
! Egrid(nx,ny) holds the band energy at fractional coordinate
! ( (i-1)/nx , (j-1)/ny ).  Traces every closed or open contour at `level`.
!
! Output: cx/cy(1:npts(c), c) fractional coordinates of contour c,
!         nc contours, closed(c) .true. if it returned to its start cell.
!=======================================================================
subroutine SemiclMarchingSquares(Egrid, nx, ny, level, cx, cy, cstart, cend, nc, closed)

   integer,  intent(in)  :: nx, ny
   real(dp), intent(in)  :: Egrid(nx,ny), level
   real(dp), intent(out) :: cx(nx*ny), cy(nx*ny)       ! flat, see note above
   integer,  intent(out) :: cstart(MAXCONT), cend(MAXCONT), nc
   logical,  intent(out) :: closed(MAXCONT)
   integer :: ptot

   logical :: below(nx,ny), visited(nx,ny)
   integer :: code(nx,ny)
   integer :: i, j, i0, j0, ic, jc, di, dj, np, s, dip, djp
   real(dp) :: px, py
   logical  :: ok

   nc   = 0
   ptot = 0
   cstart = 0 ; cend = -1
   closed = .false.

   below = (Egrid < level)

   ! cell code: corner (i,j)=bit0, (i+1,j)=bit1, (i+1,j+1)=bit2, (i,j+1)=bit3
   do j = 1, ny
      do i = 1, nx
         code(i,j) =        merge(1,0, below(i,j))                       &
                     + 2  * merge(1,0, below(iwrap(i+1,nx), j))          &
                     + 4  * merge(1,0, below(iwrap(i+1,nx), iwrap(j+1,ny))) &
                     + 8  * merge(1,0, below(i, iwrap(j+1,ny)))
      end do
   end do

   visited = .false.

   do j0 = 1, ny
      do i0 = 1, nx
         ! cells 0 and 15 have no crossing; 5 and 10 are saddles (skip as seeds)
         if (code(i0,j0) == 0 .or. code(i0,j0) == 15) cycle
         if (visited(i0,j0)) cycle
         if (code(i0,j0) == 5 .or. code(i0,j0) == 10) cycle

         if (nc >= MAXCONT) then
            call MIO_Print('WARNING: MAXCONT reached, contours dropped','semicl')
            return
         end if

         nc = nc + 1
         np = 0
         cstart(nc) = ptot + 1
         ic = i0 ; jc = j0
         ok = .true.
         dip = 0 ; djp = 0

         do s = 1, nx*ny
            visited(ic,jc) = .true.

            call march_step(code(ic,jc),                                     &
                 Egrid(ic,jc), Egrid(iwrap(ic+1,nx),jc),                       &
                 Egrid(iwrap(ic+1,nx),iwrap(jc+1,ny)), Egrid(ic,iwrap(jc+1,ny)),&
                 level, dip, djp, di, dj, ok)
            if (.not. ok) exit                      ! saddle with no incoming dir
            dip = di ; djp = dj

            call edge_point(Egrid, nx, ny, ic, jc, di, dj, level, px, py)

            if (ptot >= nx*ny) exit          ! cannot exceed the boundary-cell count
            ptot = ptot + 1
            np = np + 1
            cx(ptot) = px
            cy(ptot) = py

            ic = iwrap(ic + di, nx)
            jc = iwrap(jc + dj, ny)

            if (ic == i0 .and. jc == j0) then
               closed(nc) = .true.
               exit
            end if
            ! a cell with no crossing means the trace escaped: open contour
            if (code(ic,jc) == 0 .or. code(ic,jc) == 15) exit
            ! A SADDLE cell is legitimately traversed twice by the same
            ! contour -- once per branch -- and the two passes arrive from
            ! different directions, so the decider sends them different ways.
            ! Excluding revisits outright truncated those contours even after
            ! the saddle itself was resolved.  Only a repeat visit to a
            ! NON-saddle cell indicates a genuine loop, so stop on that.
            if (visited(ic,jc) .and. code(ic,jc) /= 5 .and. code(ic,jc) /= 10) exit
         end do

         cend(nc) = ptot
         if (np < 3) then                            ! degenerate, discard
            ptot = cstart(nc) - 1                    ! reclaim the points
            nc = nc - 1
         end if
      end do
   end do

end subroutine SemiclMarchingSquares

!-----------------------------------------------------------------------
! direction table (Perez-Piskunow MARCHING_STEPS).  Saddles (5,10) cannot be
! resolved from grid values alone; we report failure and drop the contour
! rather than guessing, so an ambiguity can never masquerade as an area.
!-----------------------------------------------------------------------
subroutine march_step(c, fa, fb, fc, fd, level, dip, djp, di, dj, ok)
   ! Corners: a=(i,j) bit0, b=(i+1,j) bit1, c=(i+1,j+1) bit2, d=(i,j+1) bit3.
   ! (dip,djp) is the direction we ARRIVED from; needed only for saddles.
   integer,  intent(in)  :: c, dip, djp
   real(dp), intent(in)  :: fa, fb, fc, fd, level
   integer,  intent(out) :: di, dj
   logical,  intent(out) :: ok
   real(dp) :: ga, gb, gc, gd, den, fsad
   ok = .true. ; di = 0 ; dj = 0
   select case (c)
   case (0)  ; di =  1
   case (1)  ; di = -1
   case (2)  ; dj = -1
   case (3)  ; di = -1
   case (4)  ; di =  1
   case (6)  ; dj = -1
   case (7)  ; di = -1
   case (8)  ; dj =  1
   case (9)  ; dj =  1
   case (11) ; dj =  1
   case (12) ; di =  1
   case (13) ; di =  1
   case (14) ; dj = -1
   case (5, 10)
      ! Saddle: the two contour branches can be connected two ways.  Resolve
      ! with the ASYMPTOTIC DECIDER, i.e. the value of the bilinear
      ! interpolant at the saddle point of the cell,
      !     f_sad = (ga*gc - gb*gd) / (ga + gc - gb - gd),
      ! which needs only the four corner values.  (Perez-Piskunow's python
      ! computes a direction here but then raises unconditionally, so his
      ! saddle contours are dropped; on a 48x48 G/hBN grid that discarded 18%
      ! of all contours.)
      ga = fa - level ; gb = fb - level ; gc = fc - level ; gd = fd - level
      den = ga + gc - gb - gd
      if (abs(den) < 1.0e-30_dp) then
         ok = .false. ; return
      end if
      fsad = (ga*gc - gb*gd)/den
      if (c == 5) then                 ! a,c below the level
         if (djp == 1) then            ! arrived moving +y
            di = merge( 1, -1, fsad < 0.0_dp)
         else if (djp == -1) then
            di = merge(-1,  1, fsad < 0.0_dp)
         else
            ok = .false. ; return      ! no incoming direction (seed): skip
         end if
      else                             ! c == 10, b,d below
         if (dip == 1) then            ! arrived moving +x
            dj = merge(-1,  1, fsad < 0.0_dp)
         else if (dip == -1) then
            dj = merge( 1, -1, fsad < 0.0_dp)
         else
            ok = .false. ; return
         end if
      end if
   case default        ! 15: fully inside, no crossing
      ok = .false.
   end select
end subroutine march_step

!-----------------------------------------------------------------------
! Linear interpolation of the crossing on the edge the contour exits through.
! Coordinates are fractional; cell (i,j) spans [(i-1)/nx, i/nx].
!-----------------------------------------------------------------------
subroutine edge_point(Egrid, nx, ny, i, j, di, dj, level, px, py)
   integer,  intent(in)  :: nx, ny, i, j, di, dj
   real(dp), intent(in)  :: Egrid(nx,ny), level
   real(dp), intent(out) :: px, py
   real(dp) :: x0, x1, y0, y1, va, vb, w

   x0 = real(i-1,dp)/real(nx,dp) ; x1 = real(i,dp)/real(nx,dp)
   y0 = real(j-1,dp)/real(ny,dp) ; y1 = real(j,dp)/real(ny,dp)

   if (di == 1) then                       ! exit through the +x edge
      va = Egrid(iwrap(i+1,nx), j)
      vb = Egrid(iwrap(i+1,nx), iwrap(j+1,ny))
      w  = safe_w(level, va, vb)
      px = x1 ; py = y0 + w*(y1-y0)
   else if (di == -1) then                 ! -x edge
      va = Egrid(i, j)
      vb = Egrid(i, iwrap(j+1,ny))
      w  = safe_w(level, va, vb)
      px = x0 ; py = y0 + w*(y1-y0)
   else if (dj == 1) then                  ! +y edge
      va = Egrid(i, iwrap(j+1,ny))
      vb = Egrid(iwrap(i+1,nx), iwrap(j+1,ny))
      w  = safe_w(level, va, vb)
      px = x0 + w*(x1-x0) ; py = y1
   else                                    ! -y edge
      va = Egrid(i, j)
      vb = Egrid(iwrap(i+1,nx), j)
      w  = safe_w(level, va, vb)
      px = x0 + w*(x1-x0) ; py = y0
   end if
end subroutine edge_point

pure real(dp) function safe_w(level, va, vb)
   real(dp), intent(in) :: level, va, vb
   real(dp) :: d
   d = vb - va
   if (abs(d) < 1.0e-30_dp) then
      safe_w = 0.5_dp
   else
      safe_w = (level - va)/d
      if (safe_w < 0.0_dp) safe_w = 0.0_dp
      if (safe_w > 1.0_dp) safe_w = 1.0_dp
   end if
end function safe_w

!=======================================================================
! Unwrap a contour across the periodic boundary (minimum image).
!
! Consecutive marching-squares points are at most one cell apart, so the
! correct periodic image is always the one with |dx|,|dy| < 1/2.  This is the
! minimum-image form of the python `connectedBandsShort` search over
! factor1,factor2 in [-3,3), and is exact rather than range-limited.
!
! Returns wound=.true. if the contour fails to close on itself, i.e. the end
! point differs from the start by a non-zero lattice vector -> OPEN ORBIT,
! which has no cyclotron area.
!=======================================================================
subroutine SemiclUnwrap(x, y, n, wound, wx, wy)
   integer,  intent(in)    :: n
   real(dp), intent(inout) :: x(n), y(n)
   logical,  intent(out)   :: wound
   integer,  intent(out)   :: wx, wy
   integer  :: k
   real(dp) :: dx, dy, x0, y0

   x0 = x(1) ; y0 = y(1)
   do k = 2, n
      dx = x(k) - x(k-1)
      dy = y(k) - y(k-1)
      dx = dx - anint(dx)
      dy = dy - anint(dy)
      x(k) = x(k-1) + dx
      y(k) = y(k-1) + dy
   end do

   ! closing segment: how many lattice vectors did we travel in total?
   dx = x(n) - x0
   dy = y(n) - y0
   wx = nint(dx)
   wy = nint(dy)
   wound = (wx /= 0 .or. wy /= 0)
end subroutine SemiclUnwrap

!=======================================================================
! Signed area of an unwrapped contour, in CARTESIAN k-space [Ang^-2].
!
! Fractional (fx,fy) -> k = fx*b1 + fy*b2 with b1 = rcell(:,1), b2 = rcell(:,2),
! then the shoelace formula WITH the closing segment included.  The sign is
! retained: a negatively oriented contour is a hole-like orbit.
!=======================================================================
subroutine SemiclAreaCart(x, y, n, rcell, area)
   integer,  intent(in)  :: n
   real(dp), intent(in)  :: x(n), y(n), rcell(3,3)
   real(dp), intent(out) :: area
   integer  :: k, kn
   real(dp) :: xa, ya, xb, yb, s

   s = 0.0_dp
   do k = 1, n
      kn = k + 1
      if (kn > n) kn = 1                     ! close the loop
      xa = x(k) *rcell(1,1) + y(k) *rcell(1,2)
      ya = x(k) *rcell(2,1) + y(k) *rcell(2,2)
      xb = x(kn)*rcell(1,1) + y(kn)*rcell(1,2)
      yb = x(kn)*rcell(2,1) + y(kn)*rcell(2,2)
      s = s + (xa*yb - xb*ya)
   end do
   area = 0.5_dp*s
end subroutine SemiclAreaCart

!=======================================================================
! Self-test against cases with an analytic answer.  Run with
!   Calculate.SemiclassicalLL .true.
!   Semicl.SelfTest           .true.
! No Hamiltonian is needed; this exercises only the contour/area/Onsager
! chain, so a failure here is unambiguously a bug in this module.
!=======================================================================
subroutine SemiclSelfTest()

   use constants, only : pi

   integer,  parameter :: NG = 256
   real(dp), allocatable :: Eg(:,:)
   real(dp), allocatable :: cx(:), cy(:), xs(:), ys(:)
   integer  :: cstart(MAXCONT), cend(MAXCONT), nc, wx, wy, np
   logical  :: closed(MAXCONT), wound
   real(dp) :: rcell(3,3), b, area, aexact, lev, hv, kr, worst, rel
   real(dp) :: aa(5), ll(5), kk(5), dAdE, mnum, mex
   integer  :: i, j, c, itest
   real(dp) :: fx, fy, kx, ky
   real(dp), allocatable :: Fg(:,:)
   real(dp) :: dlt, dA, phi, phiex
   logical  :: bmatched
   integer  :: ncell_in
   real(dp), allocatable :: Mg(:,:), gkx(:,:), gky(:,:)
   real(dp) :: delm, dAdEl, gam, gtot, worstb, aniso
   real(dp), allocatable :: E4(:,:,:), F4(:,:,:), Ediab(:,:,:)
   real(dp) :: bmn4(4), bmx4(4), al, e0, c1x, c2x, ee1, ee2, mm, hh, abz
   real(dp) :: ph4, gapw4, areah, phh, gaph
   integer  :: lab4, labh, nexact
   logical  :: okp4, okh
   ! tests 10-11: second-order term (L1c/L2)
   integer,  parameter :: NC10 = 4
   integer,  parameter :: NL10(NC10) = (/1, 1, 1, 3/)
   real(dp), parameter :: BL10(NC10) = (/0.02_dp, 0.04_dp, 0.08_dp, 0.02_dp/)
   real(dp), allocatable :: Wt(:,:,:), Wf(:,:,:), bx(:), by(:)
   real(dp) :: tt, aa0, bl, eps0, bO, rex, rco, red, rct, nu, sg, ss0, hst
   real(dp) :: psi, psied, psiLP, psihyb, eco, eed, ect, sc, relc
   real(dp) :: st0, st1, st2
   real(dp) :: l11c(5), p11c(5), e11c(5), lp11c(5), lx11c(5), t11c(3,5), h11c(5)
   real(dp) :: l11f(5), p11f(5), e11f(5), lp11f(5), lx11f(5), t11f(3,5), h11f(5)
   integer  :: ic, ieh, nn, mm2, npb, nsg, nwk, nfail

#ifdef DEBUG
   call MIO_Debug('SemiclSelfTest',0)
#endif /* DEBUG */

   ! SemiclSelfTest is a SEPARATE entry point from SemiclassicalLL (calc.F90
   ! calls one or the other), so the discretisation flags must be read here too
   ! -- otherwise the self-tests silently score the default route whatever the
   ! input asks for, and a "passes with the flag on" result means nothing.
   ! NOT inside the #ifdef DEBUG above: the production build has no -DDEBUG.
   call MIO_InputParameter('Semicl.PsiHybrid', psiHybOn, .false.)
   call MIO_InputParameter('Semicl.CoareaAnalytic', coareaAn, .false.)
   call MIO_InputParameter('Semicl.PsiTerms', psiTermsOn, .false.)

   call MIO_Print('','semicl')
   call MIO_Print('=== semiclassical self-test =====================','semicl')
   call MIO_Print('  co-area divergence: '// &
        trim(merge('analytic expansion ', 'central differences', coareaAn))// &
        trim(merge(' ; hybrid Psi route on', '                      ', psiHybOn)),'semicl')

   ! ---- a square reciprocal cell keeps the analytic answer simple ----
   b = 0.05362304107122594_dp          ! |b| of the G/hBN CellSize 55 cell
   rcell = 0.0_dp
   rcell(1,1) = b ; rcell(2,2) = b ; rcell(3,3) = 1.0_dp

   allocate(Eg(NG,NG), cx(NG*NG), cy(NG*NG), xs(NG*NG), ys(NG*NG))

   ! ---- test 1: isotropic cone E = hv*|k| centred in the cell ----------
   ! contour at E is a circle of radius E/hv, area pi*(E/hv)^2
   hv = 11.0_dp                        ! eV.Ang, ~ graphene 1.5*a*g0
   do j = 1, NG
      do i = 1, NG
         fx = real(i-1,dp)/real(NG,dp) - 0.5_dp
         fy = real(j-1,dp)/real(NG,dp) - 0.5_dp
         kx = fx*rcell(1,1) + fy*rcell(1,2)
         ky = fx*rcell(2,1) + fy*rcell(2,2)
         Eg(i,j) = hv*sqrt(kx*kx + ky*ky)
      end do
   end do

   worst = 0.0_dp
   call MIO_Print('  test 1: isotropic Dirac cone, circular contours','semicl')
   call MIO_Print('    E[eV]    k_F[Ang^-1]   A_num          A_exact        rel.err','semicl')
   do itest = 1, 5
      kr  = (0.05_dp + 0.05_dp*itest) * (b/2.0_dp)    ! stay inside the cell
      lev = hv*kr
      aexact = pi*kr*kr
      call SemiclMarchingSquares(Eg, NG, NG, lev, cx, cy, cstart, cend, nc, closed)
      area = 0.0_dp
      do c = 1, nc
         np = cend(c) - cstart(c) + 1
         if (np < 3) cycle
         xs(1:np) = cx(cstart(c):cend(c))
         ys(1:np) = cy(cstart(c):cend(c))
         call SemiclUnwrap(xs, ys, np, wound, wx, wy)
         if (wound) cycle
         call SemiclAreaCart(xs, ys, np, rcell, rel)
         area = area + abs(rel)
      end do
      rel = abs(area-aexact)/aexact
      worst = max(worst, rel)
      aa(itest) = area ; ll(itest) = lev ; kk(itest) = kr
      write(*,'(a,f8.4,3x,es12.5,3x,es13.6,2x,es13.6,2x,es10.3)') &
           '    ', lev, kr, area, aexact, rel
   end do
   call MIO_Print('    worst relative area error: '//trim(num2str(worst,6)),'semicl')

   ! ---- test 2: unit anchor, full-BZ orbit -> one flux quantum ---------
   ! A = |b1 x b2| at (n+gamma)=1 must give Phi_0/A_real
   area = abs(rcell(1,1)*rcell(2,2) - rcell(1,2)*rcell(2,1))
   call MIO_Print('','semicl')
   call MIO_Print('  test 2: Onsager unit anchor (full-cell orbit, n+gamma=1)','semicl')
   call MIO_Print('    A_BZ    = '//trim(num2str(area,8))//' Ang^-2','semicl')
   call MIO_Print('    B       = '//trim(num2str(SEMICL_C*area,8))//' T','semicl')
   call MIO_Print('    Phi0/A  = '//trim(num2str(4.135667696e-15_dp/ &
        ((2.0_dp*pi)**2/area*1.0e-20_dp),8))//' T   (must match)','semicl')

   ! ---- test 3: cyclotron mass of the same cone -----------------------
   ! m* = (hbar^2/2pi) dA/dE, and for E = hv*k the exact answer is
   ! dA/dE = 2*pi*k/hv, i.e. m* = SEMICL_M*2*pi*k_F/hv rising linearly in E.
   ! This differences the areas of test 1 exactly as SemiclMassesBand does,
   ! so it checks the constant SEMICL_M and the central difference together
   ! against a closed form.
   call MIO_Print('','semicl')
   call MIO_Print('  test 3: cyclotron mass from dA/dE (central differences)','semicl')
   call MIO_Print('    E[eV]    m*_num/m_e     m*_exact/m_e   rel.err','semicl')
   worst = 0.0_dp
   do itest = 2, 4
      dAdE = (aa(itest+1) - aa(itest-1))/(ll(itest+1) - ll(itest-1))
      mnum = SEMICL_M*dAdE
      mex  = SEMICL_M*2.0_dp*pi*kk(itest)/hv
      rel  = abs(mnum-mex)/mex
      worst = max(worst, rel)
      write(*,'(a,f8.4,3x,es13.6,2x,es13.6,2x,es10.3)') &
           '    ', ll(itest), mnum, mex, rel
   end do
   call MIO_Print('    worst relative mass error: '//trim(num2str(worst,6)),'semicl')

   ! ---- test 4: per-orbit Berry phase (E3) ----------------------------
   ! Massive Dirac cone  E(k) = sqrt((hv k)^2 + D^2)  with the analytic
   ! curvature of the two-band model,
   !     Omega(k) = -(1/2) (hv)^2 D / ((hv k)^2 + D^2)^(3/2),
   ! whose flux inside the orbit at energy E is exactly
   !     Phi_B = -pi (1 - D/E)   ->  -pi well above the gap (gamma -> 1 == 0
   ! mod 1, i.e. the Berry phase pi of a Dirac band) and 0 as the orbit closes
   ! onto the band edge.  This exercises the whole chain the real run uses:
   ! interior mask -> connected component -> plaquette sum -> gamma.
   call MIO_Print('','semicl')
   call MIO_Print('  test 4: per-orbit Berry phase from plaquette fluxes','semicl')
   dlt = 0.02_dp                        ! gap parameter, resolved by this grid
   allocate(Fg(NG,NG))
   dA = abs(rcell(1,1)*rcell(2,2) - rcell(1,2)*rcell(2,1))/real(NG*NG,dp)
   do j = 1, NG
      do i = 1, NG
         ! plaquette centre, consistent with corners (i,j)..(i+1,j+1)
         fx = (real(i-1,dp)+0.5_dp)/real(NG,dp) - 0.5_dp
         fy = (real(j-1,dp)+0.5_dp)/real(NG,dp) - 0.5_dp
         kx = fx*rcell(1,1) + fy*rcell(1,2)
         ky = fx*rcell(2,1) + fy*rcell(2,2)
         kr = hv*hv*(kx*kx + ky*ky) + dlt*dlt
         Fg(i,j) = -0.5_dp*hv*hv*dlt/(kr*sqrt(kr)) * dA
         fx = real(i-1,dp)/real(NG,dp) - 0.5_dp
         fy = real(j-1,dp)/real(NG,dp) - 0.5_dp
         kx = fx*rcell(1,1) + fy*rcell(1,2)
         ky = fx*rcell(2,1) + fy*rcell(2,2)
         Eg(i,j) = sqrt(hv*hv*(kx*kx + ky*ky) + dlt*dlt)
      end do
   end do
   call MIO_Print('    E[eV]    Phi_num[rad]   Phi_exact[rad] rel.err    gamma','semicl')
   worst = 0.0_dp
   do itest = 1, 5
      kr  = (0.05_dp + 0.05_dp*itest) * (b/2.0_dp)
      lev = sqrt(hv*hv*kr*kr + dlt*dlt)
      area = pi*kr*kr
      call SemiclOrbitBerry(Eg, Fg, NG, NG, lev, area, &
                            abs(rcell(1,1)*rcell(2,2) - rcell(1,2)*rcell(2,1)), &
                            phi, bmatched)
      phiex = -pi*(1.0_dp - dlt/lev)
      if (bmatched) then
         rel = abs(phi-phiex)/abs(phiex)
         worst = max(worst, rel)
         write(*,'(a,f8.4,3x,es13.6,2x,es13.6,2x,es10.3,2x,f8.4)') &
              '    ', lev, phi, phiex, rel, 0.5_dp - phi/(2.0_dp*pi)
      else
         write(*,'(a,f8.4,a)') '    ', lev, '   NO COMPONENT MATCHED (bug)'
         worst = 1.0_dp
      end if
   end do
   call MIO_Print('    worst relative Berry-phase error: '//trim(num2str(worst,6)),'semicl')

   ! ---- test 5: the interior summation itself -------------------------
   ! A constant flux per plaquette must integrate to (cells inside)*F, which
   ! is a pure test of the mask -> component -> plaquette bookkeeping with no
   ! physics in it: any off-by-one in the interior test shows up here.
   Fg = 1.0_dp
   kr  = 0.20_dp*(b/2.0_dp)
   lev = sqrt(hv*hv*kr*kr + dlt*dlt)
   area = pi*kr*kr
   call SemiclOrbitBerry(Eg, Fg, NG, NG, lev, area, &
                         abs(rcell(1,1)*rcell(2,2) - rcell(1,2)*rcell(2,1)), phi, bmatched)
   ncell_in = count(Eg < lev)
   call MIO_Print('','semicl')
   call MIO_Print('  test 5: interior bookkeeping with unit plaquette flux','semicl')
   call MIO_Print('    grid points below level : '//trim(num2str(ncell_in)),'semicl')
   call MIO_Print('    plaquette weight summed : '//trim(num2str(phi,8)),'semicl')
   call MIO_Print('    (with corner-fraction weighting these must be EXACTLY '// &
        'equal: every interior point is a corner of 4 plaquettes, each taking '// &
        '1/4 of it, so the weights telescope to the point count)','semicl')
   if (abs(phi - real(ncell_in,dp)) > 1.0e-9_dp*real(ncell_in,dp)) then
      call MIO_Print('    FAIL: weighted plaquette sum /= interior point count','semicl')
   else
      call MIO_Print('    PASS (exact)','semicl')
   end if
   deallocate(Fg)

   ! ---- test 6: first-order quantisation, massive Dirac ELECTRON orbit ----
   ! The Berry phase is NOT the complete first-order correction: the orbital
   ! magnetic moment enters at the same order.  For the two-band massive Dirac
   ! model both are known in closed form,
   !     Omega_c(k) = -(1/2)(hv)^2 D / eps^3 ,   M(k) = +(1/2)(hv)^2 D / eps^2
   ! (M = -Omega*(E_n-E_m)/2 for two bands, so the two share one convention),
   ! and the exact Landau levels are  E_n^2 = D^2 + 2n (hv)^2 / l_B^2, i.e.
   !     A(E) l_B^2 = 2 pi n   exactly.
   ! Since A l_B^2 = 2 pi (N + gamma + delta_m) is what this module computes,
   ! the gate is that  gamma + delta_m  is an INTEGER: the Berry phase alone
   ! leaves gamma = 1 - D/(2E), off by D/(2E) -- 11 to 28 % at these radii --
   ! and the moment must cancel that with no fitted parameter.
   call MIO_Print('','semicl')
   call MIO_Print('  test 6: Berry phase + orbital moment, massive Dirac cone','semicl')
   call MIO_Print('    SEMICL_C*(e/hbar) = '//trim(num2str(SEMICL_C*qe*1.0e-20_dp/hbar,10))// &
        '  must be 1/2pi = '//trim(num2str(1.0_dp/(2.0_dp*pi),10)),'semicl')
   allocate(Fg(NG,NG), Mg(NG,NG), gkx(NG,NG), gky(NG,NG))
   dA = abs(rcell(1,1)*rcell(2,2) - rcell(1,2)*rcell(2,1))/real(NG*NG,dp)
   do j = 1, NG
      do i = 1, NG
         fx = (real(i-1,dp)+0.5_dp)/real(NG,dp) - 0.5_dp
         fy = (real(j-1,dp)+0.5_dp)/real(NG,dp) - 0.5_dp
         kx = fx*rcell(1,1) + fy*rcell(1,2)
         ky = fx*rcell(2,1) + fy*rcell(2,2)
         kr = hv*hv*(kx*kx + ky*ky) + dlt*dlt
         Fg(i,j) = -0.5_dp*hv*hv*dlt/(kr*sqrt(kr)) * dA     ! plaquette flux
         fx = real(i-1,dp)/real(NG,dp) - 0.5_dp
         fy = real(j-1,dp)/real(NG,dp) - 0.5_dp
         kx = fx*rcell(1,1) + fy*rcell(1,2)
         ky = fx*rcell(2,1) + fy*rcell(2,2)
         kr = hv*hv*(kx*kx + ky*ky) + dlt*dlt
         Eg(i,j) = sqrt(kr)                                  ! grid point
         Mg(i,j) = 0.5_dp*hv*hv*dlt/kr
      end do
   end do
   call SemiclGradient(Eg, NG, NG, rcell, gkx, gky)
   call MIO_Print('    E[eV]   gamma     delta_m   gamma+d_m  |dev.int|   dB/B(N=1)'// &
        '   Berry-only dev','semicl')
   worst = 0.0_dp ; worstb = 0.0_dp
   do itest = 1, 5
      kr  = (0.05_dp + 0.05_dp*itest) * (b/2.0_dp)
      lev = sqrt(hv*hv*kr*kr + dlt*dlt)
      call selftest_orbit(Eg, Fg, Mg, gkx, gky, NG, rcell, lev, &
                          area, phi, delm, dAdEl, bmatched)
      if (.not. bmatched) then
         write(*,'(a,f8.4,a)') '    ', lev, '   NO CONTOUR/COMPONENT (bug)'
         worst = 1.0_dp ; cycle
      end if
      gam  = 0.5_dp - phi/(2.0_dp*pi)
      gtot = gam + delm
      rel  = abs(gtot - anint(gtot))
      worst = max(worst, rel)
      worstb = max(worstb, rel/(1.0_dp + anint(gtot)))
      write(*,'(a,f8.4,2x,f9.5,1x,f9.5,1x,f9.5,2x,es10.3,2x,es10.3,4x,es10.3)') &
           '    ', lev, gam, delm, gtot, rel, rel/(1.0_dp+anint(gtot)), &
           abs(gam - anint(gam))
   end do
   call MIO_Print('    worst |gamma+delta_m - integer| : '//trim(num2str(worst,6)),'semicl')
   call MIO_Print('    worst relative error in B (N=1) : '//trim(num2str(worstb,6)),'semicl')
   call MIO_Print('    (the last column is what Berry alone leaves behind: it is'// &
        ' 100x larger, so this test discriminates)','semicl')

   ! ---- test 7: the same for a HOLE orbit (massive Dirac valence band) -----
   ! Time reversal within the model: Omega and the traced sense both flip, the
   ! moment does not (M is the same for both bands of a two-band model).  The
   ! hole path through SemiclOrbitBerry/SemiclOrbitMoment is a different branch
   ! -- interior mask 'above', sgn = -1 -- and must land on the SAME exact
   ! Landau spectrum, here with the integer 0, i.e. this is the band edge that
   ! carries the anomalous n=0 level.
   call MIO_Print('','semicl')
   call MIO_Print('  test 7: same, HOLE orbit (valence sheet of the same model)','semicl')
   Fg = -Fg
   Eg = -Eg
   call SemiclGradient(Eg, NG, NG, rcell, gkx, gky)
   call MIO_Print('    E[eV]   gamma     delta_m   gamma+d_m  |dev.int|   signed_area','semicl')
   worst = 0.0_dp
   do itest = 1, 5
      kr  = (0.05_dp + 0.05_dp*itest) * (b/2.0_dp)
      lev = -sqrt(hv*hv*kr*kr + dlt*dlt)
      call selftest_orbit(Eg, Fg, Mg, gkx, gky, NG, rcell, lev, &
                          area, phi, delm, dAdEl, bmatched)
      if (.not. bmatched) then
         write(*,'(a,f8.4,a)') '    ', lev, '   NO CONTOUR/COMPONENT (bug)'
         worst = 1.0_dp ; cycle
      end if
      gam  = 0.5_dp - phi/(2.0_dp*pi)
      gtot = gam + delm
      rel  = abs(gtot - anint(gtot))
      worst = max(worst, rel)
      write(*,'(a,f8.4,2x,f9.5,1x,f9.5,1x,f9.5,2x,es10.3,2x,es13.6)') &
           '    ', lev, gam, delm, gtot, rel, area
      if (area > 0.0_dp) call MIO_Print('    FAIL: valence orbit came out '// &
           'electron-like (positive signed area)','semicl')
   end do
   call MIO_Print('    worst |gamma+delta_m - integer| : '//trim(num2str(worst,6)),'semicl')

   ! ---- test 8: the 1/|grad E| weighting, on an ANISOTROPIC cone ----------
   ! With a circular contour every weighting rule gives the same answer, so
   ! tests 6/7 do not exercise it.  For E = hv*sqrt(kx^2 + (ky/r)^2) the
   ! contour is an ellipse of semi-axes E/hv and r*E/hv, |grad E| varies
   ! around it, and the line integral must still return
   !     dA/dE = 2 pi r E / (hv)^2
   ! exactly.  A uniform (perimeter-like) weighting fails this by ~20 % at
   ! r = 2.  M is set constant here so that the test isolates the weight.
   call MIO_Print('','semicl')
   call MIO_Print('  test 8: line-integral dA/dE on an anisotropic cone (weighting)','semicl')
   aniso = 2.0_dp
   Mg = 1.0_dp
   do j = 1, NG
      do i = 1, NG
         fx = real(i-1,dp)/real(NG,dp) - 0.5_dp
         fy = real(j-1,dp)/real(NG,dp) - 0.5_dp
         kx = fx*rcell(1,1) + fy*rcell(1,2)
         ky = fx*rcell(2,1) + fy*rcell(2,2)
         Eg(i,j) = hv*sqrt(kx*kx + (ky/aniso)**2)
      end do
   end do
   call SemiclGradient(Eg, NG, NG, rcell, gkx, gky)
   call MIO_Print('    E[eV]    dA/dE_line     dA/dE_exact    rel.err     A_num/A_exact','semicl')
   worst = 0.0_dp
   do itest = 1, 3
      kr  = (0.04_dp + 0.03_dp*itest) * (b/2.0_dp)     ! semi-axis in x
      lev = hv*kr
      call selftest_orbit(Eg, Fg, Mg, gkx, gky, NG, rcell, lev, &
                          area, phi, delm, dAdEl, bmatched)
      if (.not. bmatched) then
         write(*,'(a,f8.4,a)') '    ', lev, '   NO CONTOUR (bug)'
         worst = 1.0_dp ; cycle
      end if
      aexact = 2.0_dp*pi*aniso*lev/(hv*hv)
      rel = abs(dAdEl - aexact)/aexact
      worst = max(worst, rel)
      write(*,'(a,f8.4,3x,es13.6,2x,es13.6,2x,es10.3,2x,f12.8)') &
           '    ', lev, dAdEl, aexact, rel, abs(area)/(pi*aniso*kr*kr)
   end do
   call MIO_Print('    worst relative dA/dE error: '//trim(num2str(worst,6)),'semicl')
   deallocate(Fg, Mg, gkx, gky)

   ! ---- test 9: strong-breakdown composite orbit (E7) ---------------------
   ! Two parabolic pockets that OVERLAP IN ENERGY, eps1 = al|k-c1|^2 and
   ! eps2 = al|k-c2|^2 + E0, coupled by a tiny D, between two flat bands so the
   ! spectrum has one full gap below and one above.  The energy-sorted sheets
   ! trace the UNION and the LENS of the two circles; in the P -> 1 limit the
   ! orbit is the pair of DIABATIC circles instead, so the composite area must be
   !     A = pi*E/al + pi*(E-E0)/al     (grid-counting accuracy),
   ! Three EXACT identities ride along, each a pure bookkeeping check:
   !   (a) the composite count equals the count of the diabatic functions
   !       themselves -- per k the two sets of energies are the same;
   !   (b) with unit plaquette flux the composite phase equals the number of
   !       points below E on the two straddling sheets (test 5, summed);
   !   (c) the hole composite from the upper gap holds the rest: |A_h|+A_e = 2 A_BZ.
   call MIO_Print('','semicl')
   call MIO_Print('  test 9: strong-breakdown composite orbit (two overlapping pockets)','semicl')
   al = 1000.0_dp ; e0 = 0.02_dp ; dlt = 1.0e-4_dp
   c1x = -0.1_dp*b ; c2x = 0.1_dp*b
   abz = abs(rcell(1,1)*rcell(2,2) - rcell(1,2)*rcell(2,1))
   allocate(E4(NG,NG,4), F4(NG,NG,4), Ediab(NG,NG,2))
   do j = 1, NG
      do i = 1, NG
         fx = real(i-1,dp)/real(NG,dp) - 0.5_dp
         fy = real(j-1,dp)/real(NG,dp) - 0.5_dp
         kx = fx*rcell(1,1) + fy*rcell(1,2)
         ky = fx*rcell(2,1) + fy*rcell(2,2)
         ee1 = al*((kx-c1x)**2 + ky*ky)
         ee2 = al*((kx-c2x)**2 + ky*ky) + e0
         mm = 0.5_dp*(ee1 + ee2)
         hh = sqrt(0.25_dp*(ee1 - ee2)**2 + dlt*dlt)
         E4(i,j,1) = -1.0_dp ; E4(i,j,2) = mm - hh
         E4(i,j,3) = mm + hh ; E4(i,j,4) = 10.0_dp
         Ediab(i,j,1) = ee1 ; Ediab(i,j,2) = ee2
      end do
   end do
   do i = 1, 4
      bmn4(i) = minval(E4(:,:,i)) ; bmx4(i) = maxval(E4(:,:,i))
   end do
   F4 = 1.0_dp
   worst = 0.0_dp ; nexact = 0
   call MIO_Print('    E[eV]    A_composite    A_diabatic     rel.err     label  hole label','semicl')
   do itest = 1, 3
      lev = 0.04_dp + 0.04_dp*real(itest,dp)
      call bd_composite(E4, NG, NG, 4, lev, bmn4, bmx4, 1.0e-4_dp, -2.0_dp, 20.0_dp, abz, &
                        F4, 1, 1, area, ph4, lab4, gapw4, okp4)
      call bd_composite(E4, NG, NG, 4, lev, bmn4, bmx4, 1.0e-4_dp, -2.0_dp, 20.0_dp, abz, &
                        F4, 1, 2, areah, phh, labh, gaph, okh)
      aexact = pi*lev/al + pi*(lev - e0)/al
      rel = abs(area - aexact)/aexact
      worst = max(worst, rel)
      if (lab4 /= 1002 .or. labh /= 2003 .or. .not. okp4) nexact = nexact + 1
      ! (a) the relabeling-free count against the diabatic functions
      if (abs(area/abz*real(NG*NG,dp) - &
              real(count(Ediab(:,:,1) < lev) + count(Ediab(:,:,2) < lev),dp)) > 0.5_dp) &
         nexact = nexact + 1
      ! (b) unit plaquette flux, summed over the straddling sheets 2 and 3
      if (abs(ph4 - real(count(E4(:,:,2) < lev) + count(E4(:,:,3) < lev),dp)) > 1.0e-9_dp*ph4) &
         nexact = nexact + 1
      ! (c) electron + hole composites cover the two middle sheets exactly
      if (abs(abs(areah) + area - 2.0_dp*abz) > 1.0e-12_dp*abz) nexact = nexact + 1
      write(*,'(a,f8.4,3x,es13.6,2x,es13.6,2x,es10.3,2x,i6,2x,i6)') &
           '    ', lev, area, aexact, rel, lab4, labh
   end do
   call MIO_Print('    worst relative area error (grid counting): '//trim(num2str(worst,6)), &
        'semicl')
   if (nexact == 0) then
      call MIO_Print('    anchors, diabatic count, unit-flux sum, e+h cover: PASS (exact)', &
           'semicl')
   else
      call MIO_Print('    FAIL: '//trim(num2str(nexact))//' exact identities violated','semicl')
   end if
   deallocate(E4, F4, Ediab)

   ! ---- test 10: second-order term (L1c/L2), square-lattice band edges ----
   ! h = -2t(cos kx a + cos ky a) has exact Landau levels near its bottom
   !     E_n = t[-4 + m b - (m^2 + 1) b^2/16] + O(b^4),   m = 2n+1,  b = eBa^2/hbar
   ! (checked against exact Hofstadter levels in figures/check_l2_limits.py):
   ! Onsager alone gives the m^2, the chi0' term the "+1".  At the energy of the
   ! exact level, Onsager alone would read the field b_O of
   ! -4 + m b_O - m^2 b_O^2/16, so the exact ratio b/b_O is a closed form, and
   ! semicl's  B_2/B_1 = 2 nu/(nu + sqrt(nu^2 + 4 dchi S))  must reproduce it.
   ! The RATIO cancels the polygon-area error, which is larger than the effect
   ! itself, and tests Psi, its sign, the unit conversion and the quadratic
   ! solve together.  Single band: w_0 = w_2 = 0 and w_1 = e_xx e_yy - e_xy^2 =
   ! 4 t^2 a^4 cos(kx a) cos(ky a).  The HOLE orbit at the band top (mirror
   ! energy, the same w_1) goes through the sgn = -1 branch and must give the
   ! same ratio.  Control: Psi -> -Psi must fail.
   call MIO_Print('','semicl')
   call MIO_Print('  test 10: second-order term, square lattice (single band: Landau-Peierls)','semicl')
   allocate(Wt(NG,NG,4), Wf(NG,NG,5), gkx(NG,NG), gky(NG,NG), bx(NG*NG), by(NG*NG))
   tt  = 0.01_dp                        ! hopping [eV]
   aa0 = 2.0_dp*pi/b                    ! lattice constant [Ang] of this square cell
   do j = 1, NG
      do i = 1, NG
         fx = real(i-1,dp)/real(NG,dp) - 0.5_dp
         fy = real(j-1,dp)/real(NG,dp) - 0.5_dp
         Eg(i,j) = -2.0_dp*tt*(cos(2.0_dp*pi*fx) + cos(2.0_dp*pi*fy))
         Wt(i,j,1) = 0.0_dp ; Wt(i,j,3) = 0.0_dp
         Wt(i,j,2) = 4.0_dp*tt*tt*aa0**4*cos(2.0_dp*pi*fx)*cos(2.0_dp*pi*fy)
         Wt(i,j,4) = Wt(i,j,2)
      end do
   end do
   nsg = 0
   call SemiclGradient(Eg, NG, NG, rcell, gkx, gky)
   call SemiclSuscFields(Wt, Eg, gkx, gky, NG, NG, rcell, psiHybOn, psiTermsOn, Wf, nsg)
   call MIO_Print('    e/h  n   b      b/b_O-1 exact  coarea        ediff'// &
        '         err_coarea  err_ediff   err(-Psi)','semicl')
   worst = 0.0_dp ; nfail = 0
   do ieh = 1, 2
      do ic = 1, NC10
         nn = NL10(ic) ; bl = BL10(ic) ; mm2 = 2*nn + 1
         eps0 = real(mm2,dp)*bl - real(mm2*mm2+1,dp)*bl*bl/16.0_dp
         bO   = 8.0_dp*(1.0_dp - sqrt(1.0_dp - eps0/4.0_dp))/real(mm2,dp)
         rex  = bl/bO
         lev  = tt*(-4.0_dp + eps0)
         if (ieh == 2) lev = -lev
         call pick_orbit(Eg, NG, rcell, lev, bx, by, npb, area)
         if (npb == 0) then
            call MIO_Print('    NO CONTOUR (bug)','semicl') ; nfail = nfail + 1 ; cycle
         end if
         ! h moves the contour by ~1/4 grid step, as dE does on the TAPW grid
         hst = 0.25_dp*(b/real(NG,dp))*2.0_dp*tt*aa0*sqrt(eps0)
         nwk = 0
         call SemiclOrbitSusc(Eg, Wt, Wf, gkx, gky, NG, NG, bx, by, npb, area, rcell, &
                              lev, hst, psiHybOn, psiTermsOn, psi, psied, psihyb, psiLP, &
                              st0, st1, st2, nwk)
         sg  = sign(1.0_dp, area)
         nu  = real(nn,dp) + 0.5_dp
         ss0 = SEMICL_C*abs(area)
         rco = semicl_onsager_B(ss0, nu, semicl_chi_perT(psi, sg))/(ss0/nu)
         red = semicl_onsager_B(ss0, nu, semicl_chi_perT(psied, sg))/(ss0/nu)
         rct = semicl_onsager_B(ss0, nu, semicl_chi_perT(-psi, sg))/(ss0/nu)
         eco = abs(rco - rex)/abs(rex - 1.0_dp)
         eed = abs(red - rex)/abs(rex - 1.0_dp)
         ect = abs(rct - rex)/abs(rex - 1.0_dp)
         worst = max(worst, eco)
         if (ect < 1.0_dp) nfail = nfail + 1
         write(*,'(a,a1,2x,i2,2x,f5.3,2x,es12.5,2x,es12.5,2x,es12.5,2x,es10.3,2x,es10.3,2x,es10.3)') &
              '     ', merge('e','h',ieh==1), nn, bl, rex-1.0_dp, rco-1.0_dp, red-1.0_dp, &
              eco, eed, ect
         if (ieh == 1 .and. ic == 1) call MIO_Print('    Psi/a^2 = '// &
              trim(num2str(psi/aa0**2,6))//' (co-area), '//trim(num2str(psied/aa0**2,6))// &
              ' (energy differences); band-bottom limit 3pi/2 = '// &
              trim(num2str(1.5_dp*pi,6)),'semicl')
      end do
   end do
   call MIO_Print('    worst relative error in b/b_O - 1 (co-area): '//trim(num2str(worst,6)),'semicl')
   if (nfail == 0 .and. worst < 0.05_dp) then
      call MIO_Print('    PASS (< 5 %; the wrong-sign control fails every case)','semicl')
   else
      call MIO_Print('    FAIL','semicl')
   end if

   ! ---- test 11: the second-order term vanishes for a massive Dirac cone ----
   ! h = hv (kx sx + ky sy) + D sz has exact levels E_n^2 = D^2 + 2n (hv)^2/l_B^2,
   ! i.e. Onsager + first order exactly, so Psi = 0 outside the gap.  The three
   ! terms Phi_0, -dPhi_1/dE, d2Phi_2/dE2 are each large and must cancel; the
   ! single-band part alone does not vanish -- the interband terms cancel it --
   ! and has the closed form  Psi_LP = 6 pi (hv)^2 D^2 / E^4.  The w_s here come
   ! from a brute-force tuple sum of Raoux Eq. 24 for the 2x2 model (X = 0),
   ! independent of the chain form in diag.F90.
   ! Run on a 128^2 and a 256^2 grid: the grid error of Psi_LP (polygon,
   ! interpolation, central differences -- all O(dk^2)) must CONVERGE at second
   ! order, error ratio ~4; a fixed threshold on one grid would only measure how
   ! small the smallest orbit is.  psi_ediff is printed as a diagnostic only: its
   ! second energy difference is noise-limited where Psi is a cancellation.
   call MIO_Print('','semicl')
   call MIO_Print('  test 11: second-order term, massive Dirac cone (must vanish)','semicl')
   deallocate(Wt, Wf, gkx, gky, bx, by)
   dlt = 0.02_dp
   call selftest11_grid(NG/2, rcell, hv, dlt, b, l11c, p11c, e11c, h11c, lp11c, lx11c, t11c)
   call selftest11_grid(NG,   rcell, hv, dlt, b, l11f, p11f, e11f, h11f, lp11f, lx11f, t11f)
   call MIO_Print('    (256^2 grid)  Phi_0        -Phi_1''       Phi_2''''       psi/scale'// &
        '   ediff/scale  Psi_LP       LP err 128^2  LP err 256^2  ratio','semicl')
   worst = 0.0_dp ; worstb = 0.0_dp ; nfail = 0
   do itest = 1, 5
      sc   = maxval(abs(t11f(:,itest)))
      rel  = abs(lp11f(itest) - lx11f(itest))/lx11f(itest)
      relc = abs(lp11c(itest) - lx11c(itest))/lx11c(itest)
      worst  = max(worst, abs(p11f(itest))/sc, abs(p11c(itest))/maxval(abs(t11c(:,itest))))
      worstb = max(worstb, rel)
      if (relc < 2.5_dp*rel) nfail = nfail + 1
      write(*,'(a,f8.4,3(2x,es12.5),2(2x,es10.3),2x,es12.5,2(2x,es10.3),2x,f6.2)') &
           '    ', l11f(itest), t11f(:,itest), p11f(itest)/sc, e11f(itest)/sc, lp11f(itest), &
           relc, rel, relc/max(rel,tiny(1.0_dp))
   end do
   call MIO_Print('    worst |Psi|/(largest term), both grids : '//trim(num2str(worst,6)),'semicl')
   call MIO_Print('    worst Landau-Peierls error vs closed form (256^2): '// &
        trim(num2str(worstb,6)),'semicl')
   if (worst < 0.01_dp .and. worstb < 0.05_dp .and. nfail == 0) then
      call MIO_Print('    PASS (Psi cancels to < 1 % of its terms on both grids; the '// &
           'single-band part does not, matches its closed form and converges at '// &
           'second order)','semicl')
   else
      call MIO_Print('    FAIL ('//trim(num2str(nfail))//' energies converge slower '// &
           'than ratio 2.5)','semicl')
   end if

   call MIO_Print('=================================================','semicl')
   deallocate(Eg, cx, cy, xs, ys)

#ifdef DEBUG
   call MIO_Debug('SemiclSelfTest',1)
#endif /* DEBUG */

end subroutine SemiclSelfTest

!-----------------------------------------------------------------------
! One orbit of the analytic self-test models.  Traces the contour at `lev`,
! keeps the closed non-winding one (these models have exactly one), and
! returns its signed area together with the per-orbit Berry phase and the
! orbital-moment offset.  Factored out so that tests 6-8 exercise EXACTLY the
! call sequence SemiclContoursAreas uses on real data -- a self-test that
! reimplements the pipeline tests nothing.
!-----------------------------------------------------------------------
subroutine selftest_orbit(Eg, Fg, Mg, gkx, gky, ng, rcell, lev, &
                          area, phase, delm, dAdE, ok)

   integer,  intent(in)  :: ng
   real(dp), intent(in)  :: Eg(ng,ng), Fg(ng,ng), Mg(ng,ng)
   real(dp), intent(in)  :: gkx(ng,ng), gky(ng,ng), rcell(3,3), lev
   real(dp), intent(out) :: area, phase, delm, dAdE
   logical,  intent(out) :: ok

   real(dp), allocatable :: bx(:), by(:)
   integer  :: nweak, npb
   real(dp) :: A_BZ

   phase = 0.0_dp ; delm = 0.0_dp ; dAdE = 0.0_dp
   ok = .false. ; nweak = 0
   A_BZ = abs(rcell(1,1)*rcell(2,2) - rcell(1,2)*rcell(2,1))

   allocate(bx(ng*ng), by(ng*ng))
   call pick_orbit(Eg, ng, rcell, lev, bx, by, npb, area)

   if (npb > 0) then
      call SemiclOrbitBerry(Eg, Fg, ng, ng, lev, area, A_BZ, phase, ok)
      call SemiclOrbitMoment(Mg, gkx, gky, ng, ng, bx, by, npb, rcell, &
                             sign(1.0_dp, area), delm, dAdE, nweak)
   end if

   deallocate(bx, by)

end subroutine selftest_orbit

!-----------------------------------------------------------------------
! The largest closed, non-winding contour at lev (the self-test models have
! exactly one), traced and unwrapped exactly as SemiclContoursAreas does.
! npb = 0 when there is none.
!-----------------------------------------------------------------------
subroutine pick_orbit(Eg, ng, rcell, lev, bx, by, npb, area)

   integer,  intent(in)  :: ng
   real(dp), intent(in)  :: Eg(ng,ng), rcell(3,3), lev
   real(dp), intent(out) :: bx(ng*ng), by(ng*ng), area
   integer,  intent(out) :: npb

   real(dp), allocatable :: cx(:), cy(:), xs(:), ys(:)
   integer  :: cstart(MAXCONT), cend(MAXCONT), nc, c, np, wx, wy
   logical  :: closed(MAXCONT), wound
   real(dp) :: ar, abest

   area = 0.0_dp
   allocate(cx(ng*ng), cy(ng*ng), xs(ng*ng), ys(ng*ng))
   call SemiclMarchingSquares(Eg, ng, ng, lev, cx, cy, cstart, cend, nc, closed)

   abest = -1.0_dp ; npb = 0
   do c = 1, nc
      np = cend(c) - cstart(c) + 1
      if (np < 3) cycle
      xs(1:np) = cx(cstart(c):cend(c))
      ys(1:np) = cy(cstart(c):cend(c))
      call SemiclUnwrap(xs, ys, np, wound, wx, wy)
      if (wound .or. .not. closed(c)) cycle
      call SemiclAreaCart(xs, ys, np, rcell, ar)
      if (abs(ar) > abest) then
         abest = abs(ar) ; area = ar ; npb = np
         bx(1:np) = xs(1:np) ; by(1:np) = ys(1:np)
      end if
   end do

   deallocate(cx, cy, xs, ys)

end subroutine pick_orbit

!-----------------------------------------------------------------------
! Self-test 11 input: the second-order weights w_0..w_2 of the UPPER band of
! h = hv (kx sx + ky sy) + D sz by brute-force enumeration of the index tuples
! of Raoux Eq. 24 (X = 0 for a linear model, so only the Q terms), and the
! single-band weight w_LP = det d2E/dk2 = hv^4 D^2/E^4.  Deliberately NOT the
! chain form of diag.F90: this is the independent reference.
!-----------------------------------------------------------------------
subroutine st_dirac_weights(kx, ky, hv, dl, w, wlp)

   use constants, only : cmplx_i
   real(dp), intent(in)  :: kx, ky, hv, dl
   real(dp), intent(out) :: w(0:2), wlp

   real(dp), parameter :: FACT(0:2) = (/1.0_dp, 1.0_dp, 2.0_dp/)
   complex(dp) :: U(2,2), sx(2,2), sy(2,2), Vx(2,2), Vy(2,2), wc(0:2), q
   real(dp) :: ep, th, ph, uu, g, p1, p2, hs(0:2)
   integer  :: n1, n2, n3, n4, idx(4), r, s, j, no

   ep = sqrt((hv*kx)**2 + (hv*ky)**2 + dl*dl)
   th = acos(dl/ep) ; ph = atan2(ky, kx)
   ! columns: lower (1) and upper (2) eigenvectors of [[D, hv k-], [hv k+, -D]]
   U(1,2) = cos(0.5_dp*th)                 ; U(2,2) = exp( cmplx_i*ph)*sin(0.5_dp*th)
   U(1,1) = -exp(-cmplx_i*ph)*sin(0.5_dp*th) ; U(2,1) = cos(0.5_dp*th)
   sx = reshape((/(0.0_dp,0.0_dp), (1.0_dp,0.0_dp), (1.0_dp,0.0_dp), (0.0_dp,0.0_dp)/), (/2,2/))*hv
   sy = reshape((/(0.0_dp,0.0_dp), (0.0_dp,1.0_dp), (0.0_dp,-1.0_dp), (0.0_dp,0.0_dp)/), (/2,2/))*hv
   Vx = matmul(conjg(transpose(U)), matmul(sx, U))
   Vy = matmul(conjg(transpose(U)), matmul(sy, U))

   uu = 1.0_dp/(2.0_dp*ep)                 ! 1/(E_a - E_other), a = 2
   wc = (0.0_dp, 0.0_dp)
   do n4 = 1, 2
      do n3 = 1, 2
         do n2 = 1, 2
            do n1 = 1, 2
               idx = (/n1, n2, n3, n4/)
               r = count(idx == 2)
               if (r == 0 .or. r == 4) cycle
               q = Vx(n1,n2)*Vx(n2,n3)*Vy(n3,n4)*Vy(n4,n1) &
                 - Vx(n1,n2)*Vy(n2,n3)*Vx(n3,n4)*Vy(n4,n1)
               no = 4 - r                   ! slots on the other band, each carrying uu
               g  = uu**no
               p1 = real(no,dp)*uu ; p2 = real(no,dp)*uu*uu
               hs(0) = 1.0_dp ; hs(1) = p1 ; hs(2) = 0.5_dp*(p1*hs(1) + p2)
               do s = 0, r-1
                  j = r - 1 - s
                  wc(s) = wc(s) - 4.0_dp*q*real((-1)**j,dp)*g*hs(j)/FACT(s)
               end do
            end do
         end do
      end do
   end do
   w   = real(wc, dp)
   wlp = hv**4*dl*dl/ep**4

end subroutine st_dirac_weights

!-----------------------------------------------------------------------
! Self-test 11 on one ng x ng grid: the massive Dirac cone with the weights of
! st_dirac_weights, at the five orbit radii of tests 4-7.  Per energy: lev,
! psi (co-area), psied (energy differences), psi_LP and its closed form, and
! the three terms Phi_0, -Phi_1', Phi_2'' separately (the scale of the
! cancellation).  Called at two resolutions by SemiclSelfTest.
!-----------------------------------------------------------------------
subroutine selftest11_grid(ng, rcell, hv, dlt, b, lev, psi, psied, psihyb, psilp, lpex, terms)

   use constants, only : pi
   integer,  intent(in)  :: ng
   real(dp), intent(in)  :: rcell(3,3), hv, dlt, b
   real(dp), intent(out) :: lev(5), psi(5), psied(5), psihyb(5), psilp(5), lpex(5), terms(3,5)

   real(dp), allocatable :: Eg(:,:), W(:,:,:), Wf(:,:,:), Wtt(:,:,:), T(:,:)
   real(dp), allocatable :: gkx(:,:), gky(:,:), bx(:), by(:)
   real(dp) :: fx, fy, kx, ky, kr, w3(0:2), wlp, area, hst, L
   real(dp) :: st0, st1, st2
   integer  :: i, j, it, nsg, npb, nwk

   allocate(Eg(ng,ng), W(ng,ng,4), Wf(ng,ng,5), Wtt(ng,ng,3), T(ng,ng), &
            gkx(ng,ng), gky(ng,ng), bx(ng*ng), by(ng*ng))
   do j = 1, ng
      do i = 1, ng
         fx = real(i-1,dp)/real(ng,dp) - 0.5_dp
         fy = real(j-1,dp)/real(ng,dp) - 0.5_dp
         kx = fx*rcell(1,1) + fy*rcell(1,2)
         ky = fx*rcell(2,1) + fy*rcell(2,2)
         Eg(i,j) = sqrt(hv*hv*(kx*kx + ky*ky) + dlt*dlt)
         call st_dirac_weights(kx, ky, hv, dlt, w3, wlp)
         W(i,j,1:3) = w3 ; W(i,j,4) = wlp
      end do
   end do
   nsg = 0
   call SemiclGradient(Eg, ng, ng, rcell, gkx, gky)
   call SemiclSuscFields(W, Eg, gkx, gky, ng, ng, rcell, psiHybOn, psiTermsOn, Wf, nsg)
   ! the three terms separately, only to give the cancellation its scale
   Wtt(:,:,1) = W(:,:,1)
   call semicl_coarea_div(W(:,:,2), Eg, gkx, gky, ng, ng, rcell, Wtt(:,:,2), nsg)
   Wtt(:,:,2) = -Wtt(:,:,2)
   call semicl_coarea_div(W(:,:,3), Eg, gkx, gky, ng, ng, rcell, T, nsg)
   call semicl_coarea_div(T, Eg, gkx, gky, ng, ng, rcell, Wtt(:,:,3), nsg)

   do it = 1, 5
      kr = (0.05_dp + 0.05_dp*it) * (b/2.0_dp)
      lev(it)  = sqrt(hv*hv*kr*kr + dlt*dlt)
      lpex(it) = 6.0_dp*pi*hv*hv*dlt*dlt/lev(it)**4
      call pick_orbit(Eg, ng, rcell, lev(it), bx, by, npb, area)
      if (npb == 0) then                       ! a bug: make the gate fail
         psi(it) = huge(1.0_dp) ; psied(it) = 0.0_dp ; psilp(it) = 0.0_dp
         psihyb(it) = 0.0_dp
         terms(:,it) = 1.0_dp
         cycle
      end if
      hst = 0.25_dp*(b/real(ng,dp))*hv*hv*kr/lev(it)
      nwk = 0
      call SemiclOrbitSusc(Eg, W, Wf, gkx, gky, ng, ng, bx, by, npb, area, rcell, &
                           lev(it), hst, psiHybOn, psiTermsOn, psi(it), psied(it), psihyb(it), &
                           psilp(it), st0, st1, st2, nwk)
      call orbit_line_int(Wtt, 3, gkx, gky, ng, ng, bx, by, npb, rcell, terms(:,it), L, nwk)
   end do

   deallocate(Eg, W, Wf, Wtt, T, gkx, gky, bx, by)

end subroutine selftest11_grid

!=======================================================================
! Read a 3-D TAPW bands file written by DiagBands with
! Diag.Calculate3DTAPWBands .true.
!
!   line 1 : Efermi
!   line 2 : 0.0  d
!   line 3 : Emin-2  Emax+2
!   line 4 : M  nspin  nk
!   then nk records:  kx ky kz  E(1..M*nspin)     (energies already in eV)
!
! The Monkhorst-Pack loop in diag.F90 is  do j=1,ny ; do i=1,nx , so record
! ip = (j-1)*nx + i and the grid unpacks with i fastest.
!=======================================================================
subroutine SemiclReadBands3D(fname, nx, ny, Egrid, nb, ok)

   character(*), intent(in)  :: fname
   integer,      intent(in)  :: nx, ny
   real(dp), allocatable, intent(out) :: Egrid(:,:,:)
   integer,      intent(out) :: nb
   logical,      intent(out) :: ok

   integer  :: u, ios, M, nspin, nk, ip, i, j, b
   real(dp) :: Ef, d0, d1, e0, e1, kx, ky, kz
   real(dp), allocatable :: row(:)

   ok = .false. ; nb = 0
   u = 771
   open(u, file=trim(fname), status='old', action='read', iostat=ios)
   if (ios /= 0) then
      call MIO_Print('ERROR: cannot open bands file '//trim(fname),'semicl')
      return
   end if

   read(u,*,iostat=ios) Ef             ; if (ios/=0) goto 900
   read(u,*,iostat=ios) d0, d1         ; if (ios/=0) goto 900
   read(u,*,iostat=ios) e0, e1         ; if (ios/=0) goto 900
   read(u,*,iostat=ios) M, nspin, nk   ; if (ios/=0) goto 900

   if (nk /= nx*ny) then
      call MIO_Print('ERROR: bands file has '//trim(num2str(nk))// &
           ' k-points but Semicl.GridX*GridY = '//trim(num2str(nx*ny)),'semicl')
      call MIO_Print('  (Semicl.GridX/GridY must match Diag.3DBandsGridX/Y)','semicl')
      close(u) ; return
   end if

   nb = M*nspin
   allocate(Egrid(nx,ny,nb), row(nb))

   do ip = 1, nk
      read(u,*,iostat=ios) kx, ky, kz, (row(b), b=1,nb)
      if (ios /= 0) then
         call MIO_Print('ERROR: truncated bands file at record '// &
              trim(num2str(ip)),'semicl')
         close(u) ; return
      end if
      j = (ip-1)/nx + 1
      i = ip - (j-1)*nx
      Egrid(i,j,:) = row(:)
   end do
   close(u)
   deallocate(row)
   ok = .true.
   call MIO_Print('Read '//trim(num2str(nk))//' k-points x '// &
        trim(num2str(nb))//' bands from '//trim(fname),'semicl')
   return

900 continue
   call MIO_Print('ERROR: malformed header in '//trim(fname),'semicl')
   close(u)

end subroutine SemiclReadBands3D

!=======================================================================
! Reads <prefix>.BerryFlux, the per-plaquette Fukui-Hatsugai field strength
! written by diag.F90 on the same k-grid as generate.bands:
!
!   # grid nx ny
!   # bands b1 b2
!   iband  i  j  flux[rad]
!
! The sum of F over any set of plaquettes is the Berry phase of the boundary
! of that region, so an orbit's Berry phase is a plain sum over its interior.
!=======================================================================
subroutine SemiclReadBerryFlux(fname, nx, ny, Bflux, b1, b2, ok)

   character(*), intent(in)  :: fname
   integer,      intent(in)  :: nx, ny
   real(dp), allocatable, intent(out) :: Bflux(:,:,:)
   integer,      intent(out) :: b1, b2
   logical,      intent(out) :: ok

   integer  :: u, ios, gx, gy, ib, i, j, nread
   real(dp) :: F
   character(len=400) :: line

   ok = .false. ; b1 = 0 ; b2 = -1
   u = 773
   open(u, file=trim(fname), status='old', action='read', iostat=ios)
   if (ios /= 0) then
      call MIO_Print('  Semicl.BerryPhase: cannot open '//trim(fname),'semicl')
      return
   end if

   gx = 0 ; gy = 0
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) /= '#') exit
      if (index(line,'# grid ') == 1)  read(line(8:),*,iostat=ios) gx, gy
      if (index(line,'# bands ') == 1) read(line(9:),*,iostat=ios) b1, b2
   end do

   if (gx /= nx .or. gy /= ny) then
      call MIO_Print('  Semicl.BerryPhase: '//trim(fname)//' is '//trim(num2str(gx))// &
           'x'//trim(num2str(gy))//' but Semicl.GridX/Y = '//trim(num2str(nx))// &
           'x'//trim(num2str(ny)),'semicl')
      close(u) ; return
   end if
   if (b2 < b1) then
      call MIO_Print('  Semicl.BerryPhase: no band window in '//trim(fname),'semicl')
      close(u) ; return
   end if

   allocate(Bflux(nx, ny, b2-b1+1))
   Bflux = 0.0_dp
   rewind(u)
   nread = 0
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
      read(line,*,iostat=ios) ib, i, j, F
      if (ios /= 0) cycle
      if (ib < b1 .or. ib > b2) cycle
      if (i < 1 .or. i > nx .or. j < 1 .or. j > ny) cycle
      Bflux(i,j,ib-b1+1) = F
      nread = nread + 1
   end do
   close(u)

   if (nread /= nx*ny*(b2-b1+1)) then
      call MIO_Print('  Semicl.BerryPhase: read '//trim(num2str(nread))//' of '// &
           trim(num2str(nx*ny*(b2-b1+1)))//' plaquettes from '//trim(fname),'semicl')
      deallocate(Bflux) ; return
   end if

   call MIO_Print('  read Berry fluxes for bands '//trim(num2str(b1))//'..'// &
        trim(num2str(b2))//' from '//trim(fname),'semicl')
   ok = .true.

end subroutine SemiclReadBerryFlux

!=======================================================================
! Connected components of a boolean mask on the k-torus (4-connectivity,
! periodic in both directions).  Returns lab(i,j) = 1..ncomp for masked
! points, 0 elsewhere.  Used to separate several pockets living on the same
! sheet at the same energy: Onsager quantises per orbit, so their Berry
! phases must not be summed together either.
!=======================================================================
subroutine SemiclLabelComponents(mask, nx, ny, lab, ncomp)

   integer, intent(in)  :: nx, ny
   logical, intent(in)  :: mask(nx,ny)
   integer, intent(out) :: lab(nx,ny), ncomp

   integer :: i, j, head, tail, ci, cj, di, ni, nj, d
   integer, allocatable :: qi(:), qj(:)
   integer :: dx(4), dy(4)

   dx = [ 1, -1,  0,  0 ]
   dy = [ 0,  0,  1, -1 ]

   allocate(qi(nx*ny), qj(nx*ny))
   lab = 0 ; ncomp = 0

   do j = 1, ny
      do i = 1, nx
         if (.not. mask(i,j) .or. lab(i,j) /= 0) cycle
         ncomp = ncomp + 1
         head = 1 ; tail = 1
         qi(1) = i ; qj(1) = j ; lab(i,j) = ncomp
         do while (head <= tail)
            ci = qi(head) ; cj = qj(head) ; head = head + 1
            do d = 1, 4
               ni = modulo(ci + dx(d) - 1, nx) + 1
               nj = modulo(cj + dy(d) - 1, ny) + 1
               if (mask(ni,nj) .and. lab(ni,nj) == 0) then
                  lab(ni,nj) = ncomp
                  tail = tail + 1
                  qi(tail) = ni ; qj(tail) = nj
               end if
            end do
         end do
      end do
   end do

   deallocate(qi, qj)

end subroutine SemiclLabelComponents

!=======================================================================
! Berry phase of one orbit: sum the plaquette fluxes over the interior the
! contour encloses.  The interior is the connected component of
! {E < level} (electron pocket, signed area > 0) or of {E > level}
! (hole pocket, signed area < 0) whose grid-counted area is closest to the
! contour's own area -- the same rank-by-area matching E1 uses to follow an
! orbit between energies.
!
! A plaquette is taken as inside when ALL FOUR of its corners are: the flux
! of the boundary ring is then missing, which is negligible for the Dirac
! orbits (their curvature is concentrated at the band touching, far from the
! boundary) and small in absolute terms for smooth parabolic pockets, where
! the phase itself is near zero.  matched = .false. means no component could
! be identified and the caller must not write a gamma.
!=======================================================================
subroutine SemiclOrbitBerry(Esheet, Bfl, nx, ny, level, area, A_BZ, phase, matched)

   integer,  intent(in)  :: nx, ny
   real(dp), intent(in)  :: Esheet(nx,ny), Bfl(nx,ny), level, area, A_BZ
   real(dp), intent(out) :: phase
   logical,  intent(out) :: matched

   logical, allocatable :: mask(:,:)
   integer, allocatable :: lab(:,:)
   integer :: ncomp, i, j, c, cbest, ip, jp, nin
   real(dp) :: acell, adiff, abest, acomp
   integer, allocatable :: ncell(:)

   phase = 0.0_dp ; matched = .false.
   allocate(mask(nx,ny), lab(nx,ny))

   if (area >= 0.0_dp) then
      mask = Esheet < level          ! electron pocket
   else
      mask = Esheet > level          ! hole pocket: the CCW interior is 'above'
   end if

   call SemiclLabelComponents(mask, nx, ny, lab, ncomp)
   if (ncomp < 1) then
      deallocate(mask, lab) ; return
   end if

   allocate(ncell(ncomp))
   ncell = 0
   do j = 1, ny
      do i = 1, nx
         if (lab(i,j) > 0) ncell(lab(i,j)) = ncell(lab(i,j)) + 1
      end do
   end do

   acell = A_BZ/real(nx*ny,dp)
   cbest = 0 ; abest = huge(1.0_dp)
   do c = 1, ncomp
      acomp = real(ncell(c),dp)*acell
      adiff = abs(acomp - abs(area))
      if (adiff < abest) then
         abest = adiff ; cbest = c
      end if
   end do

   ! Reject a match that is nowhere near the contour area: with several
   ! pockets of similar size the rank matching can be wrong, and a silently
   ! mismatched gamma is worse than none.  The tolerance has to allow the
   ! genuine O(one cell) difference between a polygon area and a cell count,
   ! which dominates for the small pockets near a fan origin.
   if (cbest == 0 .or. abest > max(0.25_dp*abs(area), 3.0_dp*acell)) then
      deallocate(mask, lab, ncell) ; return
   end if

   ! A plaquette straddling the contour is weighted by the fraction of its
   ! corners inside.  Taking only fully-enclosed plaquettes drops the whole
   ! boundary ring and cost ~5 % of the phase on the smallest self-test orbit;
   ! the fraction rule is the symmetric choice and removes that bias without
   ! pretending to sub-cell accuracy.
   do j = 1, ny
      do i = 1, nx
         ip = modulo(i, nx) + 1
         jp = modulo(j, ny) + 1
         nin = 0
         if (lab(i ,j ) == cbest) nin = nin + 1
         if (lab(ip,j ) == cbest) nin = nin + 1
         if (lab(i ,jp) == cbest) nin = nin + 1
         if (lab(ip,jp) == cbest) nin = nin + 1
         if (nin > 0) phase = phase + Bfl(i,j)*0.25_dp*real(nin,dp)
      end do
   end do

   matched = .true.
   deallocate(mask, lab, ncell)

end subroutine SemiclOrbitBerry

!=======================================================================
! Read the per-k-point orbital magnetic moment written by diag when
! Diag.OrbMoment .true.  (L1a).
!
!   # grid nx ny
!   # bands b1 b2
!   iband  i  j  M[eV.Ang^2]  Omega[Ang^2]
!
! M is the integrand of the orbital magnetic moment,
!
!   M_n(k) = Im sum_{m/=n} V^x_nm V^y_mn / (E_n - E_m)          [eV Ang^2]
!
! so that the field-modified semiclassical band energy is
! E~(k) = E(k) + (e/hbar)*M(k)*B.  Values live at the GRID POINTS, i.e. at the
! same k as generate.bands -- unlike the Berry fluxes, which live on the
! plaquettes between them.
!
! Omega is the Kubo Berry curvature built from the SAME velocity matrix
! elements.  semicl does not use it: it is written so that its sign can be
! compared with the Fukui-Hatsugai link fluxes.  The relative sign of the
! curvature and the moment is fixed by construction (same V, same loop), but
! the sign of the link product against the Kubo formula is a convention that
! only a measurement settles.
!=======================================================================
subroutine SemiclReadOrbMoment(fname, nx, ny, Mz, b1, b2, ok)

   character(*), intent(in)  :: fname
   integer,      intent(in)  :: nx, ny
   real(dp), allocatable, intent(out) :: Mz(:,:,:)
   integer,      intent(out) :: b1, b2
   logical,      intent(out) :: ok

   integer  :: u, ios, gx, gy, ib, i, j, nread
   real(dp) :: Mval, Oval
   character(len=400) :: line

   ok = .false. ; b1 = 0 ; b2 = -1
   u = 776
   open(u, file=trim(fname), status='old', action='read', iostat=ios)
   if (ios /= 0) then
      call MIO_Print('  Semicl.OrbitalMoment: cannot open '//trim(fname),'semicl')
      return
   end if

   gx = 0 ; gy = 0
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) /= '#') exit
      if (index(line,'# grid ') == 1)  read(line(8:),*,iostat=ios) gx, gy
      if (index(line,'# bands ') == 1) read(line(9:),*,iostat=ios) b1, b2
   end do

   if (gx /= nx .or. gy /= ny) then
      call MIO_Print('  Semicl.OrbitalMoment: '//trim(fname)//' is '//trim(num2str(gx))// &
           'x'//trim(num2str(gy))//' but Semicl.GridX/Y = '//trim(num2str(nx))// &
           'x'//trim(num2str(ny)),'semicl')
      close(u) ; return
   end if
   if (b2 < b1) then
      call MIO_Print('  Semicl.OrbitalMoment: no band window in '//trim(fname),'semicl')
      close(u) ; return
   end if

   allocate(Mz(nx, ny, b2-b1+1))
   Mz = 0.0_dp
   rewind(u)
   nread = 0
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
      read(line,*,iostat=ios) ib, i, j, Mval, Oval
      if (ios /= 0) cycle
      if (ib < b1 .or. ib > b2) cycle
      if (i < 1 .or. i > nx .or. j < 1 .or. j > ny) cycle
      Mz(i,j,ib-b1+1) = Mval
      nread = nread + 1
   end do
   close(u)

   if (nread /= nx*ny*(b2-b1+1)) then
      call MIO_Print('  Semicl.OrbitalMoment: read '//trim(num2str(nread))//' of '// &
           trim(num2str(nx*ny*(b2-b1+1)))//' points from '//trim(fname),'semicl')
      deallocate(Mz) ; return
   end if

   call MIO_Print('  read orbital moments for bands '//trim(num2str(b1))//'..'// &
        trim(num2str(b2))//' from '//trim(fname),'semicl')
   ok = .true.

end subroutine SemiclReadOrbMoment

!=======================================================================
! Read <prefix>.OrbSusc (Diag.OrbSusc, L1c/L2): per grid point and band the
! weights w_0, w_1, w_2 of the second-order term and the band curvature
! exx, eyy, exy.  W(:,:,4,b) = exx*eyy - exy^2 is the single-band
! (Landau-Peierls) weight, kept as the control.
!=======================================================================
subroutine SemiclReadOrbSusc(fname, nx, ny, W, b1, b2, ok)

   character(*), intent(in)  :: fname
   integer,      intent(in)  :: nx, ny
   real(dp), allocatable, intent(out) :: W(:,:,:,:)
   integer,      intent(out) :: b1, b2
   logical,      intent(out) :: ok

   integer  :: u, ios, gx, gy, ib, i, j, nread
   real(dp) :: w3(3), e3(3)
   character(len=400) :: line

   ok = .false. ; b1 = 0 ; b2 = -1
   u = 777
   open(u, file=trim(fname), status='old', action='read', iostat=ios)
   if (ios /= 0) then
      call MIO_Print('  Semicl.Susceptibility: cannot open '//trim(fname),'semicl')
      return
   end if

   gx = 0 ; gy = 0
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) /= '#') exit
      if (index(line,'# grid ') == 1)  read(line(8:),*,iostat=ios) gx, gy
      if (index(line,'# bands ') == 1) read(line(9:),*,iostat=ios) b1, b2
   end do

   if (gx /= nx .or. gy /= ny) then
      call MIO_Print('  Semicl.Susceptibility: '//trim(fname)//' is '//trim(num2str(gx))// &
           'x'//trim(num2str(gy))//' but Semicl.GridX/Y = '//trim(num2str(nx))// &
           'x'//trim(num2str(ny)),'semicl')
      close(u) ; return
   end if
   if (b2 < b1) then
      call MIO_Print('  Semicl.Susceptibility: no band window in '//trim(fname),'semicl')
      close(u) ; return
   end if

   allocate(W(nx, ny, 4, b2-b1+1))
   W = 0.0_dp
   rewind(u)
   nread = 0
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
      read(line,*,iostat=ios) ib, i, j, w3, e3
      if (ios /= 0) cycle
      if (ib < b1 .or. ib > b2) cycle
      if (i < 1 .or. i > nx .or. j < 1 .or. j > ny) cycle
      W(i,j,1:3,ib-b1+1) = w3
      W(i,j,4,ib-b1+1)   = e3(1)*e3(2) - e3(3)*e3(3)
      nread = nread + 1
   end do
   close(u)

   if (nread /= nx*ny*(b2-b1+1)) then
      call MIO_Print('  Semicl.Susceptibility: read '//trim(num2str(nread))//' of '// &
           trim(num2str(nx*ny*(b2-b1+1)))//' points from '//trim(fname),'semicl')
      deallocate(W) ; return
   end if

   call MIO_Print('  read second-order weights for bands '//trim(num2str(b1))//'..'// &
        trim(num2str(b2))//' from '//trim(fname),'semicl')
   ok = .true.

end subroutine SemiclReadOrbSusc

!=======================================================================
! Cartesian k-space gradient of one energy sheet on the periodic grid.
!
! Central differences give d/df in FRACTIONAL coordinates; the Cartesian
! gradient follows from  df_i/dk = (B^-1)_i  with B the 2x2 matrix whose ROWS
! are the reciprocal vectors b1, b2 (the columns of rcell), since
! dE/df_i = b_i . grad_k E.
!=======================================================================
subroutine SemiclGradient(Esheet, nx, ny, rcell, gkx, gky)

   integer,  intent(in)  :: nx, ny
   real(dp), intent(in)  :: Esheet(nx,ny), rcell(3,3)
   real(dp), intent(out) :: gkx(nx,ny), gky(nx,ny)

   integer  :: i, j
   real(dp) :: g1, g2, det, bi(2,2)

   ! B = [[b1x,b1y],[b2x,b2y]]  ->  inverse
   det = rcell(1,1)*rcell(2,2) - rcell(2,1)*rcell(1,2)
   if (abs(det) < tiny(1.0_dp)) then
      gkx = 0.0_dp ; gky = 0.0_dp ; return
   end if
   bi(1,1) =  rcell(2,2)/det ; bi(1,2) = -rcell(2,1)/det
   bi(2,1) = -rcell(1,2)/det ; bi(2,2) =  rcell(1,1)/det

   do j = 1, ny
      do i = 1, nx
         g1 = 0.5_dp*real(nx,dp)*(Esheet(iwrap(i+1,nx),j) - Esheet(iwrap(i-1,nx),j))
         g2 = 0.5_dp*real(ny,dp)*(Esheet(i,iwrap(j+1,ny)) - Esheet(i,iwrap(j-1,ny)))
         gkx(i,j) = bi(1,1)*g1 + bi(1,2)*g2
         gky(i,j) = bi(2,1)*g1 + bi(2,2)*g2
      end do
   end do

end subroutine SemiclGradient

!-----------------------------------------------------------------------
! Periodic bilinear interpolation of a grid quantity at fractional (fx,fy),
! with grid point (i,j) sitting at ( (i-1)/nx , (j-1)/ny ).
!-----------------------------------------------------------------------
pure real(dp) function bilin(F, nx, ny, fx, fy)
   integer,  intent(in) :: nx, ny
   real(dp), intent(in) :: F(nx,ny), fx, fy
   integer  :: i0, j0, i1, j1
   real(dp) :: tx, ty, ax, ay
   tx = (fx - floor(fx))*real(nx,dp)
   ty = (fy - floor(fy))*real(ny,dp)
   i0 = int(tx) ; ax = tx - real(i0,dp)
   j0 = int(ty) ; ay = ty - real(j0,dp)
   i0 = iwrap(i0+1, nx) ; i1 = iwrap(i0+1, nx)
   j0 = iwrap(j0+1, ny) ; j1 = iwrap(j0+1, ny)
   bilin = (1.0_dp-ax)*(1.0_dp-ay)*F(i0,j0) + ax*(1.0_dp-ay)*F(i1,j0) &
         + (1.0_dp-ax)*ay        *F(i0,j1) + ax*ay        *F(i1,j1)
end function bilin

!=======================================================================
! Per-orbit orbital-moment correction to the Onsager offset (L1b).
!
! At first order in B the semiclassical band energy carries the orbital
! magnetic moment,
!
!       E~(k) = E(k) + (e/hbar) M(k) B ,
!
! so the orbit at total energy E is the contour of E~, not of E, and its area
! is displaced by the normal shift of every boundary point:
!
!       A~(E) = A(E) - B (e/hbar) \oint M dl / |grad_k E| .
!
! Onsager, SEMICL_C*A~ = (N+gamma)*B, then keeps its old form with one extra
! per-orbit offset:
!
!       B = SEMICL_C * |A| / (N + gamma + delta_m) ,
!       delta_m = SEMICL_C*(e/hbar) * sgn * \oint M dl/|grad E|
!               = (sgn/2pi) * \oint M dl/|grad E| ,
!
! because SEMICL_C*(e/hbar) = [hbar/(2 pi e)]*[e/hbar] = 1/(2 pi) EXACTLY.
! The two unit conversions cancel, delta_m is a pure number, and no new
! physical constant enters this module -- if one ever seems necessary,
! something upstream is wrong (the same rule as for SEMICL_C itself).
!
! sgn = +1 for an electron orbit, -1 for a hole orbit: the quantised area is
! |A|, and d|A|/dE = sgn * \oint dl/|grad E|.
!
! For a massive Dirac cone this returns delta_m = Delta/(2|E-E_D|), which
! cancels the Berry-phase gamma exactly and leaves N+gamma+delta_m an INTEGER
! -- i.e. the exact Landau levels E_n = sqrt(Delta^2 + 2n(hbar v)^2/l_B^2).
! That is the analytic gate in SemiclSelfTest (tests 6 and 7), and it is the
! whole point of the correction: the Berry phase alone is NOT the complete
! first-order term.
!
! The same loop returns dA/dE as a line integral.  E1 gets dA/dE by central
! differences between neighbouring energies and has to avoid Lifshitz
! transitions; this one is local to a single contour.  They are independent,
! so a disagreement is a diagnostic rather than a duplication.
!=======================================================================
subroutine SemiclOrbitMoment(Mz, gkx, gky, nx, ny, xs, ys, np, rcell, sgn, &
                             delta_m, dAdE, nweak)

   use constants, only : pi

   integer,  intent(in)    :: nx, ny, np
   real(dp), intent(in)    :: Mz(nx,ny), gkx(nx,ny), gky(nx,ny)
   real(dp), intent(in)    :: xs(np), ys(np), rcell(3,3), sgn
   real(dp), intent(out)   :: delta_m, dAdE
   integer,  intent(inout) :: nweak

   real(dp) :: J(1), L

   call orbit_line_int(Mz, 1, gkx, gky, nx, ny, xs, ys, np, rcell, J, L, nweak)

   delta_m = sgn*J(1)/(2.0_dp*pi)
   dAdE    = sgn*L

end subroutine SemiclOrbitMoment

!-----------------------------------------------------------------------
! closed_int F_f dl/|grad E| along one unwrapped contour, for nf grid fields at
! once, and L = closed_int dl/|grad E| = |dA/dE|.  Fields are interpolated
! bilinearly at the segment midpoints.  Shared by the orbital moment (L1b) and
! the second-order term (L1c/L2).
!-----------------------------------------------------------------------
subroutine orbit_line_int(F, nf, gkx, gky, nx, ny, xs, ys, np, rcell, J, L, nweak)

   integer,  intent(in)    :: nf, nx, ny, np
   real(dp), intent(in)    :: F(nx,ny,nf), gkx(nx,ny), gky(nx,ny)
   real(dp), intent(in)    :: xs(np), ys(np), rcell(3,3)
   real(dp), intent(out)   :: J(nf), L
   integer,  intent(inout) :: nweak

   integer  :: p, q, ifl
   real(dp) :: fx, fy, dkx, dky, dl, gx, gy, gnorm, w

   J = 0.0_dp ; L = 0.0_dp
   do p = 1, np
      q = p + 1
      if (q > np) q = 1                     ! the polygon closes back on itself
      dkx = rcell(1,1)*(xs(q)-xs(p)) + rcell(1,2)*(ys(q)-ys(p))
      dky = rcell(2,1)*(xs(q)-xs(p)) + rcell(2,2)*(ys(q)-ys(p))
      dl  = sqrt(dkx*dkx + dky*dky)
      if (dl <= 0.0_dp) cycle
      fx = 0.5_dp*(xs(p)+xs(q))
      fy = 0.5_dp*(ys(p)+ys(q))
      gx = bilin(gkx, nx, ny, fx, fy)
      gy = bilin(gky, nx, ny, fx, fy)
      gnorm = sqrt(gx*gx + gy*gy)
      ! |grad E| -> 0 only at a band extremum, where the orbit has shrunk to a
      ! point; a segment sitting there would contribute an unbounded weight to
      ! both integrals.  Drop it and count, rather than clamping silently.
      if (gnorm <= 1.0e-12_dp) then
         nweak = nweak + 1
         cycle
      end if
      w = dl/gnorm
      L = L + w
      do ifl = 1, nf
         J(ifl) = J(ifl) + w*bilin(F(:,:,ifl), nx, ny, fx, fy)
      end do
   end do

end subroutine orbit_line_int

!-----------------------------------------------------------------------
! Co-area divergence  D(F) = div( F grad E / |grad E|^2 )  on the periodic grid.
!
! Integrating D(F) over {E_k < E} and applying Gauss, the flux of
! F grad E/|grad E|^2 through the contour is closed_int F dl/|grad E|, so
!       d/dE closed_int F dl/|grad E|  =  closed_int D(F) dl/|grad E| .
! An energy derivative of a per-orbit line integral is therefore ONE line
! integral of a grid field.  The identity is local to each contour, so it holds
! for electron and hole orbits alike.  Central differences in the
! SemiclGradient convention; a grid point with |grad E| <= 1e-12 (a critical
! point) gets zero flux and is counted.
!-----------------------------------------------------------------------
subroutine semicl_coarea_div(F, Esheet, gkx, gky, nx, ny, rcell, D, nsing)

   integer,  intent(in)    :: nx, ny
   real(dp), intent(in)    :: F(nx,ny), Esheet(nx,ny), gkx(nx,ny), gky(nx,ny), rcell(3,3)
   real(dp), intent(out)   :: D(nx,ny)
   integer,  intent(inout) :: nsing

   real(dp), allocatable :: ax(:,:), ay(:,:)
   integer  :: i, j
   real(dp) :: g2, det, bi(2,2), d1x, d2x, d1y, d2y
   real(dp) :: e11, e22, e12, Exx, Eyy, Exy, dF1, dF2, dFx, dFy, quad, lap
   real(dp) :: rn1, rn2

   D = 0.0_dp
   det = rcell(1,1)*rcell(2,2) - rcell(2,1)*rcell(1,2)
   if (abs(det) < tiny(1.0_dp)) return
   bi(1,1) =  rcell(2,2)/det ; bi(1,2) = -rcell(2,1)/det
   bi(2,1) = -rcell(1,2)/det ; bi(2,2) =  rcell(1,1)/det

   if (coareaAn) then
      ! ANALYTIC EXPANSION.  Writing the divergence out,
      !   D(F) = (gradF.gradE)/|gradE|^2
      !          + F*[lap(E)*|gradE|^2 - 2*gradE^T H gradE]/|gradE|^4 ,
      ! the cancellation that keeps D finite happens inside ONE numerator built
      ! from smooth derivatives of E, instead of between two large differenced
      ! fields of F*gradE/|gradE|^2.  The differenced form loses that
      ! cancellation to roundoff-scale subtraction of big numbers, which is why
      ! its integrated error does not fall with the grid.  D is still singular
      ! at a SADDLE (the continuum limit diverges as 1/r^2 there); what this
      ! buys is an integrable, converging ERROR, which is what the per-orbit
      ! line integral needs.
      rn1 = real(nx,dp) ; rn2 = real(ny,dp)
      do j = 1, ny
         do i = 1, nx
            g2 = gkx(i,j)**2 + gky(i,j)**2
            if (g2 <= 1.0e-24_dp) then
               nsing = nsing + 1 ; cycle          ! D already 0 here
            end if
            ! second derivatives in the FRACTIONAL frame (step 1/n), then
            ! rotated to Cartesian with bi = d f / d k, as the first-derivative
            ! branch below does for the gradient
            e11 = rn1*rn1*(Esheet(iwrap(i+1,nx),j) - 2.0_dp*Esheet(i,j) &
                         + Esheet(iwrap(i-1,nx),j))
            e22 = rn2*rn2*(Esheet(i,iwrap(j+1,ny)) - 2.0_dp*Esheet(i,j) &
                         + Esheet(i,iwrap(j-1,ny)))
            e12 = 0.25_dp*rn1*rn2* &
                 (Esheet(iwrap(i+1,nx),iwrap(j+1,ny)) - Esheet(iwrap(i+1,nx),iwrap(j-1,ny)) &
                - Esheet(iwrap(i-1,nx),iwrap(j+1,ny)) + Esheet(iwrap(i-1,nx),iwrap(j-1,ny)))
            Exx = bi(1,1)*bi(1,1)*e11 + 2.0_dp*bi(1,1)*bi(1,2)*e12 + bi(1,2)*bi(1,2)*e22
            Eyy = bi(2,1)*bi(2,1)*e11 + 2.0_dp*bi(2,1)*bi(2,2)*e12 + bi(2,2)*bi(2,2)*e22
            Exy = bi(1,1)*bi(2,1)*e11 &
                + (bi(1,1)*bi(2,2) + bi(1,2)*bi(2,1))*e12 + bi(1,2)*bi(2,2)*e22
            dF1 = 0.5_dp*rn1*(F(iwrap(i+1,nx),j) - F(iwrap(i-1,nx),j))
            dF2 = 0.5_dp*rn2*(F(i,iwrap(j+1,ny)) - F(i,iwrap(j-1,ny)))
            dFx = bi(1,1)*dF1 + bi(1,2)*dF2
            dFy = bi(2,1)*dF1 + bi(2,2)*dF2
            lap  = Exx + Eyy
            quad = gkx(i,j)*gkx(i,j)*Exx + 2.0_dp*gkx(i,j)*gky(i,j)*Exy &
                 + gky(i,j)*gky(i,j)*Eyy
            D(i,j) = (dFx*gkx(i,j) + dFy*gky(i,j))/g2 &
                   + F(i,j)*(lap*g2 - 2.0_dp*quad)/(g2*g2)
         end do
      end do
      return
   end if

   allocate(ax(nx,ny), ay(nx,ny))
   do j = 1, ny
      do i = 1, nx
         g2 = gkx(i,j)**2 + gky(i,j)**2
         if (g2 <= 1.0e-24_dp) then
            ax(i,j) = 0.0_dp ; ay(i,j) = 0.0_dp
            nsing = nsing + 1
         else
            ax(i,j) = F(i,j)*gkx(i,j)/g2
            ay(i,j) = F(i,j)*gky(i,j)/g2
         end if
      end do
   end do
   do j = 1, ny
      do i = 1, nx
         d1x = 0.5_dp*real(nx,dp)*(ax(iwrap(i+1,nx),j) - ax(iwrap(i-1,nx),j))
         d2x = 0.5_dp*real(ny,dp)*(ax(i,iwrap(j+1,ny)) - ax(i,iwrap(j-1,ny)))
         d1y = 0.5_dp*real(nx,dp)*(ay(iwrap(i+1,nx),j) - ay(iwrap(i-1,nx),j))
         d2y = 0.5_dp*real(ny,dp)*(ay(i,iwrap(j+1,ny)) - ay(i,iwrap(j-1,ny)))
         D(i,j) = bi(1,1)*d1x + bi(1,2)*d2x + bi(2,1)*d1y + bi(2,2)*d2y
      end do
   end do
   deallocate(ax, ay)

end subroutine semicl_coarea_div

!-----------------------------------------------------------------------
! Grid fields whose line integrals ARE the second-order quantities (L1c/L2):
!   Weff(:,:,1) = w_0 - D(w_1 - D(w_2))   ->  Psi    = Phi_0 - Phi_1' + Phi_2''
!   Weff(:,:,2) = -D(w_LP)                ->  Psi_LP = the single-band
!                                             (Landau-Peierls) part, a control
! W(:,:,1:4) = w_0, w_1, w_2, w_LP of one sheet.  Built once per sheet.
!
! With hyb, two more fields for the HYBRID route (Semicl.PsiHybrid).  The
! nested D(D(.)) above is the one construction that cannot converge: D divides
! by |grad E|^2 and then central-differences, so its peak grows as 1/h^2 near a
! critical point while truncation error only falls as h^2; nesting it squares
! that to 1/h^4, and refining the grid moves points CLOSER to the critical
! point.  Measured on an analytic cosine band against an exact D: the single
! operator's worst point grows x4 per doubling, the nested one x16, while away
! from critical points both converge as O(h^2).  The hybrid keeps ONE co-area
! divergence and takes the remaining energy derivative on the per-orbit line
! integral, where it acts on a smooth function of E:
!   Weff(:,:,3) = D(w_1)  ->  Phi_1' directly
!   Weff(:,:,4) = D(w_2)  ->  Phi_2'' as d/dE of ITS line integral
! so the singular order drops from 1/h^4 to 1/h^2.  D(w_2) is already built as
! a temporary for Weff(:,:,1), so only D(w_1) is extra work.
!-----------------------------------------------------------------------
subroutine SemiclSuscFields(W, Esheet, gkx, gky, nx, ny, rcell, hyb, terms, Weff, nsing)

   integer,  intent(in)    :: nx, ny
   real(dp), intent(in)    :: W(nx,ny,4), Esheet(nx,ny), gkx(nx,ny), gky(nx,ny), rcell(3,3)
   logical,  intent(in)    :: hyb, terms
   real(dp), intent(out)   :: Weff(nx,ny,5)
   integer,  intent(inout) :: nsing

   real(dp), allocatable :: T(:,:)

   allocate(T(nx,ny))
   call semicl_coarea_div(W(:,:,3), Esheet, gkx, gky, nx, ny, rcell, T, nsing)
   if (hyb .or. terms) Weff(:,:,4) = T       ! D(w_2), free: already built here
   ! D(D(w_2)) alone = the Phi_2'' term of the co-area route
   if (terms) call semicl_coarea_div(Weff(:,:,4), Esheet, gkx, gky, nx, ny, rcell, &
                                     Weff(:,:,5), nsing)
   T = W(:,:,2) - T
   call semicl_coarea_div(T, Esheet, gkx, gky, nx, ny, rcell, Weff(:,:,1), nsing)
   Weff(:,:,1) = W(:,:,1) - Weff(:,:,1)
   call semicl_coarea_div(W(:,:,4), Esheet, gkx, gky, nx, ny, rcell, T, nsing)
   Weff(:,:,2) = -T
   ! kept last so that with hyb off nsing counts exactly what it used to
   if (hyb .or. terms) call semicl_coarea_div(W(:,:,2), Esheet, gkx, gky, nx, ny, rcell, &
                                              Weff(:,:,3), nsing)
   deallocate(T)

end subroutine SemiclSuscFields

!=======================================================================
! Per-orbit second-order Onsager term (L1c/L2).
!
! The O(B^2) term of the Roth-Gao-Niu quantisation is (B^2/2) chi0'(E), and
! with chi0' = (e^2/12hbar^2)(1/4pi^2) Psi (weights from diag.F90's
! OrbSuscAtK, where the formula is),
!       Psi = Phi_0 - dPhi_1/dE + d2Phi_2/dE2 ,   Phi_s = closed_int w_s dl/|grad E| .
! Two independent evaluations of the same Psi:
!   psi   : ONE line integral of the co-area field w_0 - D(w_1 - D(w_2));
!   psied : Phi_s on this orbit and on the matched orbits at lev +- h, and
!           central differences in energy.  The matched orbit is the closed,
!           same-sign contour nearest in centroid whose area moved by h*dA/dE
!           to within half of that; otherwise psied = -9.99e9.
! They share nothing but the weights and the contour tracer, so their
! agreement is a diagnostic of the derivative, not a duplication.
! psiLP = the single-band part alone (co-area), the control that must NOT
! vanish where the interband terms cancel (massive Dirac).
! The line integrals are unsigned: chi0' is a Fermi-sea derivative and does
! not care which side is occupied; the hole sign enters in semicl_chi_perT.
!=======================================================================
subroutine SemiclOrbitSusc(Esheet, W, Weff, gkx, gky, nx, ny, xs, ys, np, area, rcell, &
                           lev, h, hyb, terms, psi, psied, psihyb, psiLP, t0, t1, t2, nweak)

   integer,  intent(in)    :: nx, ny, np
   real(dp), intent(in)    :: Esheet(nx,ny), W(nx,ny,4), Weff(nx,ny,5)
   real(dp), intent(in)    :: gkx(nx,ny), gky(nx,ny), xs(np), ys(np), area, rcell(3,3)
   real(dp), intent(in)    :: lev, h
   logical,  intent(in)    :: hyb, terms
   real(dp), intent(out)   :: psi, psied, psihyb, psiLP
   ! the three terms of the co-area Psi, separately: psi == t0 - t1 + t2
   real(dp), intent(out)   :: t0, t1, t2
   integer,  intent(inout) :: nweak

   real(dp), allocatable :: bx(:), by(:)
   real(dp) :: J1(1), J2(2), J3(3), L0, Lm, phi1(-1:1), phi2(-1:1), c0(2)
   real(dp) :: phi1p, dw2(-1:1)
   integer  :: sd, npb, ndum, ndum0
   logical  :: ok

   ndum0 = 0

   call orbit_line_int(Weff, 2, gkx, gky, nx, ny, xs, ys, np, rcell, J2, L0, nweak)
   psi = J2(1) ; psiLP = J2(2)

   t0 = -9.99e9_dp ; t1 = -9.99e9_dp ; t2 = -9.99e9_dp
   if (terms) then
      ! Phi_0 needs no divergence; Phi_1' = int D(w_1); Phi_2'' = int D(D(w_2))
      call orbit_line_int(W(:,:,1:1), 1, gkx, gky, nx, ny, xs, ys, np, rcell, J1, Lm, ndum0)
      t0 = J1(1)
      call orbit_line_int(Weff(:,:,3:3), 1, gkx, gky, nx, ny, xs, ys, np, rcell, J1, Lm, ndum0)
      t1 = J1(1)
      call orbit_line_int(Weff(:,:,5:5), 1, gkx, gky, nx, ny, xs, ys, np, rcell, J1, Lm, ndum0)
      t2 = J1(1)
   end if

   psied = -9.99e9_dp ; psihyb = -9.99e9_dp
   if (h <= 0.0_dp) return
   ndum = 0
   call orbit_line_int(W, 3, gkx, gky, nx, ny, xs, ys, np, rcell, J3, L0, ndum)
   phi1(0) = J3(2) ; phi2(0) = J3(3)
   ! Phi_1' on THIS orbit, one co-area divergence and no energy derivative
   phi1p = 0.0_dp
   if (hyb) then
      call orbit_line_int(Weff(:,:,3:3), 1, gkx, gky, nx, ny, xs, ys, np, rcell, J1, Lm, ndum)
      phi1p = J1(1)
   end if
   c0(1) = sum(xs(1:np))/real(np,dp) ; c0(2) = sum(ys(1:np))/real(np,dp)
   allocate(bx(nx*ny), by(nx*ny))
   do sd = -1, 1, 2
      ! the signed area always grows with E at the rate L0 (electron and hole)
      call susc_match_orbit(Esheet, nx, ny, lev + real(sd,dp)*h, area + real(sd,dp)*h*L0, &
                            0.5_dp*h*L0, c0, rcell, bx, by, npb, ok)
      if (.not. ok) then
         deallocate(bx, by) ; return
      end if
      call orbit_line_int(W(:,:,2:3), 2, gkx, gky, nx, ny, bx, by, npb, rcell, J2, Lm, ndum)
      phi1(sd) = J2(1) ; phi2(sd) = J2(2)
      ! the matched orbits are already traced: the hybrid Phi_2'' rides along
      if (hyb) then
         call orbit_line_int(Weff(:,:,4:4), 1, gkx, gky, nx, ny, bx, by, npb, rcell, J1, Lm, ndum)
         dw2(sd) = J1(1)
      end if
   end do
   deallocate(bx, by)
   psied = J3(1) - (phi1(1) - phi1(-1))/(2.0_dp*h) + (phi2(1) - 2.0_dp*phi2(0) + phi2(-1))/(h*h)
   ! Phi_0 - Phi_1' + d/dE[Phi_2'], one k-divergence per term instead of two
   if (hyb) psihyb = J3(1) - phi1p + (dw2(1) - dw2(-1))/(2.0_dp*h)

end subroutine SemiclOrbitSusc

!-----------------------------------------------------------------------
! The closed, non-winding contour at lev with the sign of aexp whose centroid
! is nearest c0 (fractional, modulo the lattice); ok if its signed area is
! within tol of aexp.
!-----------------------------------------------------------------------
subroutine susc_match_orbit(Esheet, nx, ny, lev, aexp, tol, c0, rcell, bx, by, npb, ok)

   integer,  intent(in)  :: nx, ny
   real(dp), intent(in)  :: Esheet(nx,ny), lev, aexp, tol, c0(2), rcell(3,3)
   real(dp), intent(out) :: bx(nx*ny), by(nx*ny)
   integer,  intent(out) :: npb
   logical,  intent(out) :: ok

   real(dp), allocatable :: cx(:), cy(:), xs(:), ys(:)
   integer,  allocatable :: cstart(:), cend(:)
   logical,  allocatable :: closed(:)
   integer  :: nc, c, np, wx, wy
   logical  :: wound
   real(dp) :: ar, abest, dfx, dfy, dkx, dky, d2, best

   ok = .false. ; npb = 0 ; best = huge(1.0_dp) ; abest = 0.0_dp
   allocate(cx(nx*ny), cy(nx*ny), xs(nx*ny), ys(nx*ny))
   allocate(cstart(MAXCONT), cend(MAXCONT), closed(MAXCONT))
   call SemiclMarchingSquares(Esheet, nx, ny, lev, cx, cy, cstart, cend, nc, closed)
   do c = 1, nc
      np = cend(c) - cstart(c) + 1
      if (np < 3) cycle
      xs(1:np) = cx(cstart(c):cend(c))
      ys(1:np) = cy(cstart(c):cend(c))
      call SemiclUnwrap(xs, ys, np, wound, wx, wy)
      if (wound .or. .not. closed(c)) cycle
      call SemiclAreaCart(xs, ys, np, rcell, ar)
      if (ar*aexp <= 0.0_dp) cycle
      dfx = sum(xs(1:np))/real(np,dp) - c0(1) ; dfx = dfx - anint(dfx)
      dfy = sum(ys(1:np))/real(np,dp) - c0(2) ; dfy = dfy - anint(dfy)
      dkx = rcell(1,1)*dfx + rcell(1,2)*dfy
      dky = rcell(2,1)*dfx + rcell(2,2)*dfy
      d2  = dkx*dkx + dky*dky
      if (d2 < best) then
         best = d2 ; abest = ar ; npb = np
         bx(1:np) = xs(1:np) ; by(1:np) = ys(1:np)
      end if
   end do
   ok = npb > 0 .and. abs(abest - aexp) <= tol
   deallocate(cx, cy, xs, ys, cstart, cend, closed)

end subroutine susc_match_orbit

!-----------------------------------------------------------------------
! Second-order offset per tesla (L1c/L2).  The Roth-Gao-Niu condition in
! semicl's form is
!     SEMICL_C*|A| = B * (N + gamma + delta_m + delta_chi(B)),
!     delta_chi(B) = -sgn * B * Psi / (48 pi hbar/e)  =  dchi * B ,
! and hbar/e = 2 pi SEMICL_C [T Ang^2], so dchi = -sgn Psi/(96 pi^2 SEMICL_C)
! with Psi in Ang^2 -- no new constant.  sgn = +1 electron, -1 hole: for a
! hole orbit every correction changes sign (Fuchs et al. Eq. 8), while Psi is
! built from unsigned line integrals.
!-----------------------------------------------------------------------
pure real(dp) function semicl_chi_perT(psi, sgn)
   use constants, only : pi
   real(dp), intent(in) :: psi, sgn
   semicl_chi_perT = -sgn*psi/(96.0_dp*pi*pi*SEMICL_C)
end function semicl_chi_perT

!-----------------------------------------------------------------------
! B from  dchi*B^2 + nu*B - S = 0  (S = SEMICL_C*|A|, nu = N+gamma+delta_m):
! the root that tends to S/nu as dchi -> 0, in the cancellation-free form.
! -1 when there is no positive real root -- the second-order term is then not
! a small correction and the expansion has broken down for this level.
!-----------------------------------------------------------------------
pure real(dp) function semicl_onsager_B(S, nu, dchi)
   real(dp), intent(in) :: S, nu, dchi
   real(dp) :: disc, den
   semicl_onsager_B = -1.0_dp
   disc = nu*nu + 4.0_dp*dchi*S
   if (disc < 0.0_dp) return
   den = nu + sqrt(disc)
   if (den <= 0.0_dp) return
   semicl_onsager_B = 2.0_dp*S/den
end function semicl_onsager_B

!=======================================================================
! Driver.  Stages, each skippable when its output already exists:
!
!   1. eigenvalue grid  <- generate.bands from Diag.Calculate3DTAPWBands
!                          (existing, already-validated code path)
!   2. contours + areas -> <prefix>.FermiAreas          [checkpoint]
!   3. Onsager fan      -> <prefix>.LLfan, <prefix>.Wannier
!
! Stage 3 is arithmetic, so gamma / NumLL / B-range can be changed and
! re-run from the FermiAreas checkpoint without re-diagonalising.
!=======================================================================
subroutine SemiclassicalLL()

   use cell,  only : rcell
   use name,  only : prefix

   character(len=200) :: bandsFile, fareas, bfluxFile
   integer  :: nx, ny, nb, bmin, bmax, numE, numLL
   real(dp) :: emin, emax, gamma
   logical  :: ok, restart, haveAreas, domass, useBerry, okB

   real(dp), allocatable :: Egrid(:,:,:), sheetgap(:)
   real(dp), allocatable :: Bflux(:,:,:)
   real(dp), allocatable :: Mz(:,:,:)
   integer :: b, bfb1, bfb2, mfb1, mfb2
   logical :: useMom, okM
   character(len=200) :: momFile
   logical  :: doBD, okBD
   real(dp) :: bdMinGap, bdGamma
   character(len=200) :: bdFile
   logical  :: useChi, okC
   character(len=200) :: chiFile
   logical  :: doWilson
   character(len=200) :: linksFile
   real(dp) :: wilWin, wilIso, wilMinDet
   real(dp), allocatable :: Wc(:,:,:,:)
   integer  :: cfb1, cfb2

#ifdef DEBUG
   call MIO_Debug('SemiclassicalLL',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('semicl')
#endif /* TIMER */

   call MIO_InputParameter('Semicl.BandsFile', bandsFile, 'generate.bands')
   call MIO_InputParameter('Semicl.GridX',     nx,     48)
   call MIO_InputParameter('Semicl.GridY',     ny,     48)
   call MIO_InputParameter('Semicl.BandMin',   bmin,   1)
   call MIO_InputParameter('Semicl.BandMax',   bmax,   0)     ! 0 = all
   call MIO_InputParameter('Semicl.EnergyMin', emin,  -0.20_dp)
   call MIO_InputParameter('Semicl.EnergyMax', emax,   0.20_dp)
   call MIO_InputParameter('Semicl.NumE',      numE,   256)
   call MIO_InputParameter('Semicl.NumLL',     numLL,  10)
   call MIO_InputParameter('Semicl.Gamma',     gamma,  0.5_dp)
   call MIO_InputParameter('Semicl.Restart',   restart, .true.)
   ! pure post-processing of the areas: writes an extra file, changes none
   call MIO_InputParameter('Semicl.Masses',    domass,  .true.)
   ! per-orbit gamma from the Berry phase (E3).  Needs <prefix>.BerryFlux from
   ! a Diag.BerryFlux run on the SAME k-grid; default .false. keeps the old
   ! behaviour exactly.  Semicl.Gamma < 0 then means "use the computed value".
   call MIO_InputParameter('Semicl.BerryPhase', useBerry, .false.)
   call MIO_InputParameter('Semicl.BerryFluxFile', bfluxFile, trim(prefix)//'.BerryFlux')
   ! orbital-moment correction to the Onsager offset (L1).  Needs
   ! <prefix>.OrbMoment from a Diag.OrbMoment run on the SAME k-grid; default
   ! .false. keeps the old behaviour exactly.  The Berry phase alone is NOT the
   ! complete first-order semiclassical correction -- the orbital moment enters
   ! at the same order -- so this is a correctness fix, not a refinement.
   call MIO_InputParameter('Semicl.OrbitalMoment', useMom, .false.)
   call MIO_InputParameter('Semicl.OrbMomentFile', momFile, trim(prefix)//'.OrbMoment')
   ! second-order Onsager term (L1c/L2): (B^2/2) chi0'(E) per valley, from the
   ! weights in <prefix>.OrbSusc of a Diag.OrbSusc run on the SAME k-grid.
   ! Default .false. keeps the old behaviour exactly.
   call MIO_InputParameter('Semicl.Susceptibility', useChi, .false.)
   call MIO_InputParameter('Semicl.OrbSuscFile', chiFile, trim(prefix)//'.OrbSusc')
   ! Third, independent evaluation of the same Psi, keeping ONE co-area
   ! divergence instead of two (the nested one cannot converge near a critical
   ! point -- see SemiclSuscFields).  A diagnostic column only: dchi, and so
   ! every level written, still comes from the co-area psi.
   call MIO_InputParameter('Semicl.PsiHybrid', psiHybOn, .false.)
   ! Discretisation of the co-area divergence.  The differenced form carries a
   ! GRID-INDEPENDENT integrated error (measured flat at 15.24 on an analytic
   ! cosine band from 48^2 to 384^2, against 1.1e-1 -> 2.9e-3 for this one),
   ! which is why refining the grid never converged Psi.  Default .false.
   call MIO_InputParameter('Semicl.CoareaAnalytic', coareaAn, .false.)
   call MIO_InputParameter('Semicl.PsiTerms', psiTermsOn, .false.)
   ! L5: path-ordered (non-Abelian) Wilson loop per orbit, from the multi-band
   ! link matrices a Diag.BerryLinks run writes on the SAME k-grid.  Default
   ! .false.: nothing new is read, computed or written.
   call MIO_InputParameter('Semicl.Wilson', doWilson, .false.)
   call MIO_InputParameter('Semicl.BerryLinksFile', linksFile, trim(prefix)//'.BerryLinks')
   ! half-width [eV] of the pointwise multiplet window around the orbit energy
   call MIO_InputParameter('Semicl.WilsonWindow', wilWin, 0.05_dp)
   ! the window must stay this far [eV] from the nearest state outside it, at
   ! EVERY node of the path, or the orbit is rejected
   call MIO_InputParameter('Semicl.WilsonIsolation', wilIso, 0.01_dp)
   ! |det W| gate: a unitary connection gives 1, and a shortfall means the
   ! selected subspace did not close around the loop
   call MIO_InputParameter('Semicl.WilsonMinDet', wilMinDet, 0.90_dp)
   ! strong magnetic breakdown (E7): composite orbits between full spectral
   ! gaps, written to <prefix>.BreakdownAreas and quantised into .LLfanBD /
   ! .WannierBD.  Default .false.: nothing new is computed or written.
   call MIO_InputParameter('Semicl.Breakdown', doBD, .false.)
   ! a full gap narrower than this [eV] is not taken as an anchor (grid noise)
   call MIO_InputParameter('Semicl.BreakdownMinGap', bdMinGap, 1.0e-4_dp)
   ! gamma for the breakdown fan; default = Semicl.Gamma, so < 0 means "use the
   ! composite Berry phase".  Separate because that composite phase is NOT yet
   ! the diabatic Berry phase where several sheets carry the orbit (see the
   ! SemiclBreakdownAreas header): measured on G/hBN, a fixed Dirac value 0.0
   ! scores better against the exact spectrum there.
   call MIO_InputParameter('Semicl.BreakdownGamma', bdGamma, gamma)
   bdFile = trim(prefix)//'.BreakdownAreas'

   fareas = trim(prefix)//'.FermiAreas'

   call MIO_Print('','semicl')
   call MIO_Print('=== semiclassical Landau levels ==================','semicl')
   call MIO_Print('  B[T] = '//trim(num2str(SEMICL_C,4))// &
        ' * A[Ang^-2] / (N+gamma)','semicl')

   haveAreas = .false.
   if (restart) inquire(file=trim(fareas), exist=haveAreas)

   if (haveAreas) then
      call MIO_Print('  found '//trim(fareas)//' -> skipping to Onsager stage','semicl')
   else
      call SemiclReadBands3D(bandsFile, nx, ny, Egrid, nb, ok)
      if (.not. ok) then
         call MIO_Print('  stage 1 failed; nothing written','semicl')
         return
      end if
      if (bmax <= 0 .or. bmax > nb) bmax = nb
      if (bmin < 1) bmin = 1

      ! minimum gap to the next sheet: flags magnetic-breakdown regions,
      ! where a single-sheet orbit area stops being the right quantity
      allocate(sheetgap(nb))
      sheetgap = huge(1.0_dp)
      do b = 1, nb-1
         sheetgap(b) = minval(Egrid(:,:,b+1) - Egrid(:,:,b))
      end do
      if (nb >= 1) sheetgap(nb) = huge(1.0_dp)

      bfb1 = 0 ; bfb2 = -1
      if (useBerry) then
         call SemiclReadBerryFlux(bfluxFile, nx, ny, Bflux, bfb1, bfb2, okB)
         if (.not. okB) call MIO_Print('  continuing without per-orbit gamma','semicl')
      end if

      mfb1 = 0 ; mfb2 = -1
      if (useMom) then
         call SemiclReadOrbMoment(momFile, nx, ny, Mz, mfb1, mfb2, okM)
         if (.not. okM) call MIO_Print('  continuing without the orbital-moment '// &
              'correction','semicl')
      end if

      cfb1 = 0 ; cfb2 = -1
      if (useChi) then
         call SemiclReadOrbSusc(chiFile, nx, ny, Wc, cfb1, cfb2, okC)
         if (.not. okC) call MIO_Print('  continuing without the second-order term','semicl')
      end if

      ! The checkpoints are passed as ALLOCATABLE dummies rather than as
      ! optional arguments: with independent extras the optional form needs
      ! one call site per combination.
      call SemiclContoursAreas(Egrid, nx, ny, nb, bmin, bmax, &
                               emin, emax, numE, rcell, sheetgap, fareas, &
                               Bflux, bfb1, Mz, mfb1, Wc, cfb1)
      if (doBD) call SemiclBreakdownAreas(Egrid, nx, ny, nb, emin, emax, numE, rcell, &
                                          bdFile, Bflux, bfb1, bdMinGap)
      ! L5: path-ordered Wilson loop.  Independent of the Abelian gamma above --
      ! it re-traces the contours and consumes <prefix>.BerryLinks instead.
      if (doWilson) call SemiclWilsonLoops(Egrid, nx, ny, nb, bmin, bmax, emin, emax, &
                                           numE, rcell, trim(prefix)//'.Wilson', &
                                           linksFile, wilWin, wilIso, wilMinDet)
      if (allocated(Bflux)) deallocate(Bflux)
      if (allocated(Mz))    deallocate(Mz)
      if (allocated(Wc))    deallocate(Wc)
      deallocate(Egrid, sheetgap)
   end if

   call SemiclOnsager(fareas, numLL, gamma)
   if (domass) call SemiclMasses(fareas)
   if (doBD) then
      inquire(file=trim(bdFile), exist=okBD)
      if (okBD) then
         ! the composite orbit carries no delta_m yet, so the moment is off here
         call SemiclOnsager(bdFile, numLL, bdGamma, 'BD', .true.)
      else
         call MIO_Print('  Semicl.Breakdown: no '//trim(bdFile)//' (restarted from '// &
              'an older .FermiAreas?) -- rerun with Semicl.Restart .false.','semicl')
      end if
   end if

   call MIO_Print('=================================================','semicl')

#ifdef TIMER
   call MIO_TimerStop('semicl')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('SemiclassicalLL',1)
#endif /* DEBUG */

end subroutine SemiclassicalLL

!=======================================================================
! Stage 2+3: contours and their areas, written to <prefix>.FermiAreas.
!=======================================================================
subroutine SemiclContoursAreas(Egrid, nx, ny, nb, bmin, bmax, &
                               emin, emax, numE, rcell, sheetgap, fname, &
                               Bflux, bfb1, Mz, mfb1, Wc, cfb1)

   use constants, only : pi

   integer,  intent(in) :: nx, ny, nb, bmin, bmax, numE
   real(dp), intent(in) :: Egrid(nx,ny,nb), emin, emax, rcell(3,3), sheetgap(nb)
   character(*), intent(in) :: fname
   ! per-plaquette Berry fluxes for bands bfb1 .. bfb1+size(,3)-1, and
   ! per-grid-point orbital moments for mfb1 ..  Unallocated = not available;
   ! each is independent of the other.
   real(dp), allocatable, intent(in) :: Bflux(:,:,:)
   integer,  intent(in) :: bfb1
   real(dp), allocatable, intent(in) :: Mz(:,:,:)
   integer,  intent(in) :: mfb1
   ! second-order weights (w_0, w_1, w_2, w_LP) per grid point for bands
   ! cfb1 .. (L1c/L2); unallocated = not available
   real(dp), allocatable, intent(in) :: Wc(:,:,:,:)
   integer,  intent(in) :: cfb1

   real(dp), allocatable :: cx(:), cy(:), xs(:), ys(:)
   integer,  allocatable :: cstart(:), cend(:)
   logical,  allocatable :: closed(:)
   integer :: np
   integer  :: b, ie, c, nc, wx, wy, u, ntot, nopen
   real(dp) :: lev, area, blo, bhi, perim, dE
   real(dp) :: A_BZ, occ, A_count, A_sum_signed, A_sum_abs
   real(dp), allocatable :: fill(:)
   real(dp) :: dk, kF, relerr
   integer  :: nopen_here, nbad_here, nbad
   logical  :: valid
   logical  :: wound
   logical  :: haveB, bmatch
   integer  :: ibf, ngam
   real(dp) :: phaseB, gam
   logical  :: haveM, haveExtra
   integer  :: imf, ndel, nweak
   real(dp) :: delm, dAdE_line
   real(dp), allocatable :: gkx(:,:), gky(:,:)
   logical  :: haveC, bandC
   integer  :: icf, nchi, nnoed, nsing, nnohy
   real(dp) :: psi, psied, psiLP, dchi, psihyb, pt0, pt1, pt2
   real(dp), allocatable :: Weff(:,:,:)

   allocate(cx(nx*ny), cy(nx*ny), xs(nx*ny), ys(nx*ny))
   allocate(cstart(MAXCONT), cend(MAXCONT), closed(MAXCONT))

   u = 772
   open(u, file=trim(fname), status='replace', action='write')
   write(u,'(a)') '# semiclassical Fermi contour areas'
   write(u,'(a)') '# B[T] = 10475.7686 * area / (N+gamma)   [area in Ang^-2]'
   write(u,'(a,i0,a,i0)') '# grid ', nx, ' x ', ny
   write(u,'(a,i0,a,i0)') '# bands ', bmin, ' .. ', bmax
   write(u,'(a,es14.6,a,es14.6,a,i0)') '# E from ', emin, ' to ', emax, ' in ', numE
   haveB = allocated(Bflux)
   haveM = allocated(Mz)
   ! The row format widens only when a new quantity is actually available, so a
   ! Berry-only run still writes exactly the 14-column file it used to and its
   ! checkpoints stay byte-for-byte reproducible.  Readers take the widest
   ! record that parses.
   haveC = allocated(Wc)
   haveExtra = haveM .or. haveC
   if (haveC) then
      write(u,'(a)') '# iband  E[eV]        area[Ang^-2]   signed_area    perim[Ang^-1]'// &
                     '  npts  closed  wx  wy  sheetgap[eV]   n/n0          est_rel_err'// &
                     '    berryphase[rad]  gamma            delta_m          dAdE_line'// &
                     '        psi[Ang^2]       psi_ediff[Ang^2] psi_LP[Ang^2]    dchi[1/T]'// &
                     trim(merge('        psi_hybrid[Ang^2]', '                         ', psiHybOn))
      write(u,'(a)') "#   psi = Phi_0 - dPhi_1/dE + d2Phi_2/dE2, Phi_s = closed_int w_s dl/|grad E|:"
      write(u,'(a)') "#   chi0' = (e^2/12hbar^2)(1/4pi^2) psi, the SECOND-ORDER term (L1c/L2), from the"
      write(u,'(a)') '#   co-area field w_0 - D(w_1 - D(w_2)); psi_ediff = the same by central energy'
      write(u,'(a)') '#   differences of Phi_s on the matched orbits at E +- dE (a cross-check);'
      write(u,'(a)') '#   psi_LP = the single-band (Landau-Peierls) part alone (a control).'
      write(u,'(a)') '#   dchi = -sgn*psi/(96 pi^2 SEMICL_C):  SEMICL_C*area = B*(N+gamma+delta_m) + dchi*B^2'
      if (psiHybOn) then
         write(u,'(a)') '#   psi_hybrid = Phi_0 - closed_int D(w_1) + d/dE[closed_int D(w_2)]: the same Psi with'
         write(u,'(a)') '#   ONE co-area divergence per term instead of the nested D(D(.)), whose peak grows as'
         write(u,'(a)') '#   1/h^4 near a critical point.  A DIAGNOSTIC column: dchi still comes from psi.'
      end if
      if (psiTermsOn) then
         write(u,'(a)') '#   t0, t1, t2 = Phi_0, Phi_1'' = int D(w_1), Phi_2'''' = int D(D(w_2)), the three terms'
         write(u,'(a)') '#   of the co-area psi separately (D is linear, so the split is exact and'
         write(u,'(a)') '#   psi = t0 - t1 + t2).  Their SIZE against psi is the cancellation ratio: if they'
         write(u,'(a)') '#   are orders of magnitude larger than psi, no discretisation of D can converge psi.'
      end if
   end if
   if (haveExtra) then
      if (.not. haveC) &
      write(u,'(a)') '# iband  E[eV]        area[Ang^-2]   signed_area    perim[Ang^-1]'// &
                     '  npts  closed  wx  wy  sheetgap[eV]   n/n0          est_rel_err'// &
                     '    berryphase[rad]  gamma            delta_m          dAdE_line'
      write(u,'(a)') '#   gamma = 1/2 - Phi_B/2pi, per ORBIT, sign kept (it is the valley-odd'
      write(u,'(a)') '#   quantity).  Phi_B = sum of the plaquette fluxes strictly inside the'
      write(u,'(a)') '#   orbit; -9.99e9 in a column means that quantity is NOT available for'
      write(u,'(a)') '#   this contour (no interior component matched, or the checkpoint was'
      write(u,'(a)') '#   not supplied).'
      write(u,'(a)') '#   delta_m = (sgn/2pi) * closed_int M dl/|grad E|, the ORBITAL MOMENT'
      write(u,'(a)') '#   correction: B = SEMICL_C*area/(N + gamma + delta_m).  It enters at'
      write(u,'(a)') '#   the same semiclassical order as the Berry phase, and for a massive'
      write(u,'(a)') '#   Dirac orbit it cancels gamma exactly (gamma+delta_m -> integer).'
      write(u,'(a)') '#   dAdE_line = the same line integral without M, i.e. dA/dE from THIS'
      write(u,'(a)') '#   contour alone -- independent of the central differences in .Masses.'
   else if (haveB) then
      write(u,'(a)') '# iband  E[eV]        area[Ang^-2]   signed_area    perim[Ang^-1]'// &
                     '  npts  closed  wx  wy  sheetgap[eV]   n/n0          est_rel_err'// &
                     '    berryphase[rad]  gamma'
      write(u,'(a)') '#   gamma = 1/2 - Phi_B/2pi, per ORBIT, sign kept (it is the valley-odd'
      write(u,'(a)') '#   quantity).  Phi_B = sum of the plaquette fluxes strictly inside the'
      write(u,'(a)') '#   orbit; -9.99e9 in both columns means no interior component could be'
      write(u,'(a)') '#   matched to this contour and gamma is NOT available for it.'
   else
      write(u,'(a)') '# iband  E[eV]        area[Ang^-2]   signed_area    perim[Ang^-1]'// &
                     '  npts  closed  wx  wy  sheetgap[eV]   n/n0          est_rel_err'
   end if
   write(u,'(a)') '#SUM iband  E[eV]        sum|area|      sum_signed     A_count[Ang^-2]'// &
                  '  occ_frac   ncont  nopen'
   write(u,'(a)') '#   est_rel_err ~ (dk/k_F)^2 with dk = |b1|/nx and k_F = sqrt(A/pi):'
   write(u,'(a)') '#   the polygon-discretisation error of THIS contour.  Small pockets near a'
   write(u,'(a)') '#   fan origin are traced by only a handful of cells and can carry tens of'
   write(u,'(a)') '#   percent; filter on this before trusting an individual orbit area.'
   write(u,'(a)') '#   A_count = A_BZ * (fraction of grid points with E_b<E), an independent'
   write(u,'(a)') '#   estimate of the enclosed area; it must agree with the contour area to'
   write(u,'(a)') '#   O(1/N_grid) for electron pockets (or with A_BZ - area for hole pockets).'

   dE = (emax-emin)/real(max(numE-1,1),dp)
   dk = sqrt(rcell(1,1)**2 + rcell(2,1)**2)/real(nx,dp)   ! |b1|/nx
   A_BZ = abs(rcell(1,1)*rcell(2,2) - rcell(1,2)*rcell(2,1))
   write(u,'(a,es14.6)') '# A_BZ[Ang^-2] = ', A_BZ
   write(u,'(a,i0)') '# nbands_total = ', nb

   ! Total filling n/n0 in states per moire cell, summed over EVERY band of the
   ! grid (not just bmin..bmax), measured from charge neutrality at nb/2.
   ! Per valley, per spin -- multiply by the spin/valley degeneracy downstream.
   ! phi/phi0 = B / (Phi0/A_moire) = B * A_BZ / (2*pi)^2 ... written by Onsager.
   allocate(fill(numE))
   do ie = 1, numE
      lev = emin + real(ie-1,dp)*dE
      fill(ie) = real(count(Egrid < lev),dp)/real(nx*ny,dp) - 0.5_dp*real(nb,dp)
   end do

   ntot = 0 ; nopen = 0 ; nbad = 0 ; ngam = 0 ; ndel = 0 ; nweak = 0
   nchi = 0 ; nnoed = 0 ; nsing = 0 ; nnohy = 0
   if (haveM .or. haveC) allocate(gkx(nx,ny), gky(nx,ny))
   if (haveC) allocate(Weff(nx,ny,5))

   do b = bmin, bmax
      blo = minval(Egrid(:,:,b)) ; bhi = maxval(Egrid(:,:,b))
      ! |grad_k E| weights the moment line integral; it depends only on the
      ! sheet, so it is built once per band rather than once per contour.
      if (haveM .or. haveC) call SemiclGradient(Egrid(:,:,b), nx, ny, rcell, gkx, gky)
      ! the co-area fields of the second-order term, likewise once per sheet
      icf = b - cfb1 + 1
      bandC = .false.
      if (haveC) bandC = icf >= 1 .and. icf <= size(Wc,4)
      if (bandC) call SemiclSuscFields(Wc(:,:,:,icf), Egrid(:,:,b), gkx, gky, nx, ny, rcell, &
                                       psiHybOn, psiTermsOn, Weff, nsing)
      do ie = 1, numE
         lev = emin + real(ie-1,dp)*dE
         if (lev <= blo .or. lev >= bhi) cycle     ! no contour on this sheet

         call SemiclMarchingSquares(Egrid(:,:,b), nx, ny, lev, cx, cy, cstart, cend, nc, closed)

         ! independent area estimate by counting occupied grid points
         occ     = real(count(Egrid(:,:,b) < lev),dp)/real(nx*ny,dp)
         A_count = occ*A_BZ
         A_sum_signed = 0.0_dp ; A_sum_abs = 0.0_dp
         nopen_here = 0 ; nbad_here = 0

         do c = 1, nc
            np = cend(c) - cstart(c) + 1
            if (np < 3) cycle
            xs(1:np) = cx(cstart(c):cend(c))
            ys(1:np) = cy(cstart(c):cend(c))
            call SemiclUnwrap(xs, ys, np, wound, wx, wy)
            call SemiclAreaCart(xs, ys, np, rcell, area)
            call contour_perimeter(xs, ys, np, rcell, perim)
            ntot = ntot + 1
            ! Three outcomes, only the first has a meaningful cyclotron area:
            !   closed .and. .not.wound : genuine closed orbit          -> area
            !   wound                   : open orbit across the BZ      -> no area
            !   .not.closed .and. .not.wound : TRUNCATED trace (the march hit an
            !       unresolvable saddle, code 5/10, or escaped). The endpoints are
            !       near each other so the winding test does NOT catch it, and a
            !       partial polygon would otherwise be handed a plausible but
            !       meaningless area.  Excluded and counted separately.
            valid = closed(c) .and. .not. wound
            if (wound) then
               nopen = nopen + 1 ; nopen_here = nopen_here + 1
            else if (.not. closed(c)) then
               nbad = nbad + 1 ; nbad_here = nbad_here + 1
            else
               A_sum_signed = A_sum_signed + area
               A_sum_abs    = A_sum_abs    + abs(area)
            end if
            ! polygon-discretisation error estimate for THIS orbit
            if (valid .and. abs(area) > 0.0_dp) then
               kF = sqrt(abs(area)/pi)
               relerr = (dk/kF)**2
            else
               relerr = -1.0_dp
            end if
            ! per-orbit Berry phase and Onsager offset (E3), and the
            ! orbital-moment correction that enters at the same order (L1)
            if (haveExtra) then
               phaseB = -9.99e9_dp ; gam = -9.99e9_dp
               delm   = -9.99e9_dp ; dAdE_line = -9.99e9_dp
               ibf = b - bfb1 + 1
               if (haveB .and. valid .and. ibf >= 1 .and. ibf <= size(Bflux,3)) then
                  call SemiclOrbitBerry(Egrid(:,:,b), Bflux(:,:,ibf), nx, ny, &
                                        lev, area, A_BZ, phaseB, bmatch)
                  if (bmatch) then
                     gam = 0.5_dp - phaseB/(2.0_dp*pi)
                     ngam = ngam + 1
                  else
                     phaseB = -9.99e9_dp
                  end if
               end if
               imf = b - mfb1 + 1
               if (haveM .and. valid .and. imf >= 1 .and. imf <= size(Mz,3)) then
                  call SemiclOrbitMoment(Mz(:,:,imf), gkx, gky, nx, ny, &
                                         xs, ys, np, rcell, sign(1.0_dp, area), &
                                         delm, dAdE_line, nweak)
                  ndel = ndel + 1
               end if
               if (haveC) then
                  psi = -9.99e9_dp ; psied = -9.99e9_dp ; psiLP = -9.99e9_dp ; dchi = -9.99e9_dp
                  psihyb = -9.99e9_dp
                  if (bandC .and. valid) then
                     ! h = the energy step of this file, so psi_ediff uses the
                     ! neighbouring energies' contours
                     call SemiclOrbitSusc(Egrid(:,:,b), Wc(:,:,:,icf), Weff, gkx, gky, nx, ny, &
                                          xs, ys, np, area, rcell, lev, dE, psiHybOn, &
                                          psiTermsOn, psi, psied, psihyb, psiLP, &
                                          pt0, pt1, pt2, nweak)
                     ! dchi stays on the co-area psi: the hybrid is a diagnostic
                     ! column, so enabling it changes no level that is written
                     dchi = semicl_chi_perT(psi, sign(1.0_dp, area))
                     nchi = nchi + 1
                     if (psied < -1.0e8_dp) nnoed = nnoed + 1
                     if (psiHybOn .and. psihyb < -1.0e8_dp) nnohy = nnohy + 1
                  end if
                  if (psiTermsOn) then
                  write(u,'(i6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,i6,2x,l1,2x,i3,2x,i3,'// &
                           '2x,es13.6,2x,es13.6,2x,es13.6,12(2x,es15.7))') &
                       b, lev, merge(abs(area), 0.0_dp, valid), area, perim, &
                       np, closed(c), wx, wy, sheetgap(b), fill(ie), relerr, phaseB, gam, &
                       delm, dAdE_line, psi, psied, psiLP, dchi, psihyb, pt0, pt1, pt2
                  else if (psiHybOn) then
                  write(u,'(i6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,i6,2x,l1,2x,i3,2x,i3,'// &
                           '2x,es13.6,2x,es13.6,2x,es13.6,9(2x,es15.7))') &
                       b, lev, merge(abs(area), 0.0_dp, valid), area, perim, &
                       np, closed(c), wx, wy, sheetgap(b), fill(ie), relerr, phaseB, gam, &
                       delm, dAdE_line, psi, psied, psiLP, dchi, psihyb
                  else
                  write(u,'(i6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,i6,2x,l1,2x,i3,2x,i3,'// &
                           '2x,es13.6,2x,es13.6,2x,es13.6,8(2x,es15.7))') &
                       b, lev, merge(abs(area), 0.0_dp, valid), area, perim, &
                       np, closed(c), wx, wy, sheetgap(b), fill(ie), relerr, phaseB, gam, &
                       delm, dAdE_line, psi, psied, psiLP, dchi
                  end if
               else
               write(u,'(i6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,i6,2x,l1,2x,i3,2x,i3,'// &
                        '2x,es13.6,2x,es13.6,2x,es13.6,2x,es15.7,2x,es15.7,2x,es15.7,2x,es15.7)') &
                    b, lev, merge(abs(area), 0.0_dp, valid), area, perim, &
                    np, closed(c), wx, wy, sheetgap(b), fill(ie), relerr, phaseB, gam, &
                    delm, dAdE_line
               end if
            else if (haveB) then
               phaseB = -9.99e9_dp ; gam = -9.99e9_dp
               ibf = b - bfb1 + 1
               if (valid .and. ibf >= 1 .and. ibf <= size(Bflux,3)) then
                  call SemiclOrbitBerry(Egrid(:,:,b), Bflux(:,:,ibf), nx, ny, &
                                        lev, area, A_BZ, phaseB, bmatch)
                  if (bmatch) then
                     gam = 0.5_dp - phaseB/(2.0_dp*pi)
                     ngam = ngam + 1
                  else
                     phaseB = -9.99e9_dp
                  end if
               end if
               write(u,'(i6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,i6,2x,l1,2x,i3,2x,i3,'// &
                        '2x,es13.6,2x,es13.6,2x,es13.6,2x,es15.7,2x,es15.7)') &
                    b, lev, merge(abs(area), 0.0_dp, valid), area, perim, &
                    np, closed(c), wx, wy, sheetgap(b), fill(ie), relerr, phaseB, gam
            else
               write(u,'(i6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,i6,2x,l1,2x,i3,2x,i3,2x,es13.6,2x,es13.6,2x,es13.6)') &
                    b, lev, merge(abs(area), 0.0_dp, valid), area, perim, &
                    np, closed(c), wx, wy, sheetgap(b), fill(ie), relerr
            end if
         end do

         ! per-(band,energy) consistency line
         write(u,'(a,i6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,f9.6,2x,i5,2x,i5)') &
              '#SUM ', b, lev, A_sum_abs, A_sum_signed, A_count, occ, nc, nopen_here+nbad_here
      end do
   end do

   close(u)
   call MIO_Print('  wrote '//trim(fname)//': '//trim(num2str(ntot))// &
        ' contours, '//trim(num2str(nopen))//' open/winding, '// &
        trim(num2str(nbad))//' truncated (saddle/escape)','semicl')
   if (nbad > 0) call MIO_Print('  NOTE: truncated traces are excluded from the '// &
        'fan; if there are many, refine the k-grid','semicl')
   if (haveB) call MIO_Print('  per-orbit gamma from the Berry phase: '// &
        trim(num2str(ngam))//' of '//trim(num2str(ntot))//' contours','semicl')
   if (haveM) then
      call MIO_Print('  per-orbit orbital-moment offset delta_m: '// &
           trim(num2str(ndel))//' of '//trim(num2str(ntot))//' contours','semicl')
      if (nweak > 0) call MIO_Print('  '//trim(num2str(nweak))//' contour segments '// &
           'sat at |grad E| < 1e-12 (a band extremum) and were dropped from the '// &
           'line integrals','semicl')
   end if
   if (haveC) then
      call MIO_Print('  per-orbit second-order term psi / dchi: '// &
           trim(num2str(nchi))//' of '//trim(num2str(ntot))//' contours; the energy-'// &
           'difference cross-check found no matched orbit at E +- dE for '// &
           trim(num2str(nnoed)),'semicl')
      if (nsing > 0) call MIO_Print('  '//trim(num2str(nsing))//' grid points with '// &
           '|grad E| <= 1e-12 got zero co-area flux (critical points)','semicl')
      if (psiHybOn) call MIO_Print('  hybrid route psi_hybrid (one co-area divergence '// &
           'per term): no matched orbit at E +- dE for '//trim(num2str(nnohy))// &
           ' contours','semicl')
      if (.not. haveM .and. nweak > 0) call MIO_Print('  '//trim(num2str(nweak))// &
           ' contour segments sat at |grad E| < 1e-12 and were dropped','semicl')
   end if

   if (allocated(Weff)) deallocate(Weff)
   if (allocated(gkx)) deallocate(gkx, gky)
   deallocate(cx, cy, cstart, cend, closed, xs, ys, fill)

end subroutine SemiclContoursAreas

!=======================================================================
! Strong magnetic breakdown (E7): composite orbits between spectral gaps.
!
! Breakdown couples the orbits of sheets that overlap in ENERGY: at a junction
! the carrier can cross from one sheet's contour to the other's, with the
! Landau-Zener probability P = exp(-B0/B).  When every such junction is
! transparent (P -> 1) the carrier no longer feels the avoided crossings, and
! the orbit it follows is the constant-energy contour of the DIABATIC bands --
! for a weak moire potential, the unfolded Dirac circle.  Its area needs no
! tracing: at every k the diabatic sheets hold exactly the same energies as
! the energy-sorted ones, so by Luttinger
!
!     A_e(E) = A_BZ * #{(k,b) : E_gap < E_b(k) < E} / N_k       (electron)
!     A_h(E) = A_BZ * #{(k,b) : E < E_b(k) < E_gap'} / N_k      (hole)
!
! counted from the nearest FULL spectral gap below (above) E.  A full gap is
! the one thing breakdown cannot cross -- there is no state to tunnel into --
! so the gaps are the anchors, and the Wannier fans of this limit emanate only
! from them.  Validated against the exact Hofstadter spectrum of G/hBN
! CellSize 55: at 3-10 T this fan matches the reference DOS peaks to 1.0-1.5
! meV over 0.15-0.35 eV, where the adiabatic sheet-57+ fan does not.
!
! Composite Berry phase: the sum, over every sheet that carries a contour at
! E, of its plaquette fluxes over {E_b < E} (electron) or {E_b > E} (hole),
! with the corner-fraction weighting of SemiclOrbitBerry.  Sheets lying wholly
! inside the composite contribute 2*pi*C_b and are left out, because gamma is
! only needed mod 1.  If any straddling sheet has no flux in the checkpoint
! the phase is reported as unavailable rather than summed partially.
!
! KNOWN LIMITATION of that phase (measured, G/hBN CellSize 55, K' and K): it
! is exactly valley-odd (gamma_K + gamma_K' = 1 to 3e-7) and reduces to the
! adiabatic per-orbit gamma where one sheet carries the orbit, but where
! several sheets do it drifts through non-integer values at the energies where
! a new sheet starts to carry a contour (0.17-0.20 and ~0.30 eV there), i.e.
! the adiabatic-mask sum is not yet the DIABATIC Berry phase.  Against the
! exact spectrum at 3-10 T a fixed gamma = 0 scores 0.99 / 1.67 / 2.07 meV in
! 150-218 / 218-300 / 300-345 meV vs 1.89 / 1.66 / 2.50 with this phase.
! Removing the avoided-crossing curvature hot spots by the pair half-sum made
! it WORSE at every threshold tried -- not the fix.  Hence the separate
! Semicl.BreakdownGamma.
! NOT included: the orbital-moment offset of the composite orbit.
! LIMIT: this is P -> 1 at EVERY junction.  A junction whose B0 is comparable
! to the field (partial transparency) needs the coherent network (L3b);
! figures/check_junctions.py measures B0 per seam cell to decide.
!=======================================================================
subroutine SemiclBreakdownAreas(Egrid, nx, ny, nb, emin, emax, numE, rcell, fname, &
                                Bflux, bfb1, mingap)

   use constants, only : pi

   integer,  intent(in) :: nx, ny, nb, numE, bfb1
   real(dp), intent(in) :: Egrid(nx,ny,nb), emin, emax, rcell(3,3), mingap
   character(*), intent(in) :: fname
   real(dp), allocatable, intent(in) :: Bflux(:,:,:)

   real(dp), allocatable :: bmn(:), bmx(:)
   integer  :: b, ie, kind, label, u, nrow, ngam, ngap
   real(dp) :: lev, dE, A_BZ, dk, area, phase, gapw, fillv, kF, relerr, gam, rmx, rmn
   logical  :: okp, haveB

   allocate(bmn(nb), bmx(nb))
   do b = 1, nb
      bmn(b) = minval(Egrid(:,:,b)) ; bmx(b) = maxval(Egrid(:,:,b))
   end do
   A_BZ = abs(rcell(1,1)*rcell(2,2) - rcell(1,2)*rcell(2,1))
   dk   = sqrt(rcell(1,1)**2 + rcell(2,1)**2)/real(nx,dp)
   dE   = (emax-emin)/real(max(numE-1,1),dp)
   haveB = allocated(Bflux)

   u = 776
   open(u, file=trim(fname), status='replace', action='write')
   write(u,'(a)') '# STRONG-BREAKDOWN composite orbits (Semicl.Breakdown, E7)'
   write(u,'(a)') '# every junction between sheets that overlap in energy taken as transparent'
   write(u,'(a)') '# (P -> 1): the orbit is the diabatic contour, and its area is the Luttinger'
   write(u,'(a)') '# count of ALL states between the anchoring full spectral gap and E.'
   write(u,'(a)') '#   iband = 1000+b : electron composite, anchored at the gap just below sheet b'
   write(u,'(a)') '#   iband = 2000+b : hole composite, anchored at the gap just above sheet b'
   write(u,'(a)') '# Between two anchors BOTH families are written -- they are the Streda fans'
   write(u,'(a)') '# of the two gaps.  The electron composite is a real closed orbit only while'
   write(u,'(a)') '# the states above its gap form closed pockets; near the upper gap it is the'
   write(u,'(a)') '# complement of a hole orbit and the hole composite is the physical one.'
   write(u,'(a)') '# Compare the signed areas of the adiabatic .FermiAreas to tell which.'
   write(u,'(a)') '# perim, npts, wx, wy are 0 (nothing is traced); sheetgap = anchor gap width.'
   write(u,'(a)') '# gamma, if present: composite Berry phase over every sheet with a contour'
   write(u,'(a)') '# at E, mod 2*pi.  The orbital-moment offset is NOT computed for composites.'
   write(u,'(a)') '# B[T] = 10475.7686 * area / (N+gamma)   [area in Ang^-2]'
   write(u,'(a,i0,a,i0)') '# grid ', nx, ' x ', ny
   write(u,'(a,es14.6,a,es14.6,a,i0)') '# E from ', emin, ' to ', emax, ' in ', numE
   ngap = 0
   do b = 1, nb-1
      rmx = maxval(bmx(1:b)) ; rmn = minval(bmn(b+1:nb))
      if (rmn - rmx > mingap .and. rmx < emax .and. rmn > emin) then
         write(u,'(a,i6,2x,i6,2x,es14.6,2x,es14.6,2x,es14.6)') '#GAP  below/above sheet ', &
              b, b+1, rmx, rmn, rmn - rmx
         ngap = ngap + 1
      end if
   end do
   if (haveB) then
      write(u,'(a)') '# iband  E[eV]        area[Ang^-2]   signed_area    perim[Ang^-1]'// &
                     '  npts  closed  wx  wy  sheetgap[eV]   n/n0          est_rel_err'// &
                     '    berryphase[rad]  gamma'
   else
      write(u,'(a)') '# iband  E[eV]        area[Ang^-2]   signed_area    perim[Ang^-1]'// &
                     '  npts  closed  wx  wy  sheetgap[eV]   n/n0          est_rel_err'
   end if
   write(u,'(a,es14.6)') '# A_BZ[Ang^-2] = ', A_BZ
   write(u,'(a,i0)') '# nbands_total = ', nb

   nrow = 0 ; ngam = 0
   do ie = 1, numE
      lev = emin + real(ie-1,dp)*dE
      ! same n/n0 convention as .FermiAreas: every band, from charge neutrality
      fillv = real(count(Egrid < lev),dp)/real(nx*ny,dp) - 0.5_dp*real(nb,dp)
      do kind = 1, 2
         call bd_composite(Egrid, nx, ny, nb, lev, bmn, bmx, mingap, emin, emax, A_BZ, &
                           Bflux, bfb1, kind, area, phase, label, gapw, okp)
         if (label == 0 .or. area == 0.0_dp) cycle
         kF = sqrt(abs(area)/pi)
         relerr = (dk/kF)**2
         if (haveB) then
            gam = -9.99e9_dp
            if (okp) then
               gam = 0.5_dp - phase/(2.0_dp*pi)
               ngam = ngam + 1
            else
               phase = -9.99e9_dp
            end if
            write(u,'(i6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,i6,2x,l1,2x,i3,2x,i3,'// &
                     '2x,es13.6,2x,es13.6,2x,es13.6,2x,es15.7,2x,es15.7)') &
                 label, lev, abs(area), area, 0.0_dp, 0, .true., 0, 0, gapw, fillv, relerr, &
                 phase, gam
         else
            write(u,'(i6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,i6,2x,l1,2x,i3,2x,i3,'// &
                     '2x,es13.6,2x,es13.6,2x,es13.6)') &
                 label, lev, abs(area), area, 0.0_dp, 0, .true., 0, 0, gapw, fillv, relerr
         end if
         nrow = nrow + 1
      end do
   end do
   close(u)

   call MIO_Print('  Semicl.Breakdown: '//trim(num2str(ngap))//' full spectral gap(s) in '// &
        'the energy window anchor the composites; wrote '//trim(num2str(nrow))// &
        ' composite orbits to '//trim(fname),'semicl')
   if (haveB) then
      call MIO_Print('  composite Berry phase for '//trim(num2str(ngam))//' of '// &
           trim(num2str(nrow))//' composites','semicl')
      if (ngam < nrow) call MIO_Print('  (the rest have a straddling sheet outside the '// &
           'flux checkpoint; widen Diag.BerryBandMin/Max)','semicl')
   end if
   deallocate(bmn, bmx)

end subroutine SemiclBreakdownAreas

!-----------------------------------------------------------------------
! One strong-breakdown composite at level `lev`.  kind = 1: electron-like,
! anchored at the nearest full gap BELOW lev; kind = 2: hole-like, anchored at
! the nearest full gap ABOVE.  Only gaps wider than mingap that intersect
! [elo, ehi] are anchors, so the caller decides which part of the spectrum is
! trusted (deep TAPW bands are basis-truncation artefacts).  label = 0 means
! no composite (lev inside a gap, or no anchor in range).  Factored out so the
! self-test exercises exactly the code the real run uses.
!-----------------------------------------------------------------------
subroutine bd_composite(Egrid, nx, ny, nb, lev, bmn, bmx, mingap, elo, ehi, A_BZ, &
                        Bflux, bfb1, kind, area, phase, label, gapw, okphase)

   integer,  intent(in)  :: nx, ny, nb, bfb1, kind
   real(dp), intent(in)  :: Egrid(nx,ny,nb), lev, bmn(nb), bmx(nb), mingap, elo, ehi, A_BZ
   real(dp), allocatable, intent(in) :: Bflux(:,:,:)
   real(dp), intent(out) :: area, phase, gapw
   integer,  intent(out) :: label
   logical,  intent(out) :: okphase

   integer  :: b, g, ga, bs1, bs2, bb, ibf, i, j, ip, jp, nin, nstr
   real(dp) :: rmx, rmn
   logical  :: inside

   area = 0.0_dp ; phase = 0.0_dp ; gapw = 0.0_dp ; label = 0 ; okphase = .false.

   ! a level inside a gap carries no orbit at all
   inside = .false.
   do b = 1, nb
      if (lev >= bmn(b) .and. lev <= bmx(b)) then
         inside = .true. ; exit
      end if
   end do
   if (.not. inside) return

   ! the gap between sheets g and g+1: the sheets are energy-sorted, so their
   ! intervals are ordered, but running extrema make the test robust anyway
   ga = 0
   if (kind == 1) then
      do g = nb-1, 1, -1
         rmx = maxval(bmx(1:g)) ; rmn = minval(bmn(g+1:nb))
         if (rmn - rmx > mingap .and. rmn <= lev .and. rmx < ehi .and. rmn > elo) then
            ga = g ; gapw = rmn - rmx ; exit
         end if
      end do
      if (ga == 0) return
      bs1 = ga + 1 ; bs2 = nb
      label = 1000 + bs1
      area = A_BZ*real(count(Egrid(:,:,bs1:bs2) < lev),dp)/real(nx*ny,dp)
   else
      do g = 1, nb-1
         rmx = maxval(bmx(1:g)) ; rmn = minval(bmn(g+1:nb))
         if (rmn - rmx > mingap .and. rmx >= lev .and. rmx < ehi .and. rmn > elo) then
            ga = g ; gapw = rmn - rmx ; exit
         end if
      end do
      if (ga == 0) return
      bs1 = 1 ; bs2 = ga
      label = 2000 + ga
      area = -A_BZ*real(count(Egrid(:,:,bs1:bs2) > lev),dp)/real(nx*ny,dp)
   end if

   if (.not. allocated(Bflux)) return
   nstr = 0
   do bb = bs1, bs2
      if (.not. (bmn(bb) < lev .and. bmx(bb) > lev)) cycle     ! no contour here
      ibf = bb - bfb1 + 1
      if (ibf < 1 .or. ibf > size(Bflux,3)) then
         phase = 0.0_dp ; return                               ! unavailable
      end if
      nstr = nstr + 1
      do j = 1, ny
         do i = 1, nx
            ip = modulo(i, nx) + 1
            jp = modulo(j, ny) + 1
            if (kind == 1) then
               nin = count([Egrid(i,j,bb), Egrid(ip,j,bb), Egrid(i,jp,bb), &
                            Egrid(ip,jp,bb)] < lev)
            else
               nin = count([Egrid(i,j,bb), Egrid(ip,j,bb), Egrid(i,jp,bb), &
                            Egrid(ip,jp,bb)] > lev)
            end if
            if (nin > 0) phase = phase + Bflux(i,j,ibf)*0.25_dp*real(nin,dp)
         end do
      end do
   end do
   okphase = nstr > 0

end subroutine bd_composite

!-----------------------------------------------------------------------
subroutine contour_perimeter(x, y, n, rcell, perim)
   integer,  intent(in)  :: n
   real(dp), intent(in)  :: x(n), y(n), rcell(3,3)
   real(dp), intent(out) :: perim
   integer  :: k, kn
   real(dp) :: ax, ay, bx, by
   perim = 0.0_dp
   do k = 1, n
      kn = k + 1 ; if (kn > n) kn = 1
      ax = x(k) *rcell(1,1) + y(k) *rcell(1,2)
      ay = x(k) *rcell(2,1) + y(k) *rcell(2,2)
      bx = x(kn)*rcell(1,1) + y(kn)*rcell(1,2)
      by = x(kn)*rcell(2,1) + y(kn)*rcell(2,2)
      perim = perim + sqrt((bx-ax)**2 + (by-ay)**2)
   end do
end subroutine contour_perimeter

!=======================================================================
! Stage 4: Onsager quantisation.  Cheap and parameter-dependent, so it is
! always recomputed from the FermiAreas checkpoint.
!=======================================================================
subroutine SemiclOnsager(fareas, numLL, gamma, outtag, nomom)

   use name, only : prefix

   character(*), intent(in) :: fareas
   integer,      intent(in) :: numLL
   real(dp),     intent(in) :: gamma
   ! E7: the strong-breakdown composite orbits go through this same stage.
   ! outtag is appended to the output names (.LLfan<tag>, .Wannier<tag>) so the
   ! adiabatic fan is never overwritten; nomom switches the orbital-moment
   ! offset off for a checkpoint that does not carry it.  Both absent = the
   ! original behaviour, byte for byte.
   character(*), intent(in), optional :: outtag
   logical,      intent(in), optional :: nomom
   character(len=16) :: tag

   integer  :: u, uf, uw, ios, b, np, wx, wy, nll, nwrit
   real(dp) :: lev, area, sarea, perim, gap, Bfield, nn0, B1, rerr, relerr_max
   character(len=400) :: line
   logical  :: cl
   real(dp) :: phaseB, gamOrb, g, minNg
   integer  :: nanom, nogam, nneg
   logical  :: perOrbit
   real(dp) :: delm, dAdE_line
   integer  :: nodelm
   logical  :: useMom
   ! E5: fold in a second valley's checkpoint and carry the degeneracy on n/n0
   character(len=200) :: fareas2
   real(dp) :: degen
   integer  :: ival, nval
   logical  :: have2
   ! L1c/L2: second-order term
   real(dp) :: psi, psied, psiLP, dchi
   logical  :: useChi
   integer  :: nodchi, nnochi

   call MIO_InputParameter('Semicl.MaxRelErr', relerr_max, 0.10_dp)
   ! levels with |N+gamma| below this are the anomalous one, not a field level
   call MIO_InputParameter('Semicl.MinNplusGamma', minNg, 0.10_dp)
   ! Semicl.Gamma < 0 means "use the per-orbit gamma computed from the Berry
   ! phase" (the last column of .FermiAreas, written when Semicl.BerryPhase).
   perOrbit = gamma < 0.0_dp
   ! L1: add the per-orbit orbital-moment offset to N+gamma.  Read here rather
   ! than passed down because this stage runs off the .FermiAreas checkpoint
   ! alone, exactly like Semicl.MaxRelErr and Semicl.Gamma.
   call MIO_InputParameter('Semicl.OrbitalMoment', useMom, .false.)
   ! L1c/L2: the per-orbit second-order offset dchi*B, from the same checkpoint
   call MIO_InputParameter('Semicl.Susceptibility', useChi, .false.)
   ! E5: a run is one valley, and at B=0 the two valleys have IDENTICAL areas
   ! (time reversal), so the second valley adds no new orbit -- what it adds is
   ! the opposite Berry phase, i.e. which band edge carries the anomalous
   ! N+gamma=0 level.  Point Semicl.SecondValleyAreas at the other run's
   ! .FermiAreas to fold it in; the valley is tagged in a new column.
   call MIO_InputParameter('Semicl.SecondValleyAreas', fareas2, ' ')
   ! n/n0 is written per valley per spin.  Semicl.Degeneracy multiplies it, so
   ! 4 = 2 spins x 2 valleys gives the filling a measurement would see.
   call MIO_InputParameter('Semicl.Degeneracy', degen, 1.0_dp)
   tag = ' '
   if (present(outtag)) tag = outtag
   if (present(nomom)) then
      if (nomom) useMom = .false.
      if (nomom) useChi = .false.
   end if
   ! the second-valley merge reads a .FermiAreas file; a tagged (breakdown)
   ! checkpoint is a different quantity and is quantised on its own
   if (len_trim(tag) > 0) fareas2 = ' '
   have2 = len_trim(fareas2) > 0
   nval = 1
   if (have2) then
      inquire(file=trim(fareas2), exist=have2)
      if (have2) then
         nval = 2
      else
         call MIO_Print('  Semicl.SecondValleyAreas not found: '//trim(fareas2),'semicl')
      end if
   end if
   uf = 774
   open(uf, file=trim(prefix)//'.LLfan'//trim(tag), status='replace', action='write')
   write(uf,'(a)') '# semiclassical Landau fan'
   if (len_trim(tag) > 0) write(uf,'(a)') '# STRONG-BREAKDOWN composite orbits (Semicl.Breakdown)'
   write(uf,'(a,f8.4)') '# gamma = ', gamma
   if (useMom) write(uf,'(a)') '# N+gamma includes the per-orbit orbital-moment '// &
        'offset delta_m (L1)'
   if (useChi) write(uf,'(a)') '# B solves SEMICL_C*area = B*(N+gamma+delta_m) + dchi*B^2 '// &
        '(second-order term, L1c/L2)'
   if (degen /= 1.0_dp) write(uf,'(a,f8.4)') '# degeneracy on n/n0 = ', degen
   if (nval == 2) then
      write(uf,'(a)') '# iband  E[eV]        N   B[T]           area[Ang^-2]  sheetgap[eV]  ivalley'
   else
      write(uf,'(a)') '# iband  E[eV]        N   B[T]           area[Ang^-2]  sheetgap[eV]'
   end if

   ! Wannier / Streda form: in these axes the LL branches are STRAIGHT lines
   !   n/n0 = t*(phi/phi0) + s   with integer t (Chern) and s (band filling).
   ! n/n0 is per moire cell per valley per spin, measured from charge
   ! neutrality; phi/phi0 = B*A_moire/Phi_0, and Phi_0/A_moire = B1 below.
   uw = 775
   open(uw, file=trim(prefix)//'.Wannier'//trim(tag), status='replace', action='write')
   write(uw,'(a)') '# semiclassical Wannier diagram (straight-line LL branches)'
   if (len_trim(tag) > 0) write(uw,'(a)') '# STRONG-BREAKDOWN composite orbits (Semicl.Breakdown)'
   write(uw,'(a,f8.4)') '# gamma = ', gamma
   if (useMom) write(uw,'(a)') '# N+gamma includes the per-orbit orbital-moment '// &
        'offset delta_m (L1)'
   if (useChi) write(uw,'(a)') '# B solves SEMICL_C*area = B*(N+gamma+delta_m) + dchi*B^2 '// &
        '(second-order term, L1c/L2)'
   if (nval == 2) then
      write(uw,'(a)') '# iband  n/n0          phi/phi0       B[T]           E[eV]         N   ivalley'
   else
      write(uw,'(a)') '# iband  n/n0          phi/phi0       B[T]           E[eV]         N'
   end if

   ! Phi_0/A_moire in tesla:  A_moire = (2*pi)^2/A_BZ, so B1 = Phi0*A_BZ/(2pi)^2
   B1 = SEMICL_C * A_BZ_from_header(fareas)
   write(uw,'(a,es14.6,a)') '# Phi0/A_moire = ', B1, ' T   (phi/phi0 = B/this)'

   nwrit = 0 ; nanom = 0 ; nogam = 0 ; nodelm = 0 ; nneg = 0
   nodchi = 0 ; nnochi = 0
   u = 773
   do ival = 1, nval
   if (ival == 1) then
      open(u, file=trim(fareas), status='old', action='read', iostat=ios)
   else
      open(u, file=trim(fareas2), status='old', action='read', iostat=ios)
   end if
   if (ios /= 0) then
      call MIO_Print('ERROR: cannot read the areas file for valley '// &
           trim(num2str(ival)),'semicl')
      cycle
   end if
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      gamOrb = -9.99e9_dp ; delm = -9.99e9_dp ; dchi = -9.99e9_dp
      read(line,*,iostat=ios) b, lev, area, sarea, perim, np, cl, wx, wy, gap, nn0, rerr, &
           phaseB, gamOrb, delm, dAdE_line, psi, psied, psiLP, dchi
      if (ios /= 0) then
         gamOrb = -9.99e9_dp ; delm = -9.99e9_dp ; dchi = -9.99e9_dp
         read(line,*,iostat=ios) b, lev, area, sarea, perim, np, cl, wx, wy, gap, nn0, rerr, &
              phaseB, gamOrb, delm, dAdE_line
      end if
      if (ios /= 0) then
         gamOrb = -9.99e9_dp ; delm = -9.99e9_dp
         read(line,*,iostat=ios) b, lev, area, sarea, perim, np, cl, wx, wy, gap, nn0, rerr, &
              phaseB, gamOrb
      end if
      if (ios /= 0) then
         gamOrb = -9.99e9_dp ; delm = -9.99e9_dp
         read(line,*,iostat=ios) b, lev, area, sarea, perim, np, cl, wx, wy, gap, nn0, rerr
      end if
      if (ios /= 0) cycle
      ! drop orbits whose polygon error exceeds the tolerance
      if (rerr > relerr_max) cycle
      if (area <= 0.0_dp) cycle              ! open orbit or degenerate
      g = gamma
      if (perOrbit) then
         if (gamOrb < -1.0e8_dp) then
            nogam = nogam + 1
            cycle                            ! no Berry phase for this orbit
         end if
         ! Onsager's gamma is defined MOD 1: A*l_B^2 = 2*pi*(N+gamma) with N a
         ! non-negative integer, so the integer part of gamma is absorbed into
         ! N.  It matters here because Phi_B is NOT bounded by pi -- a large
         ! orbit can enclose several curvature hotspots (measured up to
         ! Phi_B = 3.32*pi, i.e. gamma = -1.16), and using that raw value gives
         ! N+gamma < 0 and a NEGATIVE field.  The raw value stays in the
         ! .FermiAreas column, since that is the physical Berry phase and the
         ! valley relation gamma_K + gamma_K' = 1 is stated on it.
         g = gamOrb - floor(gamOrb)
      end if
      if (useMom) then
         ! delta_m is a property of the orbit, like gamma.  If the checkpoint
         ! does not carry it, the correction the user asked for cannot be
         ! applied to this orbit, so the orbit is skipped and counted rather
         ! than quietly quantised at the old, incomplete order.
         if (delm < -1.0e8_dp) then
            nodelm = nodelm + 1
            cycle
         end if
         g = g + delm
      end if
      ! same policy for the second-order offset: no silent fallback to first order
      if (useChi) then
         if (dchi < -1.0e8_dp) then
            nodchi = nodchi + 1
            cycle
         end if
      end if
      do nll = 0, numLL-1
         ! N + gamma -> 0 is the ANOMALOUS level: the orbit area vanishes, so no
         ! finite B satisfies Onsager and the level is pinned to the band edge,
         ! B-independent (for a Dirac sheet, gamma -> 0 puts it exactly there).
         ! Emitting SEMICL_C*area/0 = Infinity here was a real defect.
         !
         ! The cut is on |N+gamma| < Semicl.MinNplusGamma, not on an exact zero:
         ! a per-orbit gamma comes out at 0.0068 rather than 0, and N=0 then
         ! gives B = C*A/0.0068, a level 150x too high in field that is pure
         ! artefact.  Onsager is asymptotic in N anyway, so N+gamma well below 1
         ! carries no information.  Found by the Phase 2 metric: leaving these
         ! in moved the median residual from 0.81 to 1.37 meV and to 7.4 meV in
         ! the 8-10 T bin.
         ! gamma is taken mod 1 BEFORE delta_m is added, so a large negative
         ! orbital-moment offset can push N+gamma+delta_m below zero.  That is
         ! not a level (it would be written at B < 0); the same physical level
         ! appears at the next N.  Found on the sheet-54 hole orbit near the
         ! hSDP gap: 11 of 2289 levels at 48x48.  Without the moment
         ! N+gamma >= 0 always, so this never fires and the output is unchanged.
         if (real(nll,dp) + g <= -minNg) then
            nneg = nneg + 1
            cycle
         end if
         if (abs(real(nll,dp) + g) < minNg) then
            nanom = nanom + 1
            cycle
         end if
         if (useChi) then
            ! delta_chi = dchi*B is linear in B, so Onsager becomes a quadratic
            Bfield = semicl_onsager_B(SEMICL_C*area, real(nll,dp) + g, dchi)
            if (Bfield <= 0.0_dp) then
               nnochi = nnochi + 1
               cycle
            end if
         else
            Bfield = SEMICL_C * area / (real(nll,dp) + g)
         end if
         ! the valley column appears only when two valleys are merged, so a
         ! single-valley run keeps exactly the old file format
         if (nval == 2) then
            write(uf,'(i6,2x,es13.6,2x,i4,2x,es13.6,2x,es13.6,2x,es13.6,2x,i4)') &
                 b, lev, nll, Bfield, area, gap, ival
            write(uw,'(i6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,i4,2x,i4)') &
                 b, degen*nn0, Bfield/B1, Bfield, lev, nll, ival
         else
            write(uf,'(i6,2x,es13.6,2x,i4,2x,es13.6,2x,es13.6,2x,es13.6)') &
                 b, lev, nll, Bfield, area, gap
            write(uw,'(i6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,i4)') &
                 b, degen*nn0, Bfield/B1, Bfield, lev, nll
         end if
         nwrit = nwrit + 1
      end do
   end do
   close(u)
   end do
   close(uf) ; close(uw)
   call MIO_Print('  wrote '//trim(prefix)//'.LLfan'//trim(tag)//' and .Wannier'// &
        trim(tag)//': '// &
        trim(num2str(nwrit))//' points (Semicl.MaxRelErr = '// &
        trim(num2str(relerr_max,3))//')','semicl')
   if (perOrbit) then
      call MIO_Print('  gamma taken per orbit from the Berry phase','semicl')
      if (nogam > 0) call MIO_Print('  '//trim(num2str(nogam))// &
           ' orbits had no per-orbit gamma and were skipped','semicl')
   end if
   if (useMom) then
      call MIO_Print('  orbital-moment offset delta_m added to N+gamma '// &
           '(first-order semiclassical quantisation)','semicl')
      if (nodelm > 0) call MIO_Print('  '//trim(num2str(nodelm))// &
           ' orbits had no delta_m in the checkpoint and were skipped','semicl')
   end if
   if (useChi) then
      call MIO_Print('  second-order offset dchi*B included (quadratic Onsager, L1c/L2)','semicl')
      if (nodchi > 0) call MIO_Print('  '//trim(num2str(nodchi))// &
           ' orbits had no dchi in the checkpoint and were skipped','semicl')
      if (nnochi > 0) call MIO_Print('  '//trim(num2str(nnochi))//' levels have no '// &
           'positive root of the quadratic (second order not a small correction) '// &
           'and were dropped','semicl')
   end if
   if (nval == 2) call MIO_Print('  merged two valleys (Semicl.SecondValleyAreas); '// &
        'n/n0 multiplied by Semicl.Degeneracy = '//trim(num2str(degen,3)),'semicl')
   if (nneg > 0) call MIO_Print('  '//trim(num2str(nneg))//' levels had '// &
        'N+gamma+delta_m < 0 (a large negative orbital-moment offset) and were '// &
        'dropped: no B < 0 is written; the level appears at the next N','semicl')
   if (nanom > 0) call MIO_Print('  '//trim(num2str(nanom))//' levels have '// &
        '|N+gamma| < '//trim(num2str(minNg,3))//': the anomalous level, pinned to '// &
        'the band edge and B-independent, so no B is written for it','semicl')

end subroutine SemiclOnsager

!-----------------------------------------------------------------------
real(dp) function A_BZ_from_header(fname)
   character(*), intent(in) :: fname
   integer :: u, ios
   character(len=400) :: line
   A_BZ_from_header = 0.0_dp
   u = 776
   open(u, file=trim(fname), status='old', action='read', iostat=ios)
   if (ios /= 0) return
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (index(line,'# A_BZ[Ang^-2] =') == 1) then
         read(line(index(line,'=')+1:),*,iostat=ios) A_BZ_from_header
         exit
      end if
   end do
   close(u)
end function A_BZ_from_header

!=======================================================================
! Cyclotron mass and SdH frequency from the contour areas.
!
!    m*(E) = (hbar^2/2*pi) * dA/dE          F[T] = (hbar/2*pi*e) * A
!
! Both are pure functions of A(E), so this reads back <prefix>.FermiAreas
! and needs no eigenvalue work.  Note F = SEMICL_C*A is exactly the field of
! the (N+gamma)=1 Landau level, so .Masses and .LLfan are consistent by
! construction -- the F column is a convenience, not an independent check.
! The check that carries information is m*: for a Dirac sheet it must rise
! linearly with |E-E_D| with slope 1/v_F^2 taken from the SAME areas.
!
! dA/dE is differenced ALONG ONE ORBIT.  The number of contours on a sheet
! changes at a Lifshitz transition and differencing across one mixes two
! different orbits, returning a divergent, meaningless mass.  Each band's
! energy series is therefore cut into segments of constant valid-orbit count
! AND contiguous energy; orbits are matched between neighbouring energies by
! rank in |area|; and only symmetric (central) differences are written, since
! a one-sided difference at a segment end sits exactly where the topology is
! about to change and is the least trustworthy point there is.
!
! Rank matching is exact wherever a sheet carries a single contour -- the
! whole primary miniband window for G/hBN CellSize 55.  With several pockets
! it assumes their areas do not cross inside a segment.
!
! A hole orbit shrinks as E rises, so its m* comes out NEGATIVE.  The `hole`
! column is taken independently from the sign of the signed area, so a
! disagreement between sign(m*) and `hole` shows up as an inconsistency
! instead of being averaged away.
!=======================================================================
subroutine SemiclMasses(fareas)

   use name, only : prefix

   character(*), intent(in) :: fareas

   integer, parameter :: MAXORB = 32     ! orbits tracked per (band,energy)

   character(len=400) :: line
   integer  :: u, um, ios, b, np, wx, wy, nwrit, nover
   real(dp) :: lev, area, sarea, perim, gap, nn0, rerr
   logical  :: cl

   integer  :: nemax, nlin, bcur, ne
   integer,  allocatable :: nva(:)
   real(dp), allocatable :: Ea(:), Aa(:,:), Ra(:,:), Ga(:)
   logical,  allocatable :: Ha(:,:)

   u = 777 ; um = 778

   ! ---- pass 1: longest per-band block, to size the work arrays --------
   open(u, file=trim(fareas), status='old', action='read', iostat=ios)
   if (ios /= 0) then
      call MIO_Print('ERROR: cannot read '//trim(fareas),'semicl')
      return
   end if
   nemax = 0 ; nlin = 0 ; bcur = -huge(1)
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      read(line,*,iostat=ios) b, lev
      if (ios /= 0) cycle
      if (b /= bcur) then
         nemax = max(nemax, nlin) ; nlin = 0 ; bcur = b
      end if
      nlin = nlin + 1
   end do
   nemax = max(nemax, nlin)
   close(u)

   if (nemax < 3) then
      call MIO_Print('  fewer than 3 contour energies per band; '// &
           'no .Masses written','semicl')
      return
   end if

   allocate(Ea(nemax), Aa(MAXORB,nemax), Ra(MAXORB,nemax), Ga(nemax))
   allocate(Ha(MAXORB,nemax), nva(nemax))

   open(um, file=trim(prefix)//'.Masses', status='replace', action='write')
   write(um,'(a)') '# cyclotron mass and SdH frequency per orbit'
   write(um,'(a,es14.6)') '# m*/m_e = SEMICL_M * dA/dE,  SEMICL_M = ', SEMICL_M
   write(um,'(a,f12.4)')  '# F[T]   = SEMICL_C * A,      SEMICL_C = ', SEMICL_C
   write(um,'(a)') '# dA/dE is a CENTRAL difference along one orbit, inside a'
   write(um,'(a)') '# segment of constant orbit count (no differencing across a'
   write(um,'(a)') '# Lifshitz transition).  Orbits are ranked by |area|, iorb=1'
   write(um,'(a)') '# the largest.  m* < 0 is a hole orbit and must agree with'
   write(um,'(a)') '# the `hole` column, which comes from the signed area.'
   write(um,'(a)') '# rel_err = worst est_rel_err of the three points differenced.'
   write(um,'(a)') '# iband  E[eV]        A[Ang^-2]      dA/dE          m*/m_e'// &
                   '         F[T]           iorb  norb  hole  rel_err        sheetgap[eV]'

   ! ---- pass 2: collect per band, flush when the band index changes ----
   open(u, file=trim(fareas), status='old', action='read', iostat=ios)
   if (ios /= 0) then
      call MIO_Print('ERROR: cannot re-read '//trim(fareas),'semicl')
      close(um) ; return
   end if

   bcur = -huge(1) ; ne = 0 ; nwrit = 0 ; nover = 0
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      read(line,*,iostat=ios) b, lev, area, sarea, perim, np, cl, wx, wy, &
                              gap, nn0, rerr
      if (ios /= 0) cycle

      if (b /= bcur) then
         if (ne > 0) call SemiclMassesBand(um, bcur, ne, Ea, Aa, Ra, Ha, Ga, &
                                           nva, MAXORB, nwrit)
         bcur = b ; ne = 0
      end if

      ! area > 0 marks a genuine closed, non-winding orbit (see the writer):
      ! open and truncated traces carry 0.0 and have no cyclotron area
      if (area <= 0.0_dp) cycle

      if (ne == 0) then
         ne = 1
         Ea(1) = lev ; nva(1) = 0 ; Ga(1) = gap
      else if (lev > Ea(ne)) then
         if (ne >= nemax) cycle
         ne = ne + 1
         Ea(ne) = lev ; nva(ne) = 0 ; Ga(ne) = gap
      end if

      call SemiclOrbitInsert(Aa(:,ne), Ra(:,ne), Ha(:,ne), nva(ne), MAXORB, &
                             area, rerr, sarea < 0.0_dp)
      if (nva(ne) > MAXORB) nover = nover + 1
   end do
   if (ne > 0) call SemiclMassesBand(um, bcur, ne, Ea, Aa, Ra, Ha, Ga, &
                                     nva, MAXORB, nwrit)

   close(u) ; close(um)
   deallocate(Ea, Aa, Ra, Ga, Ha, nva)

   call MIO_Print('  wrote '//trim(prefix)//'.Masses: '// &
        trim(num2str(nwrit))//' (orbit,energy) points','semicl')
   if (nover > 0) call MIO_Print('  NOTE: '//trim(num2str(nover))// &
        ' orbits beyond MAXORB were dropped from the ranking; masses for '// &
        'those energies are unreliable','semicl')

end subroutine SemiclMasses

!-----------------------------------------------------------------------
! Insert one orbit into the per-energy list, kept sorted by DESCENDING
! |area| so orbits can be matched between neighbouring energies by rank.
! The stored list is capped at mx but nv counts EVERY valid orbit, so a
! clamped energy is visible as nv > mx (and breaks the segment) rather than
! silently renumbering the orbits.
!-----------------------------------------------------------------------
subroutine SemiclOrbitInsert(a, r, h, nv, mx, area, rerr, hole)

   integer,  intent(in)    :: mx
   real(dp), intent(inout) :: a(mx), r(mx)
   logical,  intent(inout) :: h(mx)
   integer,  intent(inout) :: nv
   real(dp), intent(in)    :: area, rerr
   logical,  intent(in)    :: hole

   integer :: i, n, j

   n = min(nv, mx)
   i = 1
   do while (i <= n)
      if (area > a(i)) exit
      i = i + 1
   end do
   nv = nv + 1
   if (i > mx) return                      ! smaller than everything kept
   do j = min(n+1, mx), i+1, -1
      a(j) = a(j-1) ; r(j) = r(j-1) ; h(j) = h(j-1)
   end do
   a(i) = area ; r(i) = rerr ; h(i) = hole

end subroutine SemiclOrbitInsert

!-----------------------------------------------------------------------
! One band: cut the energy series into segments of constant orbit count and
! contiguous energy, then central-difference each orbit inside a segment.
!-----------------------------------------------------------------------
subroutine SemiclMassesBand(um, b, ne, Ea, Aa, Ra, Ha, Ga, nva, mx, nwrit)

   integer,  intent(in)    :: um, b, ne, mx
   real(dp), intent(in)    :: Ea(ne), Aa(mx,ne), Ra(mx,ne), Ga(ne)
   logical,  intent(in)    :: Ha(mx,ne)
   integer,  intent(in)    :: nva(ne)
   integer,  intent(inout) :: nwrit

   integer  :: i, i1, i2, j, k, nk
   real(dp) :: dE0, dAdE, mst, Fsdh, rr
   logical  :: brk

   if (ne < 3) return

   ! The energy grid is uniform, so the SMALLEST step present is its spacing;
   ! a larger step means energies where this sheet carried no valid orbit and
   ! must not be differenced across.
   dE0 = huge(1.0_dp)
   do i = 2, ne
      if (Ea(i)-Ea(i-1) > 0.0_dp) dE0 = min(dE0, Ea(i)-Ea(i-1))
   end do
   if (dE0 >= huge(1.0_dp)) return

   i1 = 1
   do i = 2, ne+1
      if (i > ne) then
         brk = .true.
      else
         brk = (nva(i) /= nva(i1)) .or. (Ea(i)-Ea(i-1) > 1.5_dp*dE0)
      end if
      if (.not. brk) cycle

      i2 = i - 1
      if (i2-i1 >= 2) then
         nk = min(nva(i1), mx)
         do k = 1, nk
            do j = i1+1, i2-1
               dAdE = (Aa(k,j+1) - Aa(k,j-1))/(Ea(j+1) - Ea(j-1))
               mst  = SEMICL_M * dAdE
               Fsdh = SEMICL_C * Aa(k,j)
               rr   = max(Ra(k,j-1), Ra(k,j), Ra(k,j+1))
               write(um,'(i6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,'// &
                         'es13.6,2x,i5,2x,i5,2x,l1,2x,es13.6,2x,es13.6)') &
                    b, Ea(j), Aa(k,j), dAdE, mst, Fsdh, k, nva(j), Ha(k,j), &
                    rr, Ga(j)
               nwrit = nwrit + 1
            end do
         end do
      end if
      i1 = i
   end do

end subroutine SemiclMassesBand

!=======================================================================
! L5 -- path-ordered (non-Abelian) Wilson loop of a closed orbit.
!
! The Abelian Berry phase of an orbit is a plain SUM of plaquette fluxes over
! its interior (SemiclOrbitBerry), because Stokes' theorem abelianises for a
! U(1) connection.  For a multi-band (U(N)) connection there is no such
! formula: Stokes is surface-ordered, and only det W survives -- which is
! exactly the Abelian answer again, and so uninformative.  The eigenvalues of
! the holonomy do NOT abelianise, and getting them needs the ORDERED traversal
! of the orbit itself:
!
!     W = prod_s  M_s ,   M_s = <u_a(k_s) | u_b(k_{s+1})>  restricted to the
!                               multiplet selected at each end,
!
! walked once around a closed lattice loop that follows the contour.  The
! eigenvalues of W are exp(i*theta_a); arg(det W) reproduces the Abelian
! composite phase, so the theta_a are the new information.
!
! THE CRUX (L5_PROPOSAL.md 3.3): the multiplet must be selected POINTWISE on
! the path -- all states within Semicl.WilsonWindow of the orbit energy at THAT
! k -- and that set must be gapped from the rest at EVERY point.  Selecting it
! globally ("the sheets that carry a contour somewhere") is not adiabatically
! closed along the orbit: states enter and leave the selected space, |det W|
! falls below 1, and the number then measures the leak rather than convergence.
! With pointwise selection |det W| becomes a genuine check, which is why it is
! a gate here and not a reported diagnostic.
!=======================================================================

!-----------------------------------------------------------------------
! Reads <prefix>.BerryLinks (diag.F90, BerryLinksWrite) into
! Lnk(nbnd,nbnd,nk,2):   Lnk(m,n,ip,mu) = <u_m(k_ip) | u_n(k_ip + e_mu)>
! mu = 1 steps in i, mu = 2 steps in j; ip = (j-1)*nx + i; m,n are 1-based
! INSIDE the window, so window index m is spectrum band b1 + m - 1.
!-----------------------------------------------------------------------
subroutine SemiclReadBerryLinks(fname, nx, ny, Lnk, b1, b2, ok)

   character(*), intent(in)  :: fname
   integer,      intent(in)  :: nx, ny
   complex(dp), allocatable, intent(out) :: Lnk(:,:,:,:)
   integer,      intent(out) :: b1, b2
   logical,      intent(out) :: ok

   integer  :: u, ios, gx, gy, mu, i, j, m, n, nread, nbnd
   real(dp) :: re, im
   character(len=400) :: line

   ok = .false. ; b1 = 0 ; b2 = -1 ; gx = -1 ; gy = -1
   u = 776
   open(u, file=trim(fname), status='old', action='read', iostat=ios)
   if (ios /= 0) then
      call MIO_Print('  Semicl.Wilson: cannot open '//trim(fname),'semicl') ; return
   end if
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) /= '#') exit
      if (index(line,'# grid ') == 1)  read(line(8:),*,iostat=ios) gx, gy
      if (index(line,'# bands ') == 1) read(line(9:),*,iostat=ios) b1, b2
   end do
   if (gx /= nx .or. gy /= ny) then
      call MIO_Print('  Semicl.Wilson: '//trim(fname)//' is '//trim(num2str(gx))// &
           ' x '//trim(num2str(gy))//', the bands grid is '//trim(num2str(nx))// &
           ' x '//trim(num2str(ny)),'semicl')
      close(u) ; return
   end if
   if (b2 < b1) then
      call MIO_Print('  Semicl.Wilson: no band window in '//trim(fname),'semicl')
      close(u) ; return
   end if
   nbnd = b2 - b1 + 1
   allocate(Lnk(nbnd, nbnd, nx*ny, 2))
   Lnk = cmplx(0.0_dp, 0.0_dp, kind=dp)
   rewind(u)
   nread = 0
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
      read(line,*,iostat=ios) mu, i, j, m, n, re, im
      if (ios /= 0) cycle
      if (mu < 1 .or. mu > 2) cycle
      if (i < 1 .or. i > nx .or. j < 1 .or. j > ny) cycle
      if (m < 1 .or. m > nbnd .or. n < 1 .or. n > nbnd) cycle
      Lnk(m, n, (j-1)*nx + i, mu) = cmplx(re, im, kind=dp)
      nread = nread + 1
   end do
   close(u)
   call MIO_Print('  Semicl.Wilson: read '//trim(num2str(nread))//' of '// &
        trim(num2str(2*nbnd*nbnd*nx*ny))//' link elements from '//trim(fname),'semicl')
   ok = nread == 2*nbnd*nbnd*nx*ny

end subroutine SemiclReadBerryLinks

!-----------------------------------------------------------------------
! A closed lattice loop of nearest-neighbour steps following the contour.
! Contour points are on cell EDGES (marching squares), so they are snapped to
! the nearest node and consecutive nodes are joined by a staircase, i in
! first.  Steps are taken by minimum image, so the loop closes on the torus.
! pth(1,s), pth(2,s) = node; stp(s) = +/-1 for a step in i, +/-2 for one in j.
!-----------------------------------------------------------------------
subroutine wilson_path(xs, ys, np, nx, ny, pth, stp, nstep, ok)

   integer,  intent(in)  :: np, nx, ny
   real(dp), intent(in)  :: xs(np), ys(np)
   integer,  intent(out) :: pth(2, 8*(nx+ny)), stp(8*(nx+ny)), nstep
   logical,  intent(out) :: ok

   integer :: p, ci, cj, ti, tj, di, dj, s, maxs, i0, j0

   ok = .false. ; nstep = 0 ; maxs = 8*(nx+ny)
   ci = modulo(nint(xs(1)*real(nx,dp)), nx) + 1
   cj = modulo(nint(ys(1)*real(ny,dp)), ny) + 1
   i0 = ci ; j0 = cj
   do p = 2, np + 1
      if (p <= np) then
         ti = modulo(nint(xs(p)*real(nx,dp)), nx) + 1
         tj = modulo(nint(ys(p)*real(ny,dp)), ny) + 1
      else
         ti = i0 ; tj = j0                      ! close the loop
      end if
      ! minimum-image displacement to the target
      di = ti - ci ; if (di >  nx/2) di = di - nx ; if (di < -nx/2) di = di + nx
      dj = tj - cj ; if (dj >  ny/2) dj = dj - ny ; if (dj < -ny/2) dj = dj + ny
      do s = 1, abs(di)
         nstep = nstep + 1 ; if (nstep > maxs) return
         pth(1,nstep) = ci ; pth(2,nstep) = cj
         stp(nstep) = sign(1, di)
         ci = iwrap(ci + sign(1, di), nx)
      end do
      do s = 1, abs(dj)
         nstep = nstep + 1 ; if (nstep > maxs) return
         pth(1,nstep) = ci ; pth(2,nstep) = cj
         stp(nstep) = sign(2, dj)
         cj = iwrap(cj + sign(1, dj), ny)
      end do
   end do
   ok = nstep >= 4 .and. ci == i0 .and. cj == j0

end subroutine wilson_path

!-----------------------------------------------------------------------
! Holonomy of one closed orbit.  Ewin(:,:,m) = energy of window band m.
! Returns nsel (multiplet size), |det W|, the nsel eigen-phases [rad], the
! worst pointwise isolation along the path, and ok = passed every gate.
! ok = .false. means the orbit is REJECTED, not that its numbers are poor:
! a rejected orbit must not be reported (L5_PROPOSAL 3.3).
!-----------------------------------------------------------------------
subroutine SemiclOrbitWilson(Ewin, Lnk, nx, ny, nbnd, xs, ys, np, lev, win, isotol, &
                             mindet, detW, phases, nsel, nstep, isomin, ok)

   integer,     intent(in)  :: nx, ny, nbnd, np
   real(dp),    intent(in)  :: Ewin(nx,ny,nbnd), xs(np), ys(np), lev, win, isotol, mindet
   complex(dp), intent(in)  :: Lnk(nbnd,nbnd,nx*ny,2)
   real(dp),    intent(out) :: detW, phases(nbnd), isomin
   integer,     intent(out) :: nsel, nstep
   logical,     intent(out) :: ok

   integer,     allocatable :: pth(:,:), stp(:), sel(:,:), ns(:)
   complex(dp), allocatable :: W(:,:), M(:,:), Mf(:,:), Wtmp(:,:), ev(:), work(:), vdum(:,:)
   real(dp),    allocatable :: rwork(:)
   integer  :: maxs, s, sn, a, bq, mb, ipc, info, lwork, i2, j2
   real(dp) :: dsel, dout, iso
   complex(dp) :: dete

   ok = .false. ; detW = 0.0_dp ; nsel = 0 ; nstep = 0 ; isomin = 0.0_dp
   phases = -9.99e9_dp
   maxs = 8*(nx+ny)
   allocate(pth(2,maxs), stp(maxs))
   call wilson_path(xs, ys, np, nx, ny, pth, stp, nstep, ok)
   if (.not. ok) then
      deallocate(pth, stp) ; ok = .false. ; return
   end if

   ! ---- POINTWISE multiplet selection, and its isolation, at every node ----
   allocate(sel(nbnd, nstep), ns(nstep))
   isomin = huge(1.0_dp)
   do s = 1, nstep
      ns(s) = 0 ; dsel = 0.0_dp ; dout = huge(1.0_dp)
      do mb = 1, nbnd
         if (abs(Ewin(pth(1,s), pth(2,s), mb) - lev) <= win) then
            ns(s) = ns(s) + 1 ; sel(ns(s), s) = mb
            dsel = max(dsel, abs(Ewin(pth(1,s), pth(2,s), mb) - lev))
         else
            dout = min(dout, abs(Ewin(pth(1,s), pth(2,s), mb) - lev))
         end if
      end do
      iso = dout - dsel                   ! margin between the window and the rest
      isomin = min(isomin, iso)
      if (ns(s) == 0) then
         deallocate(pth, stp, sel, ns) ; ok = .false. ; return
      end if
   end do
   ! the multiplet must not change size, and must stay gapped, anywhere
   if (any(ns /= ns(1)) .or. isomin < isotol) then
      deallocate(pth, stp, sel, ns) ; ok = .false. ; return
   end if
   nsel = ns(1)

   ! ---- path-ordered product ----
   allocate(W(nsel,nsel), M(nsel,nsel), Wtmp(nsel,nsel), Mf(nbnd,nbnd))
   W = cmplx(0.0_dp, 0.0_dp, kind=dp)
   do a = 1, nsel
      W(a,a) = cmplx(1.0_dp, 0.0_dp, kind=dp)
   end do
   do s = 1, nstep
      sn = s + 1 ; if (sn > nstep) sn = 1
      select case (stp(s))
      case ( 1) ; Mf = Lnk(:,:,(pth(2,s)-1)*nx + pth(1,s), 1)
      case ( 2) ; Mf = Lnk(:,:,(pth(2,s)-1)*nx + pth(1,s), 2)
      case (-1)                              ! reverse of the forward link at the TARGET
         i2 = pth(1,sn) ; j2 = pth(2,sn)
         Mf = conjg(transpose(Lnk(:,:,(j2-1)*nx + i2, 1)))
      case (-2)
         i2 = pth(1,sn) ; j2 = pth(2,sn)
         Mf = conjg(transpose(Lnk(:,:,(j2-1)*nx + i2, 2)))
      end select
      do bq = 1, nsel
         do a = 1, nsel
            M(a,bq) = Mf(sel(a,s), sel(bq,sn))
         end do
      end do
      Wtmp = matmul(W, M)
      W = Wtmp
   end do

   ! ---- eigenvalues of the holonomy ----
   lwork = 8*nsel
   allocate(ev(nsel), work(lwork), rwork(2*nsel), vdum(1,1))
   call zgeev('N', 'N', nsel, W, nsel, ev, vdum, 1, vdum, 1, work, lwork, rwork, info)
   if (info /= 0) then
      deallocate(pth, stp, sel, ns, W, M, Wtmp, Mf, ev, work, rwork, vdum)
      ok = .false. ; return
   end if
   dete = cmplx(1.0_dp, 0.0_dp, kind=dp)
   do a = 1, nsel
      dete = dete*ev(a)
      phases(a) = atan2(aimag(ev(a)), real(ev(a), dp))
   end do
   detW = abs(dete)
   ! a unitary connection's holonomy has |det W| = 1; below the gate the
   ! selected subspace did not close along the loop
   ok = detW >= mindet
   deallocate(pth, stp, sel, ns, W, M, Wtmp, Mf, ev, work, rwork, vdum)

end subroutine SemiclOrbitWilson

!=======================================================================
! Driver: path-ordered Wilson loop on every closed orbit -> <prefix>.Wilson
! Re-traces the contours with the SAME marching-squares routine the areas
! stage uses, rather than widening SemiclContoursAreas further.
!=======================================================================
subroutine SemiclWilsonLoops(Egrid, nx, ny, nb, bmin, bmax, emin, emax, numE, rcell, &
                            fname, linksFile, win, isotol, mindet)

   integer,  intent(in) :: nx, ny, nb, bmin, bmax, numE
   real(dp), intent(in) :: Egrid(nx,ny,nb), emin, emax, rcell(3,3), win, isotol, mindet
   character(*), intent(in) :: fname, linksFile

   complex(dp), allocatable :: Lnk(:,:,:,:)
   real(dp), allocatable :: cx(:), cy(:), xs(:), ys(:), Ewin(:,:,:), phases(:)
   integer,  allocatable :: cstart(:), cend(:)
   logical,  allocatable :: closed(:)
   integer  :: b1, b2, nbnd, b, ie, c, nc, np, wx, wy, u, a
   integer  :: nsel, nstep, ntot, nok, nrej_sz, nrej_det
   real(dp) :: lev, dE, area, detW, isomin
   logical  :: okl, wound, ok

   call SemiclReadBerryLinks(linksFile, nx, ny, Lnk, b1, b2, okl)
   if (.not. okl) then
      call MIO_Print('  Semicl.Wilson: no usable link file, skipping','semicl')
      if (allocated(Lnk)) deallocate(Lnk)
      return
   end if
   nbnd = b2 - b1 + 1
   allocate(Ewin(nx,ny,nbnd))
   do a = 1, nbnd
      if (b1+a-1 >= 1 .and. b1+a-1 <= nb) Ewin(:,:,a) = Egrid(:,:,b1+a-1)
   end do
   allocate(cx(nx*ny), cy(nx*ny), xs(nx*ny), ys(nx*ny), phases(nbnd))
   allocate(cstart(MAXCONT), cend(MAXCONT), closed(MAXCONT))

   u = 777
   open(u, file=trim(fname), status='replace', action='write')
   write(u,'(a)') '# path-ordered (non-Abelian) Wilson loop per closed orbit (L5)'
   write(u,'(a,i0,a,i0)') '# grid ', nx, ' x ', ny
   write(u,'(a,i0,a,i0)') '# bands ', b1, ' .. ', b2
   write(u,'(a)') '# W = ordered product of the multiplet-restricted link matrices once'
   write(u,'(a)') '# around the orbit.  arg(det W) is the ABELIAN composite phase; the'
   write(u,'(a)') '# eigen-phases are the new, non-Abelian information.'
   write(u,'(a)') '# The multiplet is chosen POINTWISE (|E-lev| <= window at each node) and'
   write(u,'(a)') '# the orbit is REJECTED unless its size is constant, its isolation stays'
   write(u,'(a)') '# above the tolerance, and |det W| >= the gate.  Rejected orbits are NOT'
   write(u,'(a)') '# written: a failed gate means the number would be meaningless, not poor.'
   write(u,'(a)') '# iband  E[eV]  nsel  nstep  |detW|  isolation[eV]  phase_1 ... phase_nsel [rad]'

   dE = 0.0_dp
   if (numE > 1) dE = (emax - emin)/real(numE-1, dp)
   ntot = 0 ; nok = 0 ; nrej_sz = 0 ; nrej_det = 0
   do b = max(bmin, b1), min(bmax, b2)
      do ie = 1, numE
         lev = emin + real(ie-1,dp)*dE
         if (lev <= minval(Egrid(:,:,b)) .or. lev >= maxval(Egrid(:,:,b))) cycle
         call SemiclMarchingSquares(Egrid(:,:,b), nx, ny, lev, cx, cy, cstart, cend, nc, closed)
         do c = 1, nc
            np = cend(c) - cstart(c) + 1
            if (np < 3) cycle
            xs(1:np) = cx(cstart(c):cend(c))
            ys(1:np) = cy(cstart(c):cend(c))
            call SemiclUnwrap(xs, ys, np, wound, wx, wy)
            if (wound .or. .not. closed(c)) cycle
            call SemiclAreaCart(xs, ys, np, rcell, area)
            ntot = ntot + 1
            call SemiclOrbitWilson(Ewin, Lnk, nx, ny, nbnd, xs, ys, np, lev, win, isotol, &
                                   mindet, detW, phases, nsel, nstep, isomin, ok)
            if (.not. ok) then
               if (nsel == 0) then
                  nrej_sz = nrej_sz + 1
               else
                  nrej_det = nrej_det + 1
               end if
               cycle
            end if
            nok = nok + 1
            write(u,'(i6,2x,es13.6,2x,i4,2x,i6,2x,f10.6,2x,es13.6,16(2x,f12.8))') &
                 b, lev, nsel, nstep, detW, isomin, (phases(a), a = 1, nsel)
         end do
      end do
   end do
   close(u)
   call MIO_Print('  wrote '//trim(fname)//': '//trim(num2str(nok))//' of '// &
        trim(num2str(ntot))//' closed orbits passed every gate','semicl')
   call MIO_Print('    rejected: '//trim(num2str(nrej_sz))// &
        ' for a multiplet that changes size or loses isolation along the path, '// &
        trim(num2str(nrej_det))//' for |det W| below the gate','semicl')
   if (nok == 0) call MIO_Print('    NO orbit passed -- widen Semicl.WilsonWindow or '// &
        'lower Semicl.WilsonMinDet, but a low |det W| means the subspace did not close, '// &
        'so the phases would be meaningless','semicl')

   deallocate(Lnk, Ewin, cx, cy, xs, ys, phases, cstart, cend, closed)

end subroutine SemiclWilsonLoops

end module semicl
