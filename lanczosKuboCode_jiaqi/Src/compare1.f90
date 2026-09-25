subroutine fastNNnotsquareBulk(natoms,x,y,z,aCC,cutoff2,cutoff2bis,A1,A2,A3,maxnn)
   
   use atoms,                only : inode1, inode2, in1, in2, indxNode, nAt, Species, Rat
   use tbpar,                only : tbnn
   use cell,                 only : aG, aBN,ucell
   use atoms,                only : AtomsSetCart, frac
   use math
   implicit none
   integer, intent(in) :: natoms,maxnn
   real(dp), intent(in) :: x(natoms),y(natoms),z(natoms)
   integer :: nn(natoms,maxnn)
   real(dp), intent(in) :: aCC,cutoff2,cutoff2bis,A1(3),A2(3),A3
   integer :: num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,num10,num11,num12,num13,num14
   integer :: i,j,k,l,m,n
   real(dp) :: dx(natoms,maxnn),dy(natoms,maxnn),dz(natoms,maxnn)
   !real(dp), intent(out) :: dr(natoms,maxnn)
   
   real(dp) :: dx_t,dy_t,dz_t,d2_t, d2_t2
   real(dp) :: xmin,xmax,ymin,ymax,xcell,ycell,rho,zmin,zmax
   real(dp) :: zcell
   
   real(dp) :: xnew(9*natoms*3),ynew(9*natoms*3),znew(9*natoms*3)
   integer :: ixs(9*natoms*3),iys(9*natoms*3),vecino
   integer :: izs(9*natoms*3)
   integer :: nInteger(9*natoms*3)
   integer :: mInteger(9*natoms*3)
   
   integer :: Nx,Ny,ix,iy,iz
   integer :: Nz
   integer, allocatable :: cells(:,:,:)
   integer :: NNcount,addx,addy,addz
   
   !real(dp), parameter :: rad(3) = [1.0_dp/sqrt(3.0_dp),1.0_dp,2.0_dp/sqrt(3.0_dp)]
   real(dp), parameter :: rad(5) = [1.0_dp/sqrt(3.0_dp),1.0_dp,2.0_dp/sqrt(3.0_dp),3.0_dp/sqrt(3.0_dp),3.0_dp/sqrt(3.0_dp)]
   character(len=50) :: str
 
   integer :: ncell(3,9)
   real(dp) :: d, v(3), vtemp(3)
   real(dp) :: rmax, r0, rmaxTemp

   integer, pointer :: countN(:,:), cnt(:,:)
   integer :: np
   
   integer :: safetycounter
   logical :: prnt
   logical :: BLInPlaneInteractionRadius

   real(dp) :: distFact, dist

   integer :: sCell

   logical :: readDataFiles


   call MIO_InputParameter('ReadDataFiles',readDataFiles,.false.)
   !if (frac) call AtomsSetCart()
   !call AtomsSetFrac()
   if (readDataFiles) then
      call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')

      do i=1,tbnn
         Nradii(i,1) = rad(i)*aG*1.1_dp
      end do
      if (MIO_StringComp(str,'BoronNitride')) then
         Nradii(:,2) = Nradii(:,1)
      else
         do i=1,tbnn
            Nradii(i,2) = rad(i)*aBN*1.1_dp
         end do
      end if
      rmax = maxval(Nradii)

      maxNeigh = maxnn

      print*, "allocating NList with maxNeigh= ", maxNeigh
      call MIO_Allocate(NList,[1,inode1],[maxNeigh,inode2],'NList','neigh')
      !call MIO_Allocate(NList,[1,inode1],[maxnn,inode2],'NList','neigh')
      call MIO_Allocate(Nneigh,[inode1],[inode2],'Nneigh','neigh')
      call MIO_Allocate(NeighD,[1,1,inode1],[3,maxnn,inode2],'NeighD','neigh')
      !call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,3,inode2/),'neighCell','neigh')
      call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,maxnn,inode2/),'neighCell','neigh')
      open(1,FILE='v')
      open(2,FILE='dx')
      open(3,FILE='dy')
      !open(4,FILE='pos')
      print*, "nAt =", nAt
      do i=1,nAt
         read(1,*) Nneigh(i)
         read(1,*) (NList(j,i),j=1,Nneigh(i))
         !read(2,*) (NeighD(1,j,i),j=1,Nneigh(i))
         !read(3,*) (NeighD(2,j,i),j=1,Nneigh(i))
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
      end do
      close(1)
      close(2)
      close(3)
      close(4)

      call MIO_Allocate(NList2,[1,inode1],[maxNeigh,inode2],'NList','neigh')
      NList2 = NList

!`d + idx_1 * direction_1 + idx_2 * direction_2` is smaller than `d`, where
!`idx_1` and `idx_2` are both in `[-1,0,1] and `direction_1` and `direction_2`
!are the real space directions of the translational symmetry.`

      call MIO_InputParameter('SuperCell',sCell,1)
      !print*, "ucell", ucell
      !print*, Rat
      !print*, "jjhhh"
      !print*, Rat(:,1)
      !print*, Rat(:,NList(1,1))
      !print*, "mmmm"
      !if (sCell==1) then
         !!$OMP PARALLEL PRIVATE (v, d,ncell,i,j,ix,iy,dist)
         !do i=inode1,inode2
         do i=1,nAt
            !r0 = Nradii(tbnn,(Species(i)-1)/2 + 1) ! The cut radius. Taken as the maximum distance that the
            do j=1,Nneigh(i) ! different from ultraSmall because here we already know the neighbors
               do ix=-1,1; do iy=-1,1
                  !print*, "ixiy", ix, iy
                  if (i==NList(j,i) .and. ix==0 .and. iy==0) cycle
                  ncell(:,1) = [ix,iy,0] 
                  !print*, "matmul: " , matmul(ucell,ncell(:,1))
                  v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1)) 
                  d = norm(v)
                  !print*, "d: ", d
                  if (d < aG/sqrt(3.0_dp)*1.1_dp) then
                     neighCell(:,j,i) = ncell(:,1)
                     NeighD(:,j,i) = -v
                  end if
                  !dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                  !dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                  !if (abs(d-dist).lt.0.1_dp) then
                  !   !print*, "muak"
                  !   neighCell(:,j,i) = ncell(:,1)
                  !end if
                  
                  !v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition. 
                                                                    ! Only if it is same unit cell (ix = 0, iy = 0) or neighor cell (e.g. 0 and 1), it might be satisfied (unless the bins are very small)
                  !d = norm(v)
                  !if (d**2.0_dp .lt. cutoff2*1.1_dp) then
                     !Nneigh(i) = Nneigh(i) + 1
                     !if (Nneigh(i)>maxNeigh)  then
                     !   !print*, "Number of neighbors: ", i, Nneigh(i)
                     !   call MIO_Kill('More neighbors than expected found',&
                     !     'neigh','NeighList')
                     !end if
                     !NList(Nneigh(i),i) = j
                     !NeighD(:,Nneigh(i),i) = -v
                     !neighCell(:,j,i) = ncell(:,1)
                  !else
                  !    print*, "Are you sure this was a neighbor?", i, NList(j,i)
                  !end if
               end do; end do
            end do
         end do
         !!$OMP END PARALLEL
      !end if
   else

       call MIO_InputParameter('Neigh.LayerDistFactor',distFact,1.0_dp) ! to increase the cutoff for outerlayer neighbors

       
       dx=0.0_dp;dy=0.0_dp;dz=0.0_dp !;dr=0.0_dp
       
       PRINT*,''
       PRINT*,'Building nearest neighbor data...'
       num0 = 0
       num1 = 0
       num2 = 0
       num3 = 0
       num4 = 0
       num5 = 0
       num6 = 0
       num7 = 0
       num8 = 0
       num9 = 0
       num10= 0
       num11= 0
       num12= 0
       num13= 0
       num14= 0

       ! Divide geometry into rectangular cells of size xcell times ycell (units of Angstroms)
       ! A_celdaunitaria=3.0*dsqrt(3.0)*aCC*aCC/2.0=2.6*aCC*aCC=5.24Ang => rho=2/A_celdaunitaria
       !if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'TwistedBilayer')) then
       !   rho   = 8.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC*aCC)
       !else
       !end if
       xcell = 5.0_dp
       ycell = 5.0_dp
       zcell = A3
       !zcell = 5.0_dp
       PRINT*,'Size of binning cell', xcell, 'X', ycell
       !============================================================
       !call MIO_Allocate(Species,nAt,'Species','atoms')
       !call MIO_Allocate(NeighD,[1,1,inode1],[3,maxNeigh,inode2],'NeighD','neigh')
       !call MIO_Allocate(xnew,[1],[natoms*9],'xnew','neigh')
       !call MIO_Allocate(ynew,[1],[natoms*9],'ynew','neigh')
       !call MIO_Allocate(znew,[1],[natoms*9],'znew','neigh')
       !allocate(xnew(9*natoms),ynew(9*natoms),znew(9*natoms))
       print*, "numberOfAtoms: ", natoms, nAt
       do i=1,natoms
         xnew(i)=x(i); ynew(i)=y(i); znew(i)=z(i)
         nInteger(i) = 0
         mInteger(i) = 0
       end do
       do i=1,natoms
         n=i+natoms
         xnew(n)=x(i)+A1(1)
         ynew(n)=y(i)+A1(2)
         znew(n)=z(i)
         nInteger(n) = 1
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+2*natoms
         xnew(n)=x(i)+A1(1)+A2(1)
         ynew(n)=y(i)+A1(2)+A2(2)
         znew(n)=z(i)
         nInteger(n) = 1
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+3*natoms
         xnew(n)=x(i)+A2(1)
         ynew(n)=y(i)+A2(2)
         znew(n)=z(i)
         nInteger(n) = 0
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+4*natoms
         xnew(n)=x(i)-A1(1)+A2(1)
         ynew(n)=y(i)-A1(2)+A2(2)
         znew(n)=z(i)
         nInteger(n) = -1
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+5*natoms
         xnew(n)=x(i)-A1(1)
         ynew(n)=y(i)-A1(2)
         znew(n)=z(i)
         nInteger(n) = -1
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+6*natoms
         xnew(n)=x(i)-A1(1)-A2(1)
         ynew(n)=y(i)-A1(2)-A2(2)
         znew(n)=z(i)
         nInteger(n) = -1
         mInteger(n) = -1
       end do
       do i=1,natoms
         n=i+7*natoms
         xnew(n)=x(i)-A2(1)
         ynew(n)=y(i)-A2(2)
         znew(n)=z(i)
         nInteger(n) = 0
         mInteger(n) = -1
       end do
       do i=1,natoms
         n=i+8*natoms
         xnew(n)=x(i)+A1(1)-A2(1)
         ynew(n)=y(i)+A1(2)-A2(2)
         znew(n)=z(i)
         nInteger(n) = 1
         mInteger(n) = -1
       end do
       ! add on top for bulk
       print*, "A3 =", A3
       print*, "A1 =", A1
       do i=1,natoms
         n=i + 9*natoms
         xnew(n)=x(i); ynew(n)=y(i); znew(n)=z(i)+A3
         nInteger(n) = 0
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+natoms + 9*natoms
         xnew(n)=x(i)+A1(1)
         ynew(n)=y(i)+A1(2)
         znew(n)=z(i)+A3
         nInteger(n) = 1
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+2*natoms + 9*natoms
         xnew(n)=x(i)+A1(1)+A2(1)
         ynew(n)=y(i)+A1(2)+A2(2)
         znew(n)=z(i)+A3
         nInteger(n) = 1
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+3*natoms + 9*natoms
         xnew(n)=x(i)+A2(1)
         ynew(n)=y(i)+A2(2)
         znew(n)=z(i)+A3
         nInteger(n) = 0
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+4*natoms + 9*natoms
         xnew(n)=x(i)-A1(1)+A2(1)
         ynew(n)=y(i)-A1(2)+A2(2)
         znew(n)=z(i)+A3
         nInteger(n) = -1
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+5*natoms + 9*natoms
         xnew(n)=x(i)-A1(1)
         ynew(n)=y(i)-A1(2)
         znew(n)=z(i)+A3
         nInteger(n) = -1
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+6*natoms + 9*natoms
         xnew(n)=x(i)-A1(1)-A2(1)
         ynew(n)=y(i)-A1(2)-A2(2)
         znew(n)=z(i)+A3
         nInteger(n) = -1
         mInteger(n) = -1
       end do
       do i=1,natoms
         n=i+7*natoms + 9*natoms
         xnew(n)=x(i)-A2(1)
         ynew(n)=y(i)-A2(2)
         znew(n)=z(i)+A3
         nInteger(n) = 0
         mInteger(n) = -1
       end do
       do i=1,natoms
         n=i+8*natoms + 9*natoms
         xnew(n)=x(i)+A1(1)-A2(1)
         ynew(n)=y(i)+A1(2)-A2(2)
         znew(n)=z(i)+A3
         nInteger(n) = 1
         mInteger(n) = -1
       end do
       ! add below for bulk
       do i=1,natoms
         n=i + 9*natoms + 9*natoms
         xnew(n)=x(i); ynew(n)=y(i); znew(n)=z(i)-A3
         nInteger(n) = 0
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)+A1(1)
         ynew(n)=y(i)+A1(2)
         znew(n)=z(i)-A3
         nInteger(n) = 1
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+2*natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)+A1(1)+A2(1)
         ynew(n)=y(i)+A1(2)+A2(2)
         znew(n)=z(i)-A3
         nInteger(n) = 1
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+3*natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)+A2(1)
         ynew(n)=y(i)+A2(2)
         znew(n)=z(i)-A3
         nInteger(n) = 0
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+4*natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)-A1(1)+A2(1)
         ynew(n)=y(i)-A1(2)+A2(2)
         znew(n)=z(i)-A3
         nInteger(n) = -1
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+5*natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)-A1(1)
         ynew(n)=y(i)-A1(2)
         znew(n)=z(i)-A3
         nInteger(n) = -1
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+6*natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)-A1(1)-A2(1)
         ynew(n)=y(i)-A1(2)-A2(2)
         znew(n)=z(i)-A3
         nInteger(n) = -1
         mInteger(n) = -1
       end do
       do i=1,natoms
         n=i+7*natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)-A2(1)
         ynew(n)=y(i)-A2(2)
         znew(n)=z(i)-A3
         nInteger(n) = 0
         mInteger(n) = -1
       end do
       do i=1,natoms
         n=i+8*natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)+A1(1)-A2(1)
         ynew(n)=y(i)+A1(2)-A2(2)
         znew(n)=z(i)-A3
         nInteger(n) = 1
         mInteger(n) = -1
       end do
       !============================================================
       xmin = minval(xnew); xmax = maxval(xnew)
       ymin = minval(ynew); ymax = maxval(ynew)
       zmin = minval(znew); zmax = maxval(znew)
       
       Nx = CEILING((xmax-xmin)/xcell)
       Ny = CEILING((ymax-ymin)/ycell)
       Nz = CEILING((zmax-zmin)/zcell)
       !Nz = CEILING((zmax-zmin)/zcell)

       write(*,*) xmin,xmax,Nx
       write(*,*) ymin,ymax,Ny
       !write(*,*) zmin,zmax,zy
       rho   = 4.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC*aCC)*3.0d0
       ALLOCATE(cells(Nx, Ny, 16*CEILING(xcell*ycell*rho))) ! 2 instead of 4 if not bilayer
       write(*,*) 'Inicia llenado de cajas'
       write(*,*) 'First dimension of cells is equal to ', Nx
       write(*,*) 'Second dimension of cells is equal to ', Ny
       write(*,*) 'Third dimension of cells is equal to ', 2*CEILING(xcell*ycell*rho)
        cells = 0
       !print*, cells
       do n = 1,9*natoms*3 ! Find the index ix and iy of a binning cell and fill it with all the atoms that are inside there
       
         ix     = INT(1+(xnew(n)-xmin)/xcell)
         iy     = INT(1+(ynew(n)-ymin)/ycell)
         ixs(n) = ix
         iys(n) = iy
         izs(n) = INT(1+(znew(n)-zmin)/zcell)
       
         i = 1
         !print*, ix, iy, i
         do while (cells(ix,iy,i) .ne. 0)
           !print*, i, cells(ix,iy,i)
           i = i+1
         end do
         cells(ix,iy,i) = n
         ! HERE
         !ncell(:,1) = [ix,iy,0]
         !neighCell(:,Nneigh(i),i) = ncell(:,1)
       end do
       !print*, izs
       
       
       call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')

       do i=1,tbnn
          Nradii(i,1) = rad(i)*aG*1.1_dp
       end do
       if (MIO_StringComp(str,'BoronNitride')) then
          Nradii(:,2) = Nradii(:,1)
       else
          do i=1,tbnn
             Nradii(i,2) = rad(i)*aBN*1.1_dp
          end do
       end if
       rmax = maxval(Nradii)

       maxNeigh = maxnn
       print*, "inode1, inode2, in1, in2, maxnn", inode1, inode2, 1, natoms, in1, in2, maxnn
       call MIO_Allocate(NList,[1,inode1],[maxNeigh,inode2],'NList','neigh')
       !call MIO_Allocate(NList,[1,inode1],[maxnn,inode2],'NList','neigh')
       call MIO_Allocate(Nneigh,[inode1],[inode2],'Nneigh','neigh')
       call MIO_Allocate(NeighD,[1,1,inode1],[3,maxnn,inode2],'NeighD','neigh')
       !call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,3,inode2/),'neighCell','neigh')
       call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,maxnn,inode2/),'neighCell','neigh')
       !ALLOCATE(cells(Nx, Ny, 2*CEILING(xcell*ycell*rho)))
       !call MIO_Allocate(2*CEILING(xcell*ycell*rho),
       ! HERE
       
       write(*,*) 'Inicia busqueda'
       ! Find NNs of each atom by only searching nearby cells
       safetycounter = 0
       DO n = 1,natoms
         NNcount = 0
       
         DO addx = -2,2
         DO addy = -2,2
           ix = ixs(n) + addx
           iy = iys(n) + addy
           !iz = 0 + addz
       
           ! Periodic boundary conditions here
           IF(ix .gt. Nx) ix = ix-Nx
           IF(ix .lt. 1)  ix = ix+Nx
           IF(iy .gt. Ny) iy = iy-Ny
           IF(iy .lt. 1)  iy = iy+Ny
           !IF(iz .gt. Nz) iz = iz-Nz
           !IF(iz .lt. 1)  iz = iz+Nz
       
           i = 1
           DO WHILE(cells(ix, iy, i) .ne. 0)
             m = cells(ix, iy, i)
             IF (m .ne. n) THEN
               !write(*,*) 'hi2'
               dx_t = xnew(n) - xnew(m)
               dy_t = ynew(n) - ynew(m)
               dz_t = znew(n) - znew(m)
               ! Periodic boundary conditions here as well
               ! IF(dx_t .gt.  A1(1)-1.2_dp*aCC) dx_t = dx_t - A1(1)
               ! IF(dx_t .lt. -A1(1)+1.2_dp*aCC) dx_t = dx_t + A1(1)
               ! IF(dy_t .gt.  A2(2)-1.2_dp*aCC) dy_t = dy_t - A2(2)
               ! IF(dy_t .lt. -A2(2)+1.2_dp*aCC) dy_t = dy_t + A2(2)
               d2_t = dx_t**2.0_dp + dy_t**2.0_dp + dz_t**2.0_dp
               d2_t2 = dx_t**2.0_dp + dy_t**2.0_dp
                 !IF (n==1) print*, "distances from n=1 ", d2_t
                 !IF (n==3) print*, "distances from n=1 ", d2_t
               IF (abs(dz_t)<1.50_dp) THEN
                 !if (n .eq. 1) print*, d2_t2, cutoff2
                 IF (d2_t2 .lt. cutoff2*1.1_dp) THEN
                   NNcount = NNcount + 1
                   vecino=mod(m,natoms)
                   if (vecino==0) vecino=natoms
                   nn(n,NNcount) = vecino
                   !NList(Nneigh(i),i) = j
                   dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
                   !dr(n,NNcount)=dsqrt(d2_t)
                   NList(NNcount,n) = vecino
                   NeighD(1,NNcount,n) = -dx(n,NNcount)
                   NeighD(2,NNcount,n) = -dy(n,NNcount)
                   NeighD(3,NNcount,n) = -dz(n,NNcount)
                   Nneigh(n) = NNcount
                   !neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
                   !neighCell(2,NNcount,n) = iys(m) - iys(n)
                   neighCell(1,NNcount,n) = nInteger(m) - nInteger(n)
                   neighCell(2,NNcount,n) = mInteger(m) - mInteger(n)
                   neighCell(3,NNcount,n) = izs(m) - izs(n)
                   !neighCell(3,NNcount,n) = izs(n) - izs(m)
       !             if(n==4800) write(*,*) m,NNcount,vecino
                 ENDIf
               ELSE IF (abs(dz_t)>1.50_dp .and. abs(dz_t)<4.50_dp) THEN
                 IF (d2_t2 .lt. cutoff2bis) THEN
                   NNcount = NNcount + 1
                   vecino=mod(m,natoms)
                   if (vecino==0) vecino=natoms
                   nn(n,NNcount) = vecino
                   !NList(Nneigh(i),i) = j
                   dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
                   !dr(n,NNcount)=dsqrt(d2_t)
                   NList(NNcount,n) = vecino
                   NeighD(1,NNcount,n) = -dx(n,NNcount)
                   NeighD(2,NNcount,n) = -dy(n,NNcount)
                   NeighD(3,NNcount,n) = -dz(n,NNcount)
                   Nneigh(n) = NNcount
                   !neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
                   !neighCell(2,NNcount,n) = iys(m) - iys(n)
                   !neighCell(3,NNcount,n) = 0
                   neighCell(1,NNcount,n) = nInteger(m) - nInteger(n)
                   neighCell(2,NNcount,n) = mInteger(m) - mInteger(n)
                   neighCell(3,NNcount,n) = izs(m) - izs(n)
                   !neighCell(3,NNcount,n) = izs(n) - izs(m)
                   !neighCell(3,NNcount,n) = 
       !             if(n==4800) write(*,*) m,NNcount,vecino
                 ENDIf
               ENDIF
             ENDIF
             i = i+1
           ENDDO
       
         ENDDO
         ENDDO

         !IF(NNcount .lt. 3) THEN
         !  print*, "Are you sure you wanted to enter this safety routine?"
         !  NNcount = 0
         !  DO m = 1,natoms
         !    IF (n.ne.m) THEN
         !      dx_t = x(n) - x(m)
         !      dy_t = y(n) - y(m)
         !      dz_t = z(n) - z(m)
         !       !Periodic boundary conditions here as well
         !       IF(dx_t .gt.  A1(1)-1.2_dp*aCC) dx_t = dx_t - A1(1)
         !       IF(dx_t .lt. -A1(1)+1.2_dp*aCC) dx_t = dx_t + A1(1)
         !       IF(dy_t .gt.  A2(2)-1.2_dp*aCC) dy_t = dy_t - A2(2)
         !       IF(dy_t .lt. -A2(2)+1.2_dp*aCC) dy_t = dy_t + A2(2)
         !       d2_t = dx_t**2.0_dp + dy_t**2.0_dp + dz_t**2.0_dp
         !       IF (abs(dz_t)<0.01_dp) THEN
         !          IF (d2_t .lt. cutoff2*1.1_dp) THEN
         !               NNcount = NNcount + 1
         !               !vecino=mod(m,natoms)
         !               !if (vecino==0) vecino=natoms
         !               !nn(n,NNcount) = vecino
         !               !NList(Nneigh(i),i) = j
         !               dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
         !               !dr(n,NNcount)=dsqrt(d2_t)
         !               NList(NNcount,n) = m
         !               NeighD(1,NNcount,n) = dx(n,NNcount)
         !               NeighD(2,NNcount,n) = dy(n,NNcount)
         !               NeighD(3,NNcount,n) = dz(n,NNcount)
         !               Nneigh(n) = NNcount
         !               

         !               !ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
         !               !neighCell(:,Nneigh(i),i) = ncell(:,1)
         !               

         !               !neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
         !               !neighCell(2,NNcount,n) = iys(m) - iys(n)
         !               !neighCell(3,NNcount,n) = 0
       ! !               if(n==4800) write(*,*) m,NNcount,vecino
         !          ENDIf
         !       END IF
         !     END IF
         !  END DO
         !  safetycounter = safetycounter + 1
         !  if (safetycounter.gt.100) print*, "Carefull, you are probably going too many times through this safety loop. & 
         ! I originally implemented this because the fastNNnotsquared routine was missing neighbors for 4 atoms only..."
         !ENDIF
       
       
       ! Count the number of atoms with 1, 2, 3, 4, 5 nearest neighbors
         IF(NNcount .eq. 0) num0 = num0 + 1
         !IF(NNcount .eq. 0) print*, "NN counter = 0", n
         IF(NNcount .eq. 1) num1 = num1 + 1
         IF(NNcount .eq. 2) num2 = num2 + 1
         !IF(NNcount .eq. 2) print*, "NN counter = 2", n
         IF(NNcount .eq. 3) num3 = num3 + 1
         IF(NNcount .eq. 4) num4 = num4 + 1
         IF(NNcount .eq. 5) num5 = num5 + 1
         IF(NNcount .eq. 6) num6 = num6 + 1
         IF(NNcount .eq. 7) num7 = num7 + 1
         IF(NNcount .eq. 8) num8 = num8 + 1
         IF(NNcount .eq. 9) num9 = num9 + 1
         IF(NNcount .eq. 10) num10 = num10 + 1
         IF(NNcount .eq. 11) num11 = num11 + 1
         IF(NNcount .eq. 12) num12 = num12 + 1
         IF(NNcount .eq. 13) num13 = num13 + 1
         IF(NNcount .eq. 14) num14 = num14 + 1
       
       
       
       ENDDO

       !if (frac) call AtomsSetCart()
       !call MIO_InputParameter('SuperCell',sCell,1)
       !!if (sCell==1) then
       !   !!$OMP PARALLEL PRIVATE (v, d,ncell,i,j,ix,iy,dist)
       !   !do i=inode1,inode2
       !   do i=1,nAt
       !      !r0 = Nradii(tbnn,(Species(i)-1)/2 + 1) ! The cut radius. Taken as the maximum distance that the
       !      do j=1,Nneigh(i) ! different from ultraSmall because here we already know the neighbors
       !         do ix=-3,3; do iy=-3,3!; do iz=-1,1
       !            if (i==NList(j,i) .and. ix==0 .and. iy==0) cycle! .and. iz==0) cycle
       !            ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
       !            v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1)) 
       !            d = norm(v)
       !            dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
       !            if (abs(d-dist).lt.0.000001_dp) then
       !               !print*, "muak"
       !               !print*, i, j, NList(j,i), "iz", iz
       !               neighCell(1:2,j,i) = ncell(1:2,1)
       !            end if
       !            !v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition. 
       !                                                               ! Only if it is same unit cell (ix = 0, iy = 0) or neighor cell (e.g. 0 and 1), it might be satisfied (unless the bins are very small)
       !            !d = norm(v)
       !            !if (d**2.0_dp .lt. cutoff2*1.1_dp) then
       !            !   !Nneigh(i) = Nneigh(i) + 1
       !            !   !if (Nneigh(i)>maxNeigh)  then
       !            !   !   !print*, "Number of neighbors: ", i, Nneigh(i)
       !            !   !   call MIO_Kill('More neighbors than expected found',&
       !            !   !     'neigh','NeighList')
       !            !   !end if
       !            !   !NList(Nneigh(i),i) = j
       !            !   !NeighD(:,Nneigh(i),i) = -v
       !            !   neighCell(:,j,i) = ncell(:,1)
       !            !!else
       !            !!    print*, "Are you sure this was a neighbor?", i, NList(j,i)
       !            !end if
       !         end do; end do!; end do
       !      end do
       !   end do
       !   !!$OMP END PARALLEL
       !!end if

       write(*,*) 'hi4'
       !neighCell = cells
       DEALLOCATE(cells)
       write(*,*) 'hi6'
       !call MIO_Deallocate(xnew,'xnew','neigh')
       !call MIO_Deallocate(ynew,'ynew','neigh')
       !call MIO_Deallocate(znew,'znew','neigh')
       !deallocate(xnew,ynew,znew)
       write(*,*) 'hi5'
       
       
       PRINT*,' 0 neighbors: ',num0
       PRINT*,' 1 neighbors: ',num1
       PRINT*,' 2 neighbors: ',num2
       PRINT*,' 3 neighbors: ',num3
       PRINT*,' 4 neighbors: ',num4
       PRINT*,' 5 neighbors: ',num5
       PRINT*,' 6 neighbors: ',num6
       PRINT*,' 7 neighbors: ',num7
       PRINT*,' 8 neighbors: ',num8
       PRINT*,' 9 neighbors: ',num9
       PRINT*,' 10 neighbors: ',num10
       PRINT*,' 11 neighbors: ',num11
       PRINT*,' 12 neighbors: ',num12
       PRINT*,' 13 neighbors: ',num13
       PRINT*,' 14 neighbors: ',num14
       PRINT*,'>10 neighbors: ',natoms-(num0+num1+num2+num3+num4+num5+num6+num7+num8+num9+num10+num11+num12+num13+num14)
       
       PRINT*,''
       PRINT*,'...done'
       PRINT*,''

       call MIO_Allocate(NList2,[1,inode1],[maxNeigh,inode2],'NList','neigh')
       NList2 = NList
#ifdef MPI
       !print*, "nProc =", nProc
       !if (nProc>1) then
       !   call MIO_Allocate(countN,[nProc,numThreads],'countN','neigh')
       !   countN = 0
       !   !do i=1,natoms
       !   !$OMP PARALLEL
       !   do i=in1,in2
       !      do j=1,maxNeigh
       !         if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
       !            do k=nProc,1,-1
       !               if(NList2(j,i)>=indxNode(k)) then
       !                  countN(k,nThread+1) = countN(k,nThread+1)+1
       !                  exit
       !               end if
       !            end do
       !         end if
       !      end do
       !   end do
       !   !$OMP END PARALLEL
       !   rcvSz = sum(countN)
       !   call MIO_Allocate(rcvIndx,[2,nProc],'rcvIndx','neigh')
       !   call MIO_Allocate(rcvList,rcvSz,'rcvList','neigh')
       !   np = 1
       !   do i=1,nProc
       !      rcvIndx(1,i) = np
       !      np = np + sum(countN(i,:))
       !      rcvIndx(2,i) = np - 1
       !   end do
       !   call MIO_Allocate(cnt,[nProc,numThreads],'cnt','neigh')
       !   !do i=1,natoms
       !   !$OMP PARALLEL PRIVATE(np)
       !   do i=in1,in2
       !      do j=1,maxNeigh
       !         if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
       !            do k=nProc,1,-1
       !               if(NList2(j,i)>=indxNode(k)) then
       !                  np = rcvIndx(1,k) + sum(countN(k,1:nThread)) + cnt(k,nThread+1)
       !                  cnt(k,nThread+1) = cnt(k,nThread+1) + 1
       !                  exit
       !               end if
       !            end do
       !            rcvList(np) = NList2(j,i)
       !            NList2(j,i) = np + inode2
       !         end if
       !      end do
       !   end do
       !   !$OMP END PARALLEL
       !   call MIO_Allocate(nodeSndRcv,[nProc,nProc],'nodeSndRcv','neigh')
       !   do i=1,nProc
       !      nodeSndRcv(i,Node+1) = rcvIndx(2,i)-rcvIndx(1,i) + 1
       !   end do
       !   call MPIAllGather(MPI_IN_PLACE,0,nodeSndRcv,nProc,MPI_INTEGER)
       !   sndSz = sum(nodeSndRcv(Node+1,:))
       !   call MIO_Allocate(sndList,sndSz,'sndList','neigh')
       !   call MIO_Allocate(sndIndx,(/2,nProc/),'sndIndx','neigh')
       !   np = 1
       !   do i=1,nProc
       !      sndIndx(1,i) = np
       !      np = np + nodeSndRcv(Node+1,i)
       !      sndIndx(2,i) = np - 1
       !   end do
       !   do i=1,nProc-1
       !      j = mod(Node + i,nProc)
       !      k = mod(nProc + Node - i,nProc)
       !      call MPISendRecv(rcvList(rcvIndx(1,j+1):rcvIndx(2,j+1)),nodeSndRcv(j+1,Node+1),j,50, &
       !        sndList(sndIndx(1,k+1):sndIndx(2,k+1)),nodeSndRcv(Node+1,k+1),k,50,MPI_INTEGER)
       !      call MPISendRecv(rcvList(rcvIndx(1,k+1):rcvIndx(2,k+1)),nodeSndRcv(k+1,Node+1),k,51, &
       !        sndList(sndIndx(1,j+1):sndIndx(2,j+1)),nodeSndRcv(Node+1,j+1),j,51,MPI_INTEGER)
       !      call MPIBarrier()
       !   end do
       !   call MIO_Deallocate(cnt,'cnt','neigh')
       !   call MIO_Deallocate(countN,'countN','neigh')
       !end if
#endif /* MPI */
       call MIO_InputParameter('WriteDataFiles',prnt,.false.)
       !call AtomsSetFrac()
       if (prnt) then
          open(1,FILE='v')
          open(2,FILE='dx')
          open(3,FILE='dy')
          open(33,FILE='dz')
          open(4,FILE='pos')
          do i=1,nAt
             write(1,*) Nneigh(i)
             write(1,'(300(I8,1X))') (NList(j,i),j=1,Nneigh(i)) 
             write(2,*) (NeighD(1,j,i),j=1,Nneigh(i))
             write(3,*) (NeighD(2,j,i),j=1,Nneigh(i))
             write(33,*) (NeighD(3,j,i),j=1,Nneigh(i))
             if (Species(i)==3) then
                write(4,*) "B    ", Rat(1,i),Rat(2,i),Rat(3,i)
             else if (Species(i)==4) then
                write(4,*) "N    ", Rat(1,i),Rat(2,i),Rat(3,i)
             else
                write(4,*) "C    ", Rat(1,i),Rat(2,i),Rat(3,i)
             end if
          end do
          close(1)
          close(2)
          close(3)
          close(33)
          close(4)
       end if
   end if


end subroutine fastNNnotsquareBulk
