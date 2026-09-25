subroutine DiagBands()

   use cell,                 only : rcell, ucell, aG
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: H(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm

   logical :: MoireBS
   real(dp) :: theta, volume

   real(dp) :: gcell(3,3), grcell(3,3)

   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   complex(dp) :: HLoc(nAt, nAt)

#ifdef DEBUG
   call MIO_Debug('DiagBands',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   !call MIO_InputParameter('Bands.MoireBS',MoireBS,.false.)
   call MIO_InputParameter('Bands.NumPoints',nPts0,100)
   print*, rcell(:,1)
   print*, rcell(:,2)
   print*, rcell(:,3)
   if (MIO_InputFindBlock('Bands.Path',nPath)) then
      call MIO_Print('Band calculation','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_InputBlock('Bands.Path',path)
      do ip=1,nPath
         print*, "0: ", path(:,ip)
         path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
         print*, "1: ", path(:,ip)
         !if (MoireBS) then
         !   !print*, "theta=", theta
         !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
         !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
         !end if
         print*, "2: ", path(:,ip)
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else 
         !call MIO_Allocate(nPts,nPath-1,'nPts','diag')
         call MIO_Allocate(nPts,nPath,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1 ! coz 4 points, 3 segments
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
            nPts(ip) = nPts(ip)  + 1 ! Add the missing point at the end of last segment
         end if
      end if
      !call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      call MIO_Allocate(Kpts,[3,ptsTot+1],'Kpts','diag') ! add missing point
      print*, "ptsTot, nPath", ptsTot, nPath, nPts(1)
      Kpts(:,1) = path(:,1)
      !ip = 0
      ip = 1
      d = 0.0_dp
      do i=1,nPath-1
      !do i=2,nPath
         do j=1,nPts(i)
            ip = ip + 1
            !Kpts(:,ip) = path(:,i) + (j-1)*(path(:,i+1)-path(:,i))/nPts(i)
            Kpts(:,ip) = path(:,i) + (j)*(path(:,i+1)-path(:,i))/nPts(i)
            v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            d = d + sqrt(dot_product(v,v))
         end do
      end do
      print*, "initial Gamma: ", Kpts(:,1)
      print*, "final Gamma: ", Kpts(:,ip), ip
      call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
      call MIO_Allocate(E,[nAt,nspin,ptsTot+1],'E','diag')
      flnm = trim(prefix)//'.bands'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot+1
      d = 0.0_dp
      hv = -huge(0.0_dp)
      lc = huge(0.0_dp)
      E = 0.0_dp
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !$OMP PARALLEL DO PRIVATE(ELoc, KptsLoc, is, ip, HLoc), &
      !$OMP& SHARED(E, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, Kpts)
      do ip=1,ptsTot + 1
         do is=1,nspin
            ELoc = 0.0_dp
            HLoc = 0.0_dp
            KptsLoc = Kpts(:,ip)
            call DiagHam(nAt,nspin,is,HLoc,ELoc,KptsLoc,ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            E(:,is,ip) = ELoc
         end do
      end do
      !$OMP END PARALLEL DO
      do ip=1,ptsTot + 1
         v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
         !v = Kpts(:,ip+1) - Kpts(:,max(ip,1))
         d = d + sqrt(dot_product(v,v))
         write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is,ip)*g0, i=1,nAt),is=1,nspin)
         if (sum(nPts(:nPath))==ip-1) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ip))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
      end do
      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(H,'H','diag')
      !call MIO_Deallocate(Kpts,'Htsp','diag')
      !call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagBands',1)
#endif /* DEBUG */

end subroutine DiagBands
