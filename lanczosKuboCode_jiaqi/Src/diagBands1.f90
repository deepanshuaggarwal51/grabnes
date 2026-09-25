subroutine DiagBands()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is
   real(dp), pointer :: path(:,:), Kpts(:,:), E(:,:)
   integer, pointer :: nPts(:)
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: H(:,:,:)
   !type(cl_file) :: file
   character(len=100) :: flnm

#ifdef DEBUG
   call MIO_Debug('DiagBands',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_InputParameter('Bands.NumPoints',nPts0,100)
   if (MIO_InputFindBlock('Bands.Path',nPath)) then
      call MIO_Print('Band calculation','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_InputBlock('Bands.Path',path)
      do ip=1,nPath
         path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) +
path(3,ip)*rcell(:,3)
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else
         call MIO_Allocate(nPts,nPath-1,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
         end if
      end if
      call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      Kpts(:,1) = path(:,1)
      ip = 0
      d = 0.0_dp
      do i=1,nPath-1
         do j=1,nPts(i)
            ip = ip + 1
            Kpts(:,ip) = path(:,i) +
(j-1)*(path(:,i+1)-path(:,i))/nPts(i)
            v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            d = d + sqrt(dot_product(v,v))
         end do
      end do
      call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.bands'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot
      d = 0.0_dp
      hv = -huge(0.0_dp)
      lc = huge(0.0_dp)
      do ip=1,ptsTot
         do is=1,nspin
            call
DiagHam(nAt,nspin,is,H(:,:,is),E(:,is),Kpts(:,ip),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !do i=1,nAt
            !   if (E(i,is)<=Efermi) hv = max(hv,E(i,is))
            !   if (E(i,is)>Efermi) lc = min(lc,E(i,is))
            !end do
            hv = max(hv,maxval(E(:,is),mask=E(:,is)<=Efermi/g0))
            lc = min(lc,minval(E(:,is),mask=E(:,is)>Efermi/g0))
         end do
         v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
         d = d + sqrt(dot_product(v,v))
         write(u,'(f10.6,10f12.4,/,(10x,10f12.4))') d,((E(i,is)*g0,
i=1,nAt),is=1,nspin)
      end do
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(H,'H','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
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
