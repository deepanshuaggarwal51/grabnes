subroutine DiagHam(N,ns,is,HLoc,ELoc,KLoc,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species
   use tbpar,                 only : U

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   !complex(dp), intent(out) :: H(N,N)
   !real(dp), intent(out) :: E(N)
   complex(dp), intent(out) :: HLoc(N,N)
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   real(dp) :: R(3), zz

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   !print*, cell
   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      if (ns==2) then
         !zz = charge(1,i)*charge(2,i) ! Zch
         if (is==1) then
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
         else
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
         end if
      end if
      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         !print*, "nownow", KLoc
         !print*, "niwniw", R
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KLoc,R))
      end do
   end do
   if (edgeHopp) then
      do i=1,nQ
         do j=1,nEdgeN(i)
            in = NeI(j,i)
            R = matmul(cell,NedgeCell(:,j,i))
            HLoc(in,edgeIndx(i)) = HLoc(in,edgeIndx(i)) + edgeH(j,i)*exp(cmplx_i*dot_product(KLoc,R))
         end do
      end do
   end if
   call ZHEEV('N','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
   !call ZHEEV('V','L',N,Hts,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info) 
   if (info/=0) then
      print*, "info =", info
      call MIO_Kill('Error in diagonalization','diag','DiagHam')
   end if

end subroutine DiagHam
