subroutine DiagSpectralFunctionKGridInequivalentEnergyCut()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt, frac
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is
   integer :: ik, iee, ie
   integer :: nk(3), ptot,  i1, i2, i3 !, ik, Epts, u, is, uu
   real(dp), pointer :: Kgrid(:,:)=>NULL()
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   real(dp), pointer :: pathG(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   real(dp), pointer :: KptsGFrac(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm
   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   
   logical :: MoireBS, GaussConv
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:,:)=>NULL()
   !complex(dp), pointer :: PkcLoc(:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian1(:,:)=>NULL()
   complex(dp), pointer :: Ake1Loc(:)=>NULL()
   complex(dp), pointer :: Ake2Loc(:)=>NULL()
   complex(dp), pointer :: Ake1(:,:)=>NULL()
   complex(dp), pointer :: Ake2(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian2(:,:)=>NULL()
 
   complex(dp) :: PkcLoc(nAt,2)

   real(dp) :: GVec(3) , G01(3), G10(3), G11(3), unfoldedK(3)

   real(dp) :: eps, factor, energyGridResolution
   !integer :: i1, i2

   !real(dp) :: randu
   !integer :: randi, randj
  
   real(dp) :: area, volume, grcell(3,3), aG
   real(dp) :: gcell(3,3), vn(3)
   real(dp) :: rot(3,3)

   integer :: cellSize

   real(dp) :: rcellInv(3,3)

   real(dp) :: ll, kk

   integer :: mmm(4)

   character(len=80) :: line
   integer :: id

   real(dp) :: delta, phi, gg, aa, alignmentAngle

   logical :: rotateRefSystem, alignRefSystem

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalentEnergyCut',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   !call MIO_InputParameter('Bands.MoireBS',MoireBS,.false.)
   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)
   print*, rcell(:,1)
   print*, rcell(:,2)
   print*, rcell(:,3)

   if (MIO_InputSearchLabel('MoireCellParameters',line,id)) then
       call MIO_InputParameter('MoireCellParameters',mmm,[0,0,0,0])
       !gcell(1,:) = ucell(1,:)/mmm(1)/2.0
       !gcell(2,:) = ucell(2,:)/mmm(2)/2.0
       !gcell(3,:) = ucell(3,:)
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
       call MIO_InputParameter('Spectral.RotateReferenceSystem',rotateRefSystem,.false.)
       if (rotateRefSystem .eq. .true.) then
           gg = mmm(1)**2 + mmm(2)**2 + mmm(1)*mmm(2)
           delta = sqrt(real(mmm(3)**2 + mmm(4)**2 + mmm(3)*mmm(4))/gg)
           phi = acos((2.0_dp*mmm(1)*mmm(3)+2.0_dp*mmm(2)*mmm(4) + mmm(1)*mmm(4) + mmm(2)*mmm(3))/(2.0_dp*delta*gg))
           phi = -phi*180.0_dp/pi
           aa = phi*pi/180.0_dp
           !print*, "Angle= ", aa, phi
           rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
           !rcell(:,2) = [cos(aa+pi/3.0_dp),sin(aa+pi/3.0_dp),0.0_dp]
           rot(:,2) = [sin(aa),cos(aa),0.0_dp]
           rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
           gcell = matmul(rot,gcell)
       end if
   else
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
   end if
   call MIO_InputParameter('Spectral.AlignReferenceSystem',alignRefSystem,.false.)
   if (alignRefSystem .eq. .true.) then
       call MIO_InputParameter('Spectral.AlignmentAngle',alignmentAngle,0.0d0)
       phi = alignmentAngle
       aa = -phi*pi/180.0_dp
       !print*, "Angle= ", aa, phi
       rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
       !rcell(:,2) = [cos(aa+pi/3.0_dp),sin(aa+pi/3.0_dp),0.0_dp]
       rot(:,2) = [sin(aa),cos(aa),0.0_dp]
       rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
       gcell = matmul(rot,gcell)
   end if


   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)
   area = norm(vn)
   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume

   print*, grcell(:,1)
   print*, grcell(:,2)
   print*, grcell(:,3)

   call MIO_Print('Calculating Spectral function around K','diag')
   call MIO_InputParameter('KGrid',nk,[1,1,1])
   !call MIO_InputParameter('Epsilon',eps,0.01_dp)
   !call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
   !call MIO_InputParameter('DOS.Emin',E1,-10.0_dp)
   !call MIO_InputParameter('DOS.Emax',E2,10.0_dp)
   !call MIO_Allocate(DOS,[Epts,nspin],'DOS','diag')
   !call MIO_Allocate(E,Epts,'E','diag')
   ptot = nk(1)*nk(2)*nk(3)
   call MIO_Allocate(Kgrid,[3,ptot],'Kgrid','diag')
   ik = 0
   do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      ik = ik+1
      Kgrid(:,ik) = grcell(:,1)*(2*i1-nk(1)-1)/(2.0_dp*nk(1)) + &
        grcell(:,2)*(2*i2-nk(2)-1)/(2.0_dp*nk(2)) + grcell(:,3)*(2*i2-nk(3)-1)/(2.0_dp*nk(3))
   end do; end do; end do

   !if (MIO_InputFindBlock('Spectral.Path',nPath)) then
      !call MIO_Print('Spectral function calculation','diag')
      !call MIO_Print('Based on PRB 95, 085420 (2017)','diag')
      !call MIO_Allocate(path,[3,nPath],'path','diag')
      !call MIO_Allocate(pathG,[3,nPath],'path','diag')
      !call MIO_InputBlock('Spectral.Path',path)
      !call MIO_InputBlock('Spectral.Path',pathG)
      !do ip=1,nPath
      !   print*, "0: ", path(:,ip)
      !   path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
      !   pathG(:,ip) = pathG(1,ip)*grcell(:,1) + pathG(2,ip)*grcell(:,2) + pathG(3,ip)*grcell(:,3)
      !   print*, "1: ", path(:,ip)
      !   !if (MoireBS) then
      !   !   !print*, "theta=", theta
      !   !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
      !   !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
      !   !end if
      !   print*, "2: ", path(:,ip)
      !end do
      !if (nPath==1) then
      !   call MIO_Allocate(nPts,1,'nPts','diag')
      !   nPts(1) = 1
      !   ptsTot = 1
      !else
      !   call MIO_Allocate(nPts,nPath-1,'nPts','diag')
      !   nPts(1) = nPts0
      !   ptsTot = nPts0
      !   if (nPath > 2) then
      !      v = path(:,2) - path(:,1)
      !      d0 = sqrt(dot_product(v,v))
      !      do ip=2,nPath-1
      !         v = path(:,ip+1) - path(:,ip)
      !         d = sqrt(dot_product(v,v))
      !         nPts(ip) = nint(real(d*nPts0)/real(d0))
      !         ptsTot = ptsTot + nPts(ip)
      !      end do
      !   end if
      !end if
      !call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      !call MIO_Allocate(KptsG,[3,ptsTot],'KptsG','diag')
      !call MIO_Allocate(KptsGFrac,[3,ptsTot],'KptsGFrac','diag')
      call MIO_Allocate(Kpts,[3,ptot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptot],'KptsG','diag')
      call MIO_Allocate(KptsGFrac,[3,ptot],'KptsGFrac','diag')
      !Kpts(:,1) = path(:,1)
      !KptsG(:,1) = pathG(:,1)
      KptsG(:,1) = Kgrid(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      !KptsG(:,1) = Kpts(:,1) + GVec 
      !Kpts(:,1) = KptsG(:,1) - GVec 
      !call MIO_InputParameter('CellSize', cellSize, 1)
      !print*, "cell size =", cellSize
      print*, matmul(rcell,[1,1,0]), matmul(rcell,[1,0,0]), matmul(rcell,[0,1,0])
      G10 = matmul(rcell,[1,0,0])
      G01 = matmul(rcell,[0,1,0])
      G11 = matmul(rcell,[1,1,0])
      !Kpts(1,1) = KptsG(1,1)/cellSize
      !Kpts(2,1) = KptsG(2,1)/cellSize
      !Kpts(3,1) = KptsG(3,1) 
      !kpts(1,1) = mod(KptsG(1,1),abs(G10(1)))
      !kpts(2,1) = mod(KptsG(2,1),abs(G10(2)))
      !kpts(3,1) = KptsG(3,1) 
      ll = (KptsG(1,1)*G10(2)/G10(1) - KptsG(2,1)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
      kk = (KptsG(1,1) - ll * G01(1)) / G10(1)
      if (kk.gt.0) then
          kk = floor(kk)
      else
          kk = ceiling(kk)
      end if
      if (ll.gt.0) then
          ll = floor(ll)
      else
          ll = ceiling(ll)
      end if

      !print*, ll, kk
      kpts(1,1) = KptsG(1,1) - kk*G10(1) - ll*G01(1) 
      kpts(2,1) = KptsG(2,1) - kk*G10(2) - ll*G01(2) 
 
      !print*, "tup"
      !print*, Kpts(:,1), KptsG(:,1), GVec
      !rcellInv = matinv3(rcell)
      ip = 0 
      !do i=1,nPath-1
         !do j=1,nPts(i)
         do j=1,ptot
            ip = ip + 1
            !KptsG(:,ip) = pathG(:,i) + (j-1)*(pathG(:,i+1)-pathG(:,i))/nPts(i)
            KptsG(:,ip) = Kgrid(:,ip)
            !Kpts(:,ip) = path(:,i) + (j-1)*(path(:,i+1)-path(:,i))/nPts(i)
            !call random_number(randu)
            !randi = FLOOR(22*randu) 
            !call random_number(randu)
            !randj = FLOOR(22*randu) 
            !GVec = matmul(rcell,[randi,randj,0]) ! we only want to translate them by one reciprocal lattice vector
            !KptsG(:,ip) = Kpts(:,ip) + GVec ! Kpts is in SC, Kpts is for Graphene (PC)
            !print*, "yup"
            !print*, Kpts(:,ip), KptsG(:,ip), GVec
            !print*, KptsG(2,ip) - FLOOR(KptsG(2,ip)/G2) * G2
            !print*, mod(KptsG(2,ip),G2)
            ! 2 equations, 2 unknowns. Bring point back to SC reciprocal cell
            ! using G10 and G01. 
            ll = (KptsG(1,ip)*G10(2)/G10(1) - KptsG(2,ip)) / (G01(1)*G10(2)/G10(1) - G01(2)) 
            kk = (KptsG(1,ip) - ll * G01(1)) / G10(1)
            if (kk.gt.0) then
                kk = floor(kk)
            else
                kk = ceiling(kk)
            end if
            if (ll.gt.0) then
                ll = floor(ll)
            else
                ll = ceiling(ll)
            end if

            !print*, ll, kk
            kpts(1,ip) = KptsG(1,ip) - kk*G10(1) - ll*G01(1) 
            kpts(2,ip) = KptsG(2,ip) - kk*G10(2) - ll*G01(2) 

            !KptsGFrac(2,ip) = grcell(1,1)*(KptsG(2,ip)-KptsG(1,ip)*grcell(2,1)/grcell(1,1))/(grcell(2,2)*grcell(1,1)-grcell(1,2)*grcell(2,1))
            !KptsGFrac(1,ip) = (KptsG(1,ip)-KptsGFrac(2,ip)*grcell(1,2))/grcell(1,1)
            !KptsGFrac(3,ip) = KptsG(3,ip)/grcell(3,3)
            !!!Rat(:,i) = Rat(1,i)*ucell(:,1) + Rat(2,i)*ucell(:,2) + Rat(3,i)*ucell(:,3)
            !Kpts(:,ip) = KptsGFrac(1,ip)*rcell(:,1) + KptsGFrac(2,ip)*rcell(:,2) + KptsGFrac(3,ip)*rcell(:,3)

            !Kpts(:,ip) = KptsG(:,ip) - GVec ! this one gives the same as when starting from Kpts
            !Kpts(1,ip) = KptsG(1,ip) - FLOOR(KptsG(1,ip)/G1) * G1
            !Kpts(2,ip) = KptsG(2,ip) - FLOOR(KptsG(2,ip)/G2) * G2
            !Kpts(3,ip) = KptsG(3,ip) 
            !kpts(1,ip) = mod(kptsG(1,ip),abs(G10(1)))
            !kpts(2,ip) = mod(kptsG(2,ip),abs(G10(2)))
            !kpts(3,ip) = kptsG(3,ip) 
            !Kpts(1,ip) = matmul(rcellInv,KptsG(1,ip))
            !Kpts(2,ip) = matmul(rcellInv,KptsG(2,ip)) 
            !Kpts(3,ip) = matmul(rcellInv,KptsG(3,ip)) 
            !Kpts(1,ip) = KptsG(1,ip)/cellSize
            !Kpts(2,ip) = KptsG(2,ip)/cellSize
            !Kpts(3,ip) = KptsG(3,ip) 
            !v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            !d = d + sqrt(dot_product(v,v))
         end do
      !end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectral'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptot
      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !call MIO_Allocate(Pkc,[1,1],[ptsTot,nAt],'Pkc','diag')
      call MIO_Allocate(Pkc,[1,1,1],[ptot,nAt,2],'Pkc','diag')
      !call MIO_Allocate(Pkc,nAt,'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      !print*, "factor = ", factor
      Epts2 = CEILING(Epts/factor)
      !print*, Epts2
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptot,Epts],'AkeGaussian','diag')
      call MIO_Allocate(AkeGaussian1,[ptot,Epts],'AkeGaussian1','diag')
      call MIO_Allocate(AkeGaussian2,[ptot,Epts],'AkeGaussian2','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      AkeGaussian1 = 0.0_dp
      AkeGaussian2 = 0.0_dp
      print*, "nspin ", nspin
      is = 1
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      call MIO_InputParameter('Spectral.energyGridResolution',energyGridResolution,0.005_dp)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      !call MIO_Allocate(Ake1Loc,Epts,'Ake1Loc','diag')
      !call MIO_Allocate(Ake2Loc,Epts,'Ake2Loc','diag')
      !Ake1Loc = 0.0_dp
      !Ake2Loc = 0.0_dp
      call MIO_Allocate(Ake1,[ptot,Epts],'Ake1','diag')
      call MIO_Allocate(Ake2,[ptot,Epts],'Ake2','diag')
      Ake1 = 0.0_dp
      Ake2 = 0.0_dp
      !print*, gaussian
      !Pkc = 0.0_dp ! spectral weight PkscI(k)
      !print*, "0",  Pkc
!HERE
      !$OMP PARALLEL DO PRIVATE(iee, ie, unfoldedK, ELoc, PkcLoc, KptsLoc), REDUCTION(+:Ake1Loc, Ake2Loc), &  
      !$OMP& SHARED(KptsG, Kpts, AkeGaussian1, AkeGaussian2, AkeGaussian, nAt, nspin, is, ucell, gcell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, gaussian, Epts, Ake) 
      do ik=1,ptot ! K loop
        !Ake1(ik,1) = Ake1(ik,1) + 1
         !Ake1Loc = 0.0_dp
         !Ake2Loc = 0.0_dp
      !print*, "1",  nAt
      !print*, "2",  nspin
      !print*, "3",  is
      !print*, "4",  Pkc
      !print*, "5",  E
      !print*, "6",  Kpts
      !print*, "7",  unfoldedK
      !print*, "8",  ucell
      !print*, "9",  gcell
      !print*, "10",  H0
      !print*, "11",  maxNeigh
      !print*, "12",  hopp
      !print*, "13",  NList
      !print*, "14",  Nneigh
      !print*, "15",  neighCell
      
         !ELoc = E(:,is)
         !UnfoldedK = KptsG(:,ik)
         !PkcLoc = Pkc(ik,:,:)
         !KptsLoc = Kpts(:,ik)
         ELoc = 0.0_dp
         UnfoldedK = KptsG(:,ik)
         PkcLoc = 0.0_dp
         KptsLoc = Kpts(:,ik)
         !print*, unfoldedK(1), Kpts(1,ik)
         !do while (unfoldedK(1).gt.Kpts(1,ik))
            !i = i+1
            !print*, ik, i
            !call DiagHam(nAt,nspin,is,H(:,:,is),E(:,is),Kpts(:,ip),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !call DiagSpectralWeight(nAt,nspin,is,Pkc(ik,:),E(:,is),Kpts(:,ip),KptsG(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !call DiagSpectralWeightNishi(nAt,nspin,is,Pkc(ik,:),E(:,is),Kpts(:,ik),KptsG(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !print*, "Start this"
            !print*, frac
            call DiagSpectralWeightWeiKuInequivalent(nAt,nspin,is,PkcLoc,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !print*, "Finsish this"
            !unfoldedK(1) = unfoldedK(1) - G10(1)
            !unfoldedK(2) = unfoldedK(2) - G10(2)
            !print*, unfoldedK(1), Kpts(1,ik), G1(1)
            !print*, "1 :", ik, Pkc(ik,nAt,1)
            !print*, "2 :", ik, Pkc(ik,nAt,2)
            !do i=1,nAt
            !!   if (E(i,is)<=Efermi) hv = max(hv,E(i,is))
            !!   if (E(i,is)>Efermi) lc = min(lc,E(i,is))
            !!end do
            !hv = max(hv,maxval(E(:,is),mask=E(:,is)<=Efermi/g0)) ! mask restrict search for E smaller than Efermi
            !lc = min(lc,minval(E(:,is),mask=E(:,is)>Efermi/g0))
            !v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
            !d = d + sqrt(dot_product(v,v))
            !write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is)*g0, i=1,nAt),is=1,nspin)
            !do iee=1,Epts  ! epsilon
            !    do ie=1,nAt   ! epsilonIksc
            !       is = 1
            !       !if (ie.eq.1) then
            !           if(abs(E(ie,is) - Energy(iee)).lt.(0.005/g0)) then
            !                Ake(ik,iee) = Ake(ik,iee) + abs(Pkc(ik,ie))**2
            !           end if
            !       !end if
            !    end do
            !end do
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                   !if (ie.eq.1) then
                       !if(abs(E(ie,is) - Energy(iee)).lt.(energyGridResolution/g0)) then
                       if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                            Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLoc(ie,1))**2
                            Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLoc(ie,2))**2
                            !Ake1Loc(iee) = Ake1Loc(iee) + abs(PkcLoc(ie,1))**2
                            !Ake2Loc(iee) = Ake2Loc(iee) + abs(PkcLoc(ie,1))**2
                       end if
                   !end if
                end do
            end do
         !end do
         Ake(ik,:) = Ake1(ik,:) + Ake2(ik,:)
         !Ake(ik,:) = Ake1Loc(:) + Ake2Loc(:)
         !do ie=1,nAt   ! epsilonIksc
         !   is = 1
         !   !if(abs(E(ie,is) - Energy(iee)).lt.(0.1/g0)) then
         !   iee = floor(E(ie,is)/(E2-E1)*Epts)+Epts/2
         !   if(iee.lt.1 .or. iee.gt.Epts) then
         !       cycle
         !   else
         !       Ake(ik,iee) = Ake(ik,iee) + Pkc(ik,ie)
         !   end if
         !end do
         !do ie=1,nAt   ! epsilonIksc
         !end do
         !print*, SIZE(real(Ake(ik,:))), SIZE(gaussian)
         AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
         !AkeGaussian(ik,:) = convolve(gaussian(:),real(Ake(ik,:)),Epts)
         !do i1=1,nAt
         !   do i2=1,Epts
         !      Energy(i2) = E1 + (E2-E1)*(i2-1)/(Epts-1)
         !      !DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-Eig(i1,is))**2/(2.0_dp*eps**2))
         !      Ake(ik,i2) = Ake(ik,i2) + Pkc(ik,i1)*exp(-(Energy(i2)-E(i1,is))**2/(2.0_dp*eps**2)) 
         !      !Ake(ik,i2) = Ake(ik,i2) + Pkc(ik,i1)
         !   end do
         !end do
         !end do ! this loop has to be finished first as we are still adding Pkcs to Ake
         !do ip=1,ptsTot ! kc loop
         !   v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
         !   d = d + sqrt(dot_product(v,v))
         !   do iee=1,Epts  ! epsilon
         !       ! write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is)*g0, i=1,nAt),is=1,nspin)
         !       ! float of lenght 12 with 6 after the comma
         !       ! repeat float of lengt 14 with 6 after the commq 10 times
         !       ! go to next line
         !       ! 10 empty spaces
         !       ! repeat float of lengt 14 with 6 after the commq 10 times, as many times as needed because of ()
         !       ! Ill make it simpler, but maybe bigger file
         !       write(u,'(f12.6,f12.6,f12.6)') d, Energy(iee), REAL(Ake(ik,iee))
         !   end do
         !end do
      end do
      !$OMP END PARALLEL DO
      !ik = 0
      !!do ik=1,ptsTot ! K loop
      !do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      !   ik = ik+1
      do ik=1,ptot ! K loop
         !v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
         !d = d + sqrt(dot_product(v,v))
         !if (sum(nPts(:nPath))==ik) then
         !   nPath = nPath+1
         !   call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ik))// &
         !     '   '//trim(num2str(d,6)),'diag')
         !end if
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                print*, REAL(Ake(ik,iee))
                write(u,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake(ik,iee))
            end do
         end if
      end do
      !end do; end do; end do
             

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   !end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalentEnergyCut',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunctionKGridInequivalentEnergyCut

