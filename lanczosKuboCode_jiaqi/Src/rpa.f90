subroutine KuboInitWF(Psi)

   use constants,            only : twopi, cmplx_i, cmplx_0, cmplx_1
   use parallel,             only : nDiv, procID
   !use random,               only : RandNum, rand_t, RandSeed

   complex(dp), intent(out) :: Psi(inode1:)

   integer :: i, clock, n
   real(dp) :: r
   complex(dp) :: cnum, s
   !type(rand_t), save :: rng
   integer, pointer :: seed(:)
   logical PDOS
   integer PDOSAtomNumber

#ifdef DEBUG
   call MIO_Debug('KuboInitWF',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

   call random_seed(size = n)
   allocate(seed(n))
   call system_clock(COUNT=clock)
   seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   call random_seed(PUT = seed)

   !!$OMP PARALLEL DO 
   !call RandSeed(rng,nThread)
   !do i=1,nAt
   do i=inode1, inode2
      call random_number(r)
      !Psi(i) = exp(twopi*cmplx_i*RandNum(rng))/sqrt(real(nAt))
      Psi(i) = exp(twopi*cmplx_i*r)/sqrt(real(nAt))
   end do
   !!$OMP END PARALLEL DO

   cnum = 0.0_dp
   !$OMP PARALLEL DO REDUCTION(+:cnum)
   do i=1,nAt
      cnum = cnum + Psi(i)*conjg(Psi(i))
   end do
   !$OMP END PARALLEL DO
#ifdef MPI
   call MPIRedSum(cnum,s,1,MPI_DOUBLE_COMPLEX)
   cnum = s
#endif /* MPI */
   call MIO_Print('Initial norm of wave function: '//num2str(sqrt(abs(cnum)),12),'kubo')

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboInitWF',1)
#endif /* DEBUG */

end subroutine KuboInitWF
