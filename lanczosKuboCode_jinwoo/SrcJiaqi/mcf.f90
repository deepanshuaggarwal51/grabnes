Program stuff
   Implicit None
   !$omp parallel do
   CALL mainRoutine(  )
   !$omp end parallel
Contains
   Subroutine mainRoutine()
      complex, pointer :: Pkc(:,:,:)=>NULL()
      integer :: ik, M, N, O
      M = 100
      N = 1000
      O = 2
      allocate(Pkc(M,N,O))
      Pkc = 0.0
      !$omp parallel do
      Do ik = 1, M
            call otherRoutine(Pkc(ik,:,:),N)
      End Do
      !$omp end parallel do
      !print*, Pkc(1,:,:)
      print*, "done"
   End Subroutine mainRoutine

   Subroutine otherRoutine(PkcLoc, N)
       complex, intent(out) :: PkcLoc(N,2)
       integer :: N, i, j
       PkcLoc = 0.0
       do i =1, N
          do j = 1,N
             PkcLoc(i,1) = PkcLoc(i,1) + 1.0*j
             PkcLoc(i,2) = PkcLoc(i,2) + 1.0*j
          end do
       end do

   End Subroutine otherRoutine
End Program stuff
