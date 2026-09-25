!---------------------------------------------------------!
!             MEMORY, INPUT & OUTPUT LIBRARY              !
!                                                         !
!     Copyright (C) 2012 Rafael Martinez-Gordillo         !
!                                                         !
!   This file is distributed under the terms of the       !
!   GNU General Public License. See the file `LICENSE'    !
!   in the root directory of the present distribution,    !
!   or http://www.gnu.org/copyleft/gpl.txt                !
!                                                         !
!---------------------------------------------------------!

module mem

   use precision,            only : dp, sp
   use sys,                  only : SysIOErr, SysKill

   implicit none

   PRIVATE

   integer, public, parameter :: Isz = 4                                                            
   integer, public, parameter :: Rsz = 4                                                            
   integer, public, parameter :: Dsz = 8                                                            
   integer, public, parameter :: Csz = 8                                                            
   integer, public, parameter :: Zsz = 16                                                           
   integer, public, parameter :: Ssz = 1                                                            
   integer, public, parameter :: Lsz = 4                                                            

   integer, save :: totalMem = 0, peakMem=0
   integer :: istat

   public :: MemCount
   public :: MemPrint
   public :: MemGetName

   public :: MemAlloc
   public :: MemDealloc
 
   interface MemAlloc
      module procedure Alloc_i1_s, Alloc_i1_c
      module procedure Alloc_i2_s, Alloc_i2_c
      module procedure Alloc_i3_s, Alloc_i3_c
      module procedure Alloc_i4_s, Alloc_i4_c
      module procedure Alloc_r1_s, Alloc_r1_c
      module procedure Alloc_r2_s, Alloc_r2_c
      module procedure Alloc_r3_s, Alloc_r3_c
      module procedure Alloc_r4_s, Alloc_r4_c
      module procedure Alloc_d1_s, Alloc_d1_c
      module procedure Alloc_d2_s, Alloc_d2_c
      module procedure Alloc_d3_s, Alloc_d3_c
      module procedure Alloc_d4_s, Alloc_d4_c
      module procedure Alloc_l1_s, Alloc_l1_c
      module procedure Alloc_l2_s, Alloc_l2_c
      module procedure Alloc_l3_s, Alloc_l3_c
      module procedure Alloc_l4_s, Alloc_l4_c
      module procedure Alloc_c1_s, Alloc_c1_c
      module procedure Alloc_c2_s, Alloc_c2_c
      module procedure Alloc_c3_s, Alloc_c3_c
      module procedure Alloc_c4_s, Alloc_c4_c
      module procedure Alloc_z1_s, Alloc_z1_c
      module procedure Alloc_z2_s, Alloc_z2_c
      module procedure Alloc_z3_s, Alloc_z3_c
      module procedure Alloc_z4_s, Alloc_z4_c
      module procedure Alloc_s1_s, Alloc_s1_c
      module procedure Alloc_s2_s, Alloc_s2_c
      module procedure Alloc_s3_s, Alloc_s3_c
      module procedure Alloc_s4_s, Alloc_s4_c
      module procedure Alloc_is
      module procedure Alloc_rs, Alloc_ds
      module procedure Alloc_cs, Alloc_zs
      module procedure Alloc_ls, Alloc_ss
   end interface
   interface MemDealloc
      module procedure Dealloc_i1
      module procedure Dealloc_i2
      module procedure Dealloc_i3
      module procedure Dealloc_i4
      module procedure Dealloc_r1
      module procedure Dealloc_r2
      module procedure Dealloc_r3
      module procedure Dealloc_r4
      module procedure Dealloc_d1
      module procedure Dealloc_d2
      module procedure Dealloc_d3
      module procedure Dealloc_d4
      module procedure Dealloc_l1
      module procedure Dealloc_l2
      module procedure Dealloc_l3
      module procedure Dealloc_l4
      module procedure Dealloc_c1
      module procedure Dealloc_c2
      module procedure Dealloc_c3
      module procedure Dealloc_c4
      module procedure Dealloc_z1
      module procedure Dealloc_z2
      module procedure Dealloc_z3
      module procedure Dealloc_z4
      module procedure Dealloc_s1
      module procedure Dealloc_s2
      module procedure Dealloc_s3
      module procedure Dealloc_s4
   end interface

contains

!****** Subroutine: MemCount **************************************************
!******************************************************************************
!
!  Counts the memory in use
!
! INPUT -----------------------------------------------------------------------
!
! integer sz                      : Size in bytes of the array allocated
!                                   (if sz>0) or deallocated (if sz<0)
! character arrname               : Name of the array
! character rname                 : Name of the routine that ask for memory
!
!******************************************************************************
subroutine MemCount(sz,arrname,rname)

    implicit none

    integer, intent(in) :: sz
    character(len=*), intent(in) :: arrname, rname

    totalMem = totalMem + sz
    if (totalMem > peakMem) then
        peakMem = totalMem
    end if

end subroutine MemCount
!****** End subroutine: MemCount **********************************************
!******************************************************************************


!****** Subroutine: MemPrint **************************************************
!******************************************************************************
!
!  Print the memory usage
!
!
!******************************************************************************
subroutine MemPrint()

    use format,              only : num2bytes
    use sys,                 only : SysPrint

    implicit none

    character(len=50) :: mess

    call SysPrint('')
    mess = 'Peak memory usage: '//trim(num2bytes(peakMem))
    call SysPrint(mess,'mem')

end subroutine MemPrint
!****** End subroutine: MemPrint **********************************************
!******************************************************************************


!****** Subroutine: Alloc_i1 **********************************************
subroutine Alloc_i1_c(array,szl,szu,name,routine,save)

   implicit none

   integer, pointer :: array(:)
   integer, intent(in) :: szl(1), szu(1)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer, pointer :: tmp(:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(1*2), cpBnd(1*2), indx(1*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:1) = szl
   indx(1+1:2*1) = szu
   if (associated(array)) then
      oldSz = size(array)*4
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:1) = lbound(array)
      prevBounds(1+1:2*1) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(arrSz,arrname,rname)
      array = 0
   end if
   if (copy) then
      do i=1,1
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,1
         cpBnd(i+1) = min(prevBounds(i+1),szu(i))
      end do
      array(cpBnd(1):cpBnd(2)) = &
        tmp(cpBnd(1):cpBnd(2))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_i1_c
subroutine Alloc_i1_s(array,szu,name,routine,save)

   implicit none

   integer, pointer :: array(:)
   integer, intent(in) :: szu(1)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(1)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_i1_s
!******************************************************************************
!****** Subroutine: Dealloc_i1 ********************************************
subroutine Dealloc_i1(array,name,routine)

   implicit none

   integer, pointer :: array(:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_i1
!******************************************************************************

!****** Subroutine: Alloc_i2 **********************************************
subroutine Alloc_i2_c(array,szl,szu,name,routine,save)

   implicit none

   integer, pointer :: array(:,:)
   integer, intent(in) :: szl(2), szu(2)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer, pointer :: tmp(:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(2*2), cpBnd(2*2), indx(2*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:2) = szl
   indx(2+1:2*2) = szu
   if (associated(array)) then
      oldSz = size(array)*4
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:2) = lbound(array)
      prevBounds(2+1:2*2) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(arrSz,arrname,rname)
      array = 0
   end if
   if (copy) then
      do i=1,2
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,2
         cpBnd(i+2) = min(prevBounds(i+2),szu(i))
      end do
      array(cpBnd(1):cpBnd(3),cpBnd(2):cpBnd(4)) = &
        tmp(cpBnd(1):cpBnd(3),cpBnd(2):cpBnd(4))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_i2_c
subroutine Alloc_i2_s(array,szu,name,routine,save)

   implicit none

   integer, pointer :: array(:,:)
   integer, intent(in) :: szu(2)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(2)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_i2_s
!******************************************************************************
!****** Subroutine: Dealloc_i2 ********************************************
subroutine Dealloc_i2(array,name,routine)

   implicit none

   integer, pointer :: array(:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_i2
!******************************************************************************

!****** Subroutine: Alloc_i3 **********************************************
subroutine Alloc_i3_c(array,szl,szu,name,routine,save)

   implicit none

   integer, pointer :: array(:,:,:)
   integer, intent(in) :: szl(3), szu(3)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer, pointer :: tmp(:,:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(3*2), cpBnd(3*2), indx(3*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:3) = szl
   indx(3+1:2*3) = szu
   if (associated(array)) then
      oldSz = size(array)*4
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:3) = lbound(array)
      prevBounds(3+1:2*3) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2),szl(3):szu(3)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(arrSz,arrname,rname)
      array = 0
   end if
   if (copy) then
      do i=1,3
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,3
         cpBnd(i+3) = min(prevBounds(i+3),szu(i))
      end do
      array(cpBnd(1):cpBnd(4),cpBnd(2):cpBnd(5),cpBnd(3):cpBnd(6)) = &
        tmp(cpBnd(1):cpBnd(4),cpBnd(2):cpBnd(5),cpBnd(3):cpBnd(6))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_i3_c
subroutine Alloc_i3_s(array,szu,name,routine,save)

   implicit none

   integer, pointer :: array(:,:,:)
   integer, intent(in) :: szu(3)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(3)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_i3_s
!******************************************************************************
!****** Subroutine: Dealloc_i3 ********************************************
subroutine Dealloc_i3(array,name,routine)

   implicit none

   integer, pointer :: array(:,:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_i3
!******************************************************************************

!****** Subroutine: Alloc_i4 **********************************************
subroutine Alloc_i4_c(array,szl,szu,name,routine,save)

   implicit none

   integer, pointer :: array(:,:,:,:)
   integer, intent(in) :: szl(4), szu(4)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer, pointer :: tmp(:,:,:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(4*2), cpBnd(4*2), indx(4*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:4) = szl
   indx(4+1:2*4) = szu
   if (associated(array)) then
      oldSz = size(array)*4
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:4) = lbound(array)
      prevBounds(4+1:2*4) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2),szl(3):szu(3),szl(4):szu(4)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(arrSz,arrname,rname)
      array = 0
   end if
   if (copy) then
      do i=1,4
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,4
         cpBnd(i+4) = min(prevBounds(i+4),szu(i))
      end do
      array(cpBnd(1):cpBnd(5),cpBnd(2):cpBnd(6),cpBnd(3):cpBnd(7),cpBnd(4):cpBnd(8)) = &
        tmp(cpBnd(1):cpBnd(5),cpBnd(2):cpBnd(6),cpBnd(3):cpBnd(7),cpBnd(4):cpBnd(8))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_i4_c
subroutine Alloc_i4_s(array,szu,name,routine,save)

   implicit none

   integer, pointer :: array(:,:,:,:)
   integer, intent(in) :: szu(4)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(4)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_i4_s
!******************************************************************************
!****** Subroutine: Dealloc_i4 ********************************************
subroutine Dealloc_i4(array,name,routine)

   implicit none

   integer, pointer :: array(:,:,:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_i4
!******************************************************************************

!****** Subroutine: Alloc_r1 **********************************************
subroutine Alloc_r1_c(array,szl,szu,name,routine,save)

   implicit none

   real(sp), pointer :: array(:)
   integer, intent(in) :: szl(1), szu(1)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   real(sp), pointer :: tmp(:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(1*2), cpBnd(1*2), indx(1*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:1) = szl
   indx(1+1:2*1) = szu
   if (associated(array)) then
      oldSz = size(array)*4
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:1) = lbound(array)
      prevBounds(1+1:2*1) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(arrSz,arrname,rname)
      array = 0.0_sp
   end if
   if (copy) then
      do i=1,1
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,1
         cpBnd(i+1) = min(prevBounds(i+1),szu(i))
      end do
      array(cpBnd(1):cpBnd(2)) = &
        tmp(cpBnd(1):cpBnd(2))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_r1_c
subroutine Alloc_r1_s(array,szu,name,routine,save)

   implicit none

   real(sp), pointer :: array(:)
   integer, intent(in) :: szu(1)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(1)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_r1_s
!******************************************************************************
!****** Subroutine: Dealloc_r1 ********************************************
subroutine Dealloc_r1(array,name,routine)

   implicit none

   real(sp), pointer :: array(:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_r1
!******************************************************************************

!****** Subroutine: Alloc_r2 **********************************************
subroutine Alloc_r2_c(array,szl,szu,name,routine,save)

   implicit none

   real(sp), pointer :: array(:,:)
   integer, intent(in) :: szl(2), szu(2)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   real(sp), pointer :: tmp(:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(2*2), cpBnd(2*2), indx(2*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:2) = szl
   indx(2+1:2*2) = szu
   if (associated(array)) then
      oldSz = size(array)*4
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:2) = lbound(array)
      prevBounds(2+1:2*2) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(arrSz,arrname,rname)
      array = 0.0_sp
   end if
   if (copy) then
      do i=1,2
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,2
         cpBnd(i+2) = min(prevBounds(i+2),szu(i))
      end do
      array(cpBnd(1):cpBnd(3),cpBnd(2):cpBnd(4)) = &
        tmp(cpBnd(1):cpBnd(3),cpBnd(2):cpBnd(4))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_r2_c
subroutine Alloc_r2_s(array,szu,name,routine,save)

   implicit none

   real(sp), pointer :: array(:,:)
   integer, intent(in) :: szu(2)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(2)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_r2_s
!******************************************************************************
!****** Subroutine: Dealloc_r2 ********************************************
subroutine Dealloc_r2(array,name,routine)

   implicit none

   real(sp), pointer :: array(:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_r2
!******************************************************************************

!****** Subroutine: Alloc_r3 **********************************************
subroutine Alloc_r3_c(array,szl,szu,name,routine,save)

   implicit none

   real(sp), pointer :: array(:,:,:)
   integer, intent(in) :: szl(3), szu(3)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   real(sp), pointer :: tmp(:,:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(3*2), cpBnd(3*2), indx(3*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:3) = szl
   indx(3+1:2*3) = szu
   if (associated(array)) then
      oldSz = size(array)*4
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:3) = lbound(array)
      prevBounds(3+1:2*3) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2),szl(3):szu(3)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(arrSz,arrname,rname)
      array = 0.0_sp
   end if
   if (copy) then
      do i=1,3
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,3
         cpBnd(i+3) = min(prevBounds(i+3),szu(i))
      end do
      array(cpBnd(1):cpBnd(4),cpBnd(2):cpBnd(5),cpBnd(3):cpBnd(6)) = &
        tmp(cpBnd(1):cpBnd(4),cpBnd(2):cpBnd(5),cpBnd(3):cpBnd(6))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_r3_c
subroutine Alloc_r3_s(array,szu,name,routine,save)

   implicit none

   real(sp), pointer :: array(:,:,:)
   integer, intent(in) :: szu(3)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(3)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_r3_s
!******************************************************************************
!****** Subroutine: Dealloc_r3 ********************************************
subroutine Dealloc_r3(array,name,routine)

   implicit none

   real(sp), pointer :: array(:,:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_r3
!******************************************************************************

!****** Subroutine: Alloc_r4 **********************************************
subroutine Alloc_r4_c(array,szl,szu,name,routine,save)

   implicit none

   real(sp), pointer :: array(:,:,:,:)
   integer, intent(in) :: szl(4), szu(4)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   real(sp), pointer :: tmp(:,:,:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(4*2), cpBnd(4*2), indx(4*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:4) = szl
   indx(4+1:2*4) = szu
   if (associated(array)) then
      oldSz = size(array)*4
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:4) = lbound(array)
      prevBounds(4+1:2*4) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2),szl(3):szu(3),szl(4):szu(4)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(arrSz,arrname,rname)
      array = 0.0_sp
   end if
   if (copy) then
      do i=1,4
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,4
         cpBnd(i+4) = min(prevBounds(i+4),szu(i))
      end do
      array(cpBnd(1):cpBnd(5),cpBnd(2):cpBnd(6),cpBnd(3):cpBnd(7),cpBnd(4):cpBnd(8)) = &
        tmp(cpBnd(1):cpBnd(5),cpBnd(2):cpBnd(6),cpBnd(3):cpBnd(7),cpBnd(4):cpBnd(8))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_r4_c
subroutine Alloc_r4_s(array,szu,name,routine,save)

   implicit none

   real(sp), pointer :: array(:,:,:,:)
   integer, intent(in) :: szu(4)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(4)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_r4_s
!******************************************************************************
!****** Subroutine: Dealloc_r4 ********************************************
subroutine Dealloc_r4(array,name,routine)

   implicit none

   real(sp), pointer :: array(:,:,:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_r4
!******************************************************************************

!****** Subroutine: Alloc_d1 **********************************************
subroutine Alloc_d1_c(array,szl,szu,name,routine,save)

   implicit none

   real(dp), pointer :: array(:)
   integer, intent(in) :: szl(1), szu(1)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   real(dp), pointer :: tmp(:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(1*2), cpBnd(1*2), indx(1*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:1) = szl
   indx(1+1:2*1) = szu
   if (associated(array)) then
      oldSz = size(array)*8
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:1) = lbound(array)
      prevBounds(1+1:2*1) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(arrSz,arrname,rname)
      array = 0.0_dp
   end if
   if (copy) then
      do i=1,1
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,1
         cpBnd(i+1) = min(prevBounds(i+1),szu(i))
      end do
      array(cpBnd(1):cpBnd(2)) = &
        tmp(cpBnd(1):cpBnd(2))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_d1_c
subroutine Alloc_d1_s(array,szu,name,routine,save)

   implicit none

   real(dp), pointer :: array(:)
   integer, intent(in) :: szu(1)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(1)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_d1_s
!******************************************************************************
!****** Subroutine: Dealloc_d1 ********************************************
subroutine Dealloc_d1(array,name,routine)

   implicit none

   real(dp), pointer :: array(:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_d1
!******************************************************************************

!****** Subroutine: Alloc_d2 **********************************************
subroutine Alloc_d2_c(array,szl,szu,name,routine,save)

   implicit none

   real(dp), pointer :: array(:,:)
   integer, intent(in) :: szl(2), szu(2)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   real(dp), pointer :: tmp(:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(2*2), cpBnd(2*2), indx(2*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:2) = szl
   indx(2+1:2*2) = szu
   if (associated(array)) then
      oldSz = size(array)*8
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:2) = lbound(array)
      prevBounds(2+1:2*2) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(arrSz,arrname,rname)
      array = 0.0_dp
   end if
   if (copy) then
      do i=1,2
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,2
         cpBnd(i+2) = min(prevBounds(i+2),szu(i))
      end do
      array(cpBnd(1):cpBnd(3),cpBnd(2):cpBnd(4)) = &
        tmp(cpBnd(1):cpBnd(3),cpBnd(2):cpBnd(4))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_d2_c
subroutine Alloc_d2_s(array,szu,name,routine,save)

   implicit none

   real(dp), pointer :: array(:,:)
   integer, intent(in) :: szu(2)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(2)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_d2_s
!******************************************************************************
!****** Subroutine: Dealloc_d2 ********************************************
subroutine Dealloc_d2(array,name,routine)

   implicit none

   real(dp), pointer :: array(:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_d2
!******************************************************************************

!****** Subroutine: Alloc_d3 **********************************************
subroutine Alloc_d3_c(array,szl,szu,name,routine,save)

   implicit none

   real(dp), pointer :: array(:,:,:)
   integer, intent(in) :: szl(3), szu(3)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   real(dp), pointer :: tmp(:,:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(3*2), cpBnd(3*2), indx(3*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:3) = szl
   indx(3+1:2*3) = szu
   if (associated(array)) then
      oldSz = size(array)*8
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:3) = lbound(array)
      prevBounds(3+1:2*3) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2),szl(3):szu(3)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(arrSz,arrname,rname)
      array = 0.0_dp
   end if
   if (copy) then
      do i=1,3
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,3
         cpBnd(i+3) = min(prevBounds(i+3),szu(i))
      end do
      array(cpBnd(1):cpBnd(4),cpBnd(2):cpBnd(5),cpBnd(3):cpBnd(6)) = &
        tmp(cpBnd(1):cpBnd(4),cpBnd(2):cpBnd(5),cpBnd(3):cpBnd(6))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_d3_c
subroutine Alloc_d3_s(array,szu,name,routine,save)

   implicit none

   real(dp), pointer :: array(:,:,:)
   integer, intent(in) :: szu(3)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(3)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_d3_s
!******************************************************************************
!****** Subroutine: Dealloc_d3 ********************************************
subroutine Dealloc_d3(array,name,routine)

   implicit none

   real(dp), pointer :: array(:,:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_d3
!******************************************************************************

!****** Subroutine: Alloc_d4 **********************************************
subroutine Alloc_d4_c(array,szl,szu,name,routine,save)

   implicit none

   real(dp), pointer :: array(:,:,:,:)
   integer, intent(in) :: szl(4), szu(4)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   real(dp), pointer :: tmp(:,:,:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(4*2), cpBnd(4*2), indx(4*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:4) = szl
   indx(4+1:2*4) = szu
   if (associated(array)) then
      oldSz = size(array)*8
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:4) = lbound(array)
      prevBounds(4+1:2*4) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2),szl(3):szu(3),szl(4):szu(4)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(arrSz,arrname,rname)
      array = 0.0_dp
   end if
   if (copy) then
      do i=1,4
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,4
         cpBnd(i+4) = min(prevBounds(i+4),szu(i))
      end do
      array(cpBnd(1):cpBnd(5),cpBnd(2):cpBnd(6),cpBnd(3):cpBnd(7),cpBnd(4):cpBnd(8)) = &
        tmp(cpBnd(1):cpBnd(5),cpBnd(2):cpBnd(6),cpBnd(3):cpBnd(7),cpBnd(4):cpBnd(8))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_d4_c
subroutine Alloc_d4_s(array,szu,name,routine,save)

   implicit none

   real(dp), pointer :: array(:,:,:,:)
   integer, intent(in) :: szu(4)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(4)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_d4_s
!******************************************************************************
!****** Subroutine: Dealloc_d4 ********************************************
subroutine Dealloc_d4(array,name,routine)

   implicit none

   real(dp), pointer :: array(:,:,:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_d4
!******************************************************************************

!****** Subroutine: Alloc_l1 **********************************************
subroutine Alloc_l1_c(array,szl,szu,name,routine,save)

   implicit none

   logical, pointer :: array(:)
   integer, intent(in) :: szl(1), szu(1)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   logical, pointer :: tmp(:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(1*2), cpBnd(1*2), indx(1*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:1) = szl
   indx(1+1:2*1) = szu
   if (associated(array)) then
      oldSz = size(array)*4
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:1) = lbound(array)
      prevBounds(1+1:2*1) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(arrSz,arrname,rname)
      array = .false.
   end if
   if (copy) then
      do i=1,1
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,1
         cpBnd(i+1) = min(prevBounds(i+1),szu(i))
      end do
      array(cpBnd(1):cpBnd(2)) = &
        tmp(cpBnd(1):cpBnd(2))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_l1_c
subroutine Alloc_l1_s(array,szu,name,routine,save)

   implicit none

   logical, pointer :: array(:)
   integer, intent(in) :: szu(1)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(1)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_l1_s
!******************************************************************************
!****** Subroutine: Dealloc_l1 ********************************************
subroutine Dealloc_l1(array,name,routine)

   implicit none

   logical, pointer :: array(:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_l1
!******************************************************************************

!****** Subroutine: Alloc_l2 **********************************************
subroutine Alloc_l2_c(array,szl,szu,name,routine,save)

   implicit none

   logical, pointer :: array(:,:)
   integer, intent(in) :: szl(2), szu(2)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   logical, pointer :: tmp(:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(2*2), cpBnd(2*2), indx(2*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:2) = szl
   indx(2+1:2*2) = szu
   if (associated(array)) then
      oldSz = size(array)*4
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:2) = lbound(array)
      prevBounds(2+1:2*2) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(arrSz,arrname,rname)
      array = .false.
   end if
   if (copy) then
      do i=1,2
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,2
         cpBnd(i+2) = min(prevBounds(i+2),szu(i))
      end do
      array(cpBnd(1):cpBnd(3),cpBnd(2):cpBnd(4)) = &
        tmp(cpBnd(1):cpBnd(3),cpBnd(2):cpBnd(4))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_l2_c
subroutine Alloc_l2_s(array,szu,name,routine,save)

   implicit none

   logical, pointer :: array(:,:)
   integer, intent(in) :: szu(2)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(2)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_l2_s
!******************************************************************************
!****** Subroutine: Dealloc_l2 ********************************************
subroutine Dealloc_l2(array,name,routine)

   implicit none

   logical, pointer :: array(:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_l2
!******************************************************************************

!****** Subroutine: Alloc_l3 **********************************************
subroutine Alloc_l3_c(array,szl,szu,name,routine,save)

   implicit none

   logical, pointer :: array(:,:,:)
   integer, intent(in) :: szl(3), szu(3)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   logical, pointer :: tmp(:,:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(3*2), cpBnd(3*2), indx(3*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:3) = szl
   indx(3+1:2*3) = szu
   if (associated(array)) then
      oldSz = size(array)*4
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:3) = lbound(array)
      prevBounds(3+1:2*3) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2),szl(3):szu(3)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(arrSz,arrname,rname)
      array = .false.
   end if
   if (copy) then
      do i=1,3
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,3
         cpBnd(i+3) = min(prevBounds(i+3),szu(i))
      end do
      array(cpBnd(1):cpBnd(4),cpBnd(2):cpBnd(5),cpBnd(3):cpBnd(6)) = &
        tmp(cpBnd(1):cpBnd(4),cpBnd(2):cpBnd(5),cpBnd(3):cpBnd(6))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_l3_c
subroutine Alloc_l3_s(array,szu,name,routine,save)

   implicit none

   logical, pointer :: array(:,:,:)
   integer, intent(in) :: szu(3)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(3)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_l3_s
!******************************************************************************
!****** Subroutine: Dealloc_l3 ********************************************
subroutine Dealloc_l3(array,name,routine)

   implicit none

   logical, pointer :: array(:,:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_l3
!******************************************************************************

!****** Subroutine: Alloc_l4 **********************************************
subroutine Alloc_l4_c(array,szl,szu,name,routine,save)

   implicit none

   logical, pointer :: array(:,:,:,:)
   integer, intent(in) :: szl(4), szu(4)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   logical, pointer :: tmp(:,:,:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(4*2), cpBnd(4*2), indx(4*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:4) = szl
   indx(4+1:2*4) = szu
   if (associated(array)) then
      oldSz = size(array)*4
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:4) = lbound(array)
      prevBounds(4+1:2*4) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2),szl(3):szu(3),szl(4):szu(4)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(arrSz,arrname,rname)
      array = .false.
   end if
   if (copy) then
      do i=1,4
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,4
         cpBnd(i+4) = min(prevBounds(i+4),szu(i))
      end do
      array(cpBnd(1):cpBnd(5),cpBnd(2):cpBnd(6),cpBnd(3):cpBnd(7),cpBnd(4):cpBnd(8)) = &
        tmp(cpBnd(1):cpBnd(5),cpBnd(2):cpBnd(6),cpBnd(3):cpBnd(7),cpBnd(4):cpBnd(8))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_l4_c
subroutine Alloc_l4_s(array,szu,name,routine,save)

   implicit none

   logical, pointer :: array(:,:,:,:)
   integer, intent(in) :: szu(4)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(4)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_l4_s
!******************************************************************************
!****** Subroutine: Dealloc_l4 ********************************************
subroutine Dealloc_l4(array,name,routine)

   implicit none

   logical, pointer :: array(:,:,:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*4
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_l4
!******************************************************************************

!****** Subroutine: Alloc_c1 **********************************************
subroutine Alloc_c1_c(array,szl,szu,name,routine,save)

   implicit none

   complex(sp), pointer :: array(:)
   integer, intent(in) :: szl(1), szu(1)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   complex(sp), pointer :: tmp(:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(1*2), cpBnd(1*2), indx(1*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:1) = szl
   indx(1+1:2*1) = szu
   if (associated(array)) then
      oldSz = size(array)*8
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:1) = lbound(array)
      prevBounds(1+1:2*1) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(arrSz,arrname,rname)
      array = (0.0_sp,0.0_sp)
   end if
   if (copy) then
      do i=1,1
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,1
         cpBnd(i+1) = min(prevBounds(i+1),szu(i))
      end do
      array(cpBnd(1):cpBnd(2)) = &
        tmp(cpBnd(1):cpBnd(2))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_c1_c
subroutine Alloc_c1_s(array,szu,name,routine,save)

   implicit none

   complex(sp), pointer :: array(:)
   integer, intent(in) :: szu(1)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(1)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_c1_s
!******************************************************************************
!****** Subroutine: Dealloc_c1 ********************************************
subroutine Dealloc_c1(array,name,routine)

   implicit none

   complex(sp), pointer :: array(:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_c1
!******************************************************************************

!****** Subroutine: Alloc_c2 **********************************************
subroutine Alloc_c2_c(array,szl,szu,name,routine,save)

   implicit none

   complex(sp), pointer :: array(:,:)
   integer, intent(in) :: szl(2), szu(2)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   complex(sp), pointer :: tmp(:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(2*2), cpBnd(2*2), indx(2*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:2) = szl
   indx(2+1:2*2) = szu
   if (associated(array)) then
      oldSz = size(array)*8
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:2) = lbound(array)
      prevBounds(2+1:2*2) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(arrSz,arrname,rname)
      array = (0.0_sp,0.0_sp)
   end if
   if (copy) then
      do i=1,2
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,2
         cpBnd(i+2) = min(prevBounds(i+2),szu(i))
      end do
      array(cpBnd(1):cpBnd(3),cpBnd(2):cpBnd(4)) = &
        tmp(cpBnd(1):cpBnd(3),cpBnd(2):cpBnd(4))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_c2_c
subroutine Alloc_c2_s(array,szu,name,routine,save)

   implicit none

   complex(sp), pointer :: array(:,:)
   integer, intent(in) :: szu(2)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(2)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_c2_s
!******************************************************************************
!****** Subroutine: Dealloc_c2 ********************************************
subroutine Dealloc_c2(array,name,routine)

   implicit none

   complex(sp), pointer :: array(:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_c2
!******************************************************************************

!****** Subroutine: Alloc_c3 **********************************************
subroutine Alloc_c3_c(array,szl,szu,name,routine,save)

   implicit none

   complex(sp), pointer :: array(:,:,:)
   integer, intent(in) :: szl(3), szu(3)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   complex(sp), pointer :: tmp(:,:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(3*2), cpBnd(3*2), indx(3*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:3) = szl
   indx(3+1:2*3) = szu
   if (associated(array)) then
      oldSz = size(array)*8
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:3) = lbound(array)
      prevBounds(3+1:2*3) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2),szl(3):szu(3)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(arrSz,arrname,rname)
      array = (0.0_sp,0.0_sp)
   end if
   if (copy) then
      do i=1,3
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,3
         cpBnd(i+3) = min(prevBounds(i+3),szu(i))
      end do
      array(cpBnd(1):cpBnd(4),cpBnd(2):cpBnd(5),cpBnd(3):cpBnd(6)) = &
        tmp(cpBnd(1):cpBnd(4),cpBnd(2):cpBnd(5),cpBnd(3):cpBnd(6))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_c3_c
subroutine Alloc_c3_s(array,szu,name,routine,save)

   implicit none

   complex(sp), pointer :: array(:,:,:)
   integer, intent(in) :: szu(3)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(3)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_c3_s
!******************************************************************************
!****** Subroutine: Dealloc_c3 ********************************************
subroutine Dealloc_c3(array,name,routine)

   implicit none

   complex(sp), pointer :: array(:,:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_c3
!******************************************************************************

!****** Subroutine: Alloc_c4 **********************************************
subroutine Alloc_c4_c(array,szl,szu,name,routine,save)

   implicit none

   complex(sp), pointer :: array(:,:,:,:)
   integer, intent(in) :: szl(4), szu(4)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   complex(sp), pointer :: tmp(:,:,:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(4*2), cpBnd(4*2), indx(4*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:4) = szl
   indx(4+1:2*4) = szu
   if (associated(array)) then
      oldSz = size(array)*8
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:4) = lbound(array)
      prevBounds(4+1:2*4) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2),szl(3):szu(3),szl(4):szu(4)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(arrSz,arrname,rname)
      array = (0.0_sp,0.0_sp)
   end if
   if (copy) then
      do i=1,4
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,4
         cpBnd(i+4) = min(prevBounds(i+4),szu(i))
      end do
      array(cpBnd(1):cpBnd(5),cpBnd(2):cpBnd(6),cpBnd(3):cpBnd(7),cpBnd(4):cpBnd(8)) = &
        tmp(cpBnd(1):cpBnd(5),cpBnd(2):cpBnd(6),cpBnd(3):cpBnd(7),cpBnd(4):cpBnd(8))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_c4_c
subroutine Alloc_c4_s(array,szu,name,routine,save)

   implicit none

   complex(sp), pointer :: array(:,:,:,:)
   integer, intent(in) :: szu(4)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(4)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_c4_s
!******************************************************************************
!****** Subroutine: Dealloc_c4 ********************************************
subroutine Dealloc_c4(array,name,routine)

   implicit none

   complex(sp), pointer :: array(:,:,:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*8
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_c4
!******************************************************************************

!****** Subroutine: Alloc_z1 **********************************************
subroutine Alloc_z1_c(array,szl,szu,name,routine,save)

   implicit none

   complex(dp), pointer :: array(:)
   integer, intent(in) :: szl(1), szu(1)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   complex(dp), pointer :: tmp(:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(1*2), cpBnd(1*2), indx(1*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:1) = szl
   indx(1+1:2*1) = szu
   if (associated(array)) then
      oldSz = size(array)*16
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:1) = lbound(array)
      prevBounds(1+1:2*1) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*16
      call MemCount(arrSz,arrname,rname)
      array = (0.0_dp,0.0_dp)
   end if
   if (copy) then
      do i=1,1
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,1
         cpBnd(i+1) = min(prevBounds(i+1),szu(i))
      end do
      array(cpBnd(1):cpBnd(2)) = &
        tmp(cpBnd(1):cpBnd(2))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_z1_c
subroutine Alloc_z1_s(array,szu,name,routine,save)

   implicit none

   complex(dp), pointer :: array(:)
   integer, intent(in) :: szu(1)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(1)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_z1_s
!******************************************************************************
!****** Subroutine: Dealloc_z1 ********************************************
subroutine Dealloc_z1(array,name,routine)

   implicit none

   complex(dp), pointer :: array(:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*16
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_z1
!******************************************************************************

!****** Subroutine: Alloc_z2 **********************************************
subroutine Alloc_z2_c(array,szl,szu,name,routine,save)

   implicit none

   complex(dp), pointer :: array(:,:)
   integer, intent(in) :: szl(2), szu(2)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   complex(dp), pointer :: tmp(:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(2*2), cpBnd(2*2), indx(2*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:2) = szl
   indx(2+1:2*2) = szu
   if (associated(array)) then
      oldSz = size(array)*16
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:2) = lbound(array)
      prevBounds(2+1:2*2) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*16
      call MemCount(arrSz,arrname,rname)
      array = (0.0_dp,0.0_dp)
   end if
   if (copy) then
      do i=1,2
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,2
         cpBnd(i+2) = min(prevBounds(i+2),szu(i))
      end do
      array(cpBnd(1):cpBnd(3),cpBnd(2):cpBnd(4)) = &
        tmp(cpBnd(1):cpBnd(3),cpBnd(2):cpBnd(4))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_z2_c
subroutine Alloc_z2_s(array,szu,name,routine,save)

   implicit none

   complex(dp), pointer :: array(:,:)
   integer, intent(in) :: szu(2)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(2)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_z2_s
!******************************************************************************
!****** Subroutine: Dealloc_z2 ********************************************
subroutine Dealloc_z2(array,name,routine)

   implicit none

   complex(dp), pointer :: array(:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*16
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_z2
!******************************************************************************

!****** Subroutine: Alloc_z3 **********************************************
subroutine Alloc_z3_c(array,szl,szu,name,routine,save)

   implicit none

   complex(dp), pointer :: array(:,:,:)
   integer, intent(in) :: szl(3), szu(3)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   complex(dp), pointer :: tmp(:,:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(3*2), cpBnd(3*2), indx(3*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:3) = szl
   indx(3+1:2*3) = szu
   if (associated(array)) then
      oldSz = size(array)*16
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:3) = lbound(array)
      prevBounds(3+1:2*3) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2),szl(3):szu(3)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*16
      call MemCount(arrSz,arrname,rname)
      array = (0.0_dp,0.0_dp)
   end if
   if (copy) then
      do i=1,3
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,3
         cpBnd(i+3) = min(prevBounds(i+3),szu(i))
      end do
      array(cpBnd(1):cpBnd(4),cpBnd(2):cpBnd(5),cpBnd(3):cpBnd(6)) = &
        tmp(cpBnd(1):cpBnd(4),cpBnd(2):cpBnd(5),cpBnd(3):cpBnd(6))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_z3_c
subroutine Alloc_z3_s(array,szu,name,routine,save)

   implicit none

   complex(dp), pointer :: array(:,:,:)
   integer, intent(in) :: szu(3)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(3)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_z3_s
!******************************************************************************
!****** Subroutine: Dealloc_z3 ********************************************
subroutine Dealloc_z3(array,name,routine)

   implicit none

   complex(dp), pointer :: array(:,:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*16
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_z3
!******************************************************************************

!****** Subroutine: Alloc_z4 **********************************************
subroutine Alloc_z4_c(array,szl,szu,name,routine,save)

   implicit none

   complex(dp), pointer :: array(:,:,:,:)
   integer, intent(in) :: szl(4), szu(4)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   complex(dp), pointer :: tmp(:,:,:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(4*2), cpBnd(4*2), indx(4*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:4) = szl
   indx(4+1:2*4) = szu
   if (associated(array)) then
      oldSz = size(array)*16
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:4) = lbound(array)
      prevBounds(4+1:2*4) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2),szl(3):szu(3),szl(4):szu(4)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*16
      call MemCount(arrSz,arrname,rname)
      array = (0.0_dp,0.0_dp)
   end if
   if (copy) then
      do i=1,4
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,4
         cpBnd(i+4) = min(prevBounds(i+4),szu(i))
      end do
      array(cpBnd(1):cpBnd(5),cpBnd(2):cpBnd(6),cpBnd(3):cpBnd(7),cpBnd(4):cpBnd(8)) = &
        tmp(cpBnd(1):cpBnd(5),cpBnd(2):cpBnd(6),cpBnd(3):cpBnd(7),cpBnd(4):cpBnd(8))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_z4_c
subroutine Alloc_z4_s(array,szu,name,routine,save)

   implicit none

   complex(dp), pointer :: array(:,:,:,:)
   integer, intent(in) :: szu(4)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(4)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_z4_s
!******************************************************************************
!****** Subroutine: Dealloc_z4 ********************************************
subroutine Dealloc_z4(array,name,routine)

   implicit none

   complex(dp), pointer :: array(:,:,:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*16
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_z4
!******************************************************************************

!****** Subroutine: Alloc_s1 **********************************************
subroutine Alloc_s1_c(array,szl,szu,name,routine,save)

   implicit none

   character(*), pointer :: array(:)
   integer, intent(in) :: szl(1), szu(1)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   character(len=len(array)), pointer :: tmp(:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(1*2), cpBnd(1*2), indx(1*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:1) = szl
   indx(1+1:2*1) = szu
   if (associated(array)) then
      oldSz = size(array)*len(array)*1
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:1) = lbound(array)
      prevBounds(1+1:2*1) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*len(array)*1
      call MemCount(arrSz,arrname,rname)
      array = ''
   end if
   if (copy) then
      do i=1,1
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,1
         cpBnd(i+1) = min(prevBounds(i+1),szu(i))
      end do
      array(cpBnd(1):cpBnd(2)) = &
        tmp(cpBnd(1):cpBnd(2))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_s1_c
subroutine Alloc_s1_s(array,szu,name,routine,save)

   implicit none

   character(*), pointer :: array(:)
   integer, intent(in) :: szu(1)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(1)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_s1_s
!******************************************************************************
!****** Subroutine: Dealloc_s1 ********************************************
subroutine Dealloc_s1(array,name,routine)

   implicit none

   character(*), pointer :: array(:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*len(array)*1
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_s1
!******************************************************************************

!****** Subroutine: Alloc_s2 **********************************************
subroutine Alloc_s2_c(array,szl,szu,name,routine,save)

   implicit none

   character(*), pointer :: array(:,:)
   integer, intent(in) :: szl(2), szu(2)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   character(len=len(array)), pointer :: tmp(:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(2*2), cpBnd(2*2), indx(2*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:2) = szl
   indx(2+1:2*2) = szu
   if (associated(array)) then
      oldSz = size(array)*len(array)*1
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:2) = lbound(array)
      prevBounds(2+1:2*2) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*len(array)*1
      call MemCount(arrSz,arrname,rname)
      array = ''
   end if
   if (copy) then
      do i=1,2
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,2
         cpBnd(i+2) = min(prevBounds(i+2),szu(i))
      end do
      array(cpBnd(1):cpBnd(3),cpBnd(2):cpBnd(4)) = &
        tmp(cpBnd(1):cpBnd(3),cpBnd(2):cpBnd(4))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_s2_c
subroutine Alloc_s2_s(array,szu,name,routine,save)

   implicit none

   character(*), pointer :: array(:,:)
   integer, intent(in) :: szu(2)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(2)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_s2_s
!******************************************************************************
!****** Subroutine: Dealloc_s2 ********************************************
subroutine Dealloc_s2(array,name,routine)

   implicit none

   character(*), pointer :: array(:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*len(array)*1
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_s2
!******************************************************************************

!****** Subroutine: Alloc_s3 **********************************************
subroutine Alloc_s3_c(array,szl,szu,name,routine,save)

   implicit none

   character(*), pointer :: array(:,:,:)
   integer, intent(in) :: szl(3), szu(3)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   character(len=len(array)), pointer :: tmp(:,:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(3*2), cpBnd(3*2), indx(3*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:3) = szl
   indx(3+1:2*3) = szu
   if (associated(array)) then
      oldSz = size(array)*len(array)*1
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:3) = lbound(array)
      prevBounds(3+1:2*3) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2),szl(3):szu(3)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*len(array)*1
      call MemCount(arrSz,arrname,rname)
      array = ''
   end if
   if (copy) then
      do i=1,3
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,3
         cpBnd(i+3) = min(prevBounds(i+3),szu(i))
      end do
      array(cpBnd(1):cpBnd(4),cpBnd(2):cpBnd(5),cpBnd(3):cpBnd(6)) = &
        tmp(cpBnd(1):cpBnd(4),cpBnd(2):cpBnd(5),cpBnd(3):cpBnd(6))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_s3_c
subroutine Alloc_s3_s(array,szu,name,routine,save)

   implicit none

   character(*), pointer :: array(:,:,:)
   integer, intent(in) :: szu(3)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(3)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_s3_s
!******************************************************************************
!****** Subroutine: Dealloc_s3 ********************************************
subroutine Dealloc_s3(array,name,routine)

   implicit none

   character(*), pointer :: array(:,:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*len(array)*1
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_s3
!******************************************************************************

!****** Subroutine: Alloc_s4 **********************************************
subroutine Alloc_s4_c(array,szl,szu,name,routine,save)

   implicit none

   character(*), pointer :: array(:,:,:,:)
   integer, intent(in) :: szl(4), szu(4)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   character(len=len(array)), pointer :: tmp(:,:,:,:)
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds(4*2), cpBnd(4*2), indx(4*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:4) = szl
   indx(4+1:2*4) = szu
   if (associated(array)) then
      oldSz = size(array)*len(array)*1
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:4) = lbound(array)
      prevBounds(4+1:2*4) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array(szl(1):szu(1),szl(2):szu(2),szl(3):szu(3),szl(4):szu(4)),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*len(array)*1
      call MemCount(arrSz,arrname,rname)
      array = ''
   end if
   if (copy) then
      do i=1,4
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,4
         cpBnd(i+4) = min(prevBounds(i+4),szu(i))
      end do
      array(cpBnd(1):cpBnd(5),cpBnd(2):cpBnd(6),cpBnd(3):cpBnd(7),cpBnd(4):cpBnd(8)) = &
        tmp(cpBnd(1):cpBnd(5),cpBnd(2):cpBnd(6),cpBnd(3):cpBnd(7),cpBnd(4):cpBnd(8))
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc_s4_c
subroutine Alloc_s4_s(array,szu,name,routine,save)

   implicit none

   character(*), pointer :: array(:,:,:,:)
   integer, intent(in) :: szu(4)
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl(4)

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_s4_s
!******************************************************************************
!****** Subroutine: Dealloc_s4 ********************************************
subroutine Dealloc_s4(array,name,routine)

   implicit none

   character(*), pointer :: array(:,:,:,:)
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*len(array)*1
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc_s4
!******************************************************************************
!****** Subroutine: Alloc_is***************************************************
subroutine Alloc_is(array,sz,name,routine,save)

   implicit none

   integer, pointer, intent(inout) :: array(:)
   integer, intent(in) :: sz
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szu(1), szl(1)

   szu = sz
   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_is
!******************************************************************************
!****** Subroutine: Alloc_rs **************************************************
subroutine Alloc_rs(array,sz,name,routine,save)

   implicit none

   real(sp), pointer, intent(inout) :: array(:)
   integer, intent(in) :: sz
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szu(1), szl(1)

   szu = sz
   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_rs
!******************************************************************************
!****** Subroutine: Alloc_ds **************************************************
subroutine Alloc_ds(array,sz,name,routine,save)

   implicit none

   real(dp), pointer, intent(inout) :: array(:)
   integer, intent(in) :: sz
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szu(1), szl(1)

   szu = sz
   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_ds
!******************************************************************************
!****** Subroutine: Alloc_cs **************************************************
subroutine Alloc_cs(array,sz,name,routine,save)

   implicit none

   complex(sp), pointer, intent(inout) :: array(:)
   integer, intent(in) :: sz
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szu(1), szl(1)

   szu = sz
   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_cs
!******************************************************************************
!****** Subroutine: Alloc_zs **************************************************
subroutine Alloc_zs(array,sz,name,routine,save)

   implicit none

   complex(dp), pointer, intent(inout) :: array(:)
   integer, intent(in) :: sz
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szu(1), szl(1)

   szu = sz
   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_zs
!******************************************************************************
!****** Subroutine: Alloc_ls **************************************************
subroutine Alloc_ls(array,sz,name,routine,save)

   implicit none

   logical, pointer, intent(inout) :: array(:)
   integer, intent(in) :: sz
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szu(1), szl(1)

   szu = sz
   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_ls
!******************************************************************************
!****** Subroutine: Alloc_ss **************************************************
subroutine Alloc_ss(array,sz,name,routine,save)

   implicit none

   character(*), pointer, intent(inout) :: array(:)
   integer, intent(in) :: sz
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szu(1), szl(1)

   szu = sz
   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_ss
!******************************************************************************

!****** Subroutine: MemGetName ************************************************
!******************************************************************************
!
!  Get names for the variable and the routine
!
! INPUT -----------------------------------------------------------------------
!
! character name                 : Name of the variable
! character routine              : Calling routine
!
! OUTPUT ----------------------------------------------------------------------
!
! character arrname              : Final name of the array
! character rname                : Final name of the routine
!
!******************************************************************************
subroutine MemGetName(arrname,rname,name,routine)

   implicit none

   character(len=32), intent(out) :: arrname, rname
   character(len=*), intent(in), optional :: name, routine

   if (present(name)) then
      arrname = name
   else
      arrname = 'unknown'
   end if
   if (present(routine)) then
      rname = routine
   else
      rname = 'unknown'
   end if

end subroutine MemGetName
!****** End subroutine: MemGetName ********************************************
!******************************************************************************


end module mem
