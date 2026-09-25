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
                                                                                                    
module mpi_i                                                                                        
                                                                                                    
   use mpitime,              only : MPITimer                                                        
   use mpi_inc                                                                                      
   use mpikinds                                                                                     
                                                                                                    
   implicit none                                                                                    
                                                                                                    
   PRIVATE                                                                                          
                                                                                                    
   integer :: ierr                                                                                  
                                                                                                    
   public :: MPIBCast                                                                               
   interface MPIBCast                                                                               
      module procedure MPIBCast_i, MPIBCast_i_v                                                     
   end interface                                                                                    
                                                                                                    
   public :: MPIRedSum                                                                              
   interface MPIRedSum                                                                              
      module procedure MPIRedSum_i,  MPIRedSum_i_v                                                  
   end interface                                                                                    
                                                                                                    
   public :: MPIAllRedSum                                                                           
   interface MPIAllRedSum                                                                           
      module procedure MPIAllRedSum_i,  MPIAllRedSum_i_v                                            
   end interface                                                                                    
                                                                                                    
   public :: MPIGather                                                                              
   interface MPIGather                                                                              
      module procedure MPIGather_i,  MPIGather_i_v                                                  
   end interface                                                                                    
                                                                                                    
   public :: MPIAllGather                                                                           
   interface MPIAllGather                                                                           
      module procedure MPIAllGather_i,  MPIAllGather_i_v                                            
      module procedure MPIAllGather2_i, MPIAllGather2_i_v                                           
   end interface                                                                                    
                                                                                                    
   public :: MPIAllGatherV                                                                          
   interface MPIAllGatherV                                                                          
      module procedure MPIAllGatherV_i, MPIAllGatherV_i_v                                           
   end interface                                                                                    
                                                                                                    
   public :: MPISendRecv                                                                            
   interface MPISendRecv                                                                            
      module procedure MPISendRecv_i                                                                
   end interface                                                                                    
                                                                                                    
contains                                                                                            
                                                                                                    
subroutine MPIBCast_i(buffer,sz,type)                                                               
   integer, intent(inout) :: buffer                                                                 
   integer :: sz                                                                                    
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)                                     
   call MPITimer(0)                                                                                 
end subroutine MPIBCast_i                                                                           
subroutine MPIBCast_i_v(buffer,sz,type)                                                             
   integer, intent(inout) :: buffer(*)                                                              
   integer :: sz                                                                                    
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)                                     
   call MPITimer(0)                                                                                 
end subroutine MPIBCast_i_v                                                                         
                                                                                                    
subroutine MPIRedSum_i(sendbuf,recvbuf,sz,type)                                                     
   integer, intent(in) :: sendbuf                                                                   
   integer, intent(out) :: recvbuf                                                                  
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)                   
   call MPITimer(0)                                                                                 
end subroutine MPIRedSum_i                                                                          
subroutine MPIRedSum_i_v(sendbuf,recvbuf,sz,type)                                                   
   integer, intent(in) :: sendbuf(*)                                                                
   integer, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)                   
   call MPITimer(0)                                                                                 
end subroutine MPIRedSum_i_v                                                                        
                                                                                                    
subroutine MPIAllRedSum_i(sendbuf,recvbuf,sz,type)                                                  
   integer, intent(in) :: sendbuf                                                                   
   integer, intent(out) :: recvbuf                                                                  
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)                          
   call MPITimer(0)                                                                                 
end subroutine MPIAllRedSum_i                                                                       
subroutine MPIAllRedSum_i_v(sendbuf,recvbuf,sz,type)                                                
   integer, intent(in) :: sendbuf(*)                                                                
   integer, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)                          
   call MPITimer(0)                                                                                 
end subroutine MPIAllRedSum_i_v                                                                     
                                                                                                    
subroutine MPIGather_i(sendbuf,sendsz,recvbuf,recvsz,type)                                          
   integer, intent(in) :: sendbuf                                                                   
   integer, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)           
   call MPITimer(0)                                                                                 
end subroutine MPIGather_i                                                                          
subroutine MPIGather_i_v(sendbuf,sendsz,recvbuf,recvsz,type)                                        
   integer, intent(in) :: sendbuf(*)                                                                
   integer, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)           
   call MPITimer(0)                                                                                 
end subroutine MPIGather_i_v                                                                        
                                                                                                    
subroutine MPIAllGather_i(sendbuf,sendsz,recvbuf,recvsz,type)                                       
   integer, intent(in) :: sendbuf                                                                   
   integer, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather_i                                                                       
subroutine MPIAllGather_i_v(sendbuf,sendsz,recvbuf,recvsz,type)                                     
   integer, intent(in) :: sendbuf(*)                                                                
   integer, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather_i_v                                                                     
subroutine MPIAllGather2_i(sendbuf,sendsz,recvbuf,recvsz,type)                                      
   integer, intent(in) :: sendbuf                                                                   
   integer, intent(out) :: recvbuf(:,:)                                                             
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather2_i                                                                      
subroutine MPIAllGather2_i_v(sendbuf,sendsz,recvbuf,recvsz,type)                                    
   integer, intent(in) :: sendbuf(*)                                                                
   integer, intent(out) :: recvbuf(:,:)                                                             
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather2_i_v                                                                    
                                                                                                    
                                                                                                    
subroutine MPIAllGatherV_i(sendbuf,sendsz,recvbuf,recvsz,displs,type)                               
   integer, intent(in) :: sendbuf                                                                   
   integer, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sendsz, recvsz(*), displs(*)                                              
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)          
   call MPITimer(0)                                                                                 
end subroutine MPIAllGatherV_i                                                                      
subroutine MPIAllGatherV_i_v(sendbuf,sendsz,recvbuf,recvsz,displs,type)                             
   integer, intent(in) :: sendbuf(*)                                                                
   integer, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sendsz, recvsz(*), displs(*)                                              
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)          
   call MPITimer(0)                                                                                 
end subroutine MPIAllGatherV_i_v                                                                    
                                                                                                    
subroutine MPISendRecv_i(sendbuf,sendsz,dest,sendtag,recvbuf,recvsz,source,recvtag,type)            
   integer, intent(in) :: sendbuf(*)                                                                
   integer, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sendsz, recvsz, dest, source                                              
   integer, intent(in) :: sendtag, recvtag, type                                                    
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_SendRecv(sendbuf,sendsz,type,dest,sendtag,recvbuf,recvsz,type,source,recvtag,&          
     MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)                                                         
   call MPITimer(0)                                                                                 
end subroutine MPISendRecv_i                                                                        
                                                                                                    
end module mpi_i                                                                                    
                                                                                                    
module mpi_r                                                                                        
                                                                                                    
   use mpitime,              only : MPITimer                                                        
   use mpi_inc                                                                                      
   use mpikinds                                                                                     
                                                                                                    
   implicit none                                                                                    
                                                                                                    
   PRIVATE                                                                                          
                                                                                                    
   integer :: ierr                                                                                  
                                                                                                    
   public :: MPIBCast                                                                               
   interface MPIBCast                                                                               
      module procedure MPIBCast_r, MPIBCast_r_v                                                     
   end interface                                                                                    
                                                                                                    
   public :: MPIRedSum                                                                              
   interface MPIRedSum                                                                              
      module procedure MPIRedSum_r,  MPIRedSum_r_v                                                  
   end interface                                                                                    
                                                                                                    
   public :: MPIAllRedSum                                                                           
   interface MPIAllRedSum                                                                           
      module procedure MPIAllRedSum_r,  MPIAllRedSum_r_v                                            
   end interface                                                                                    
                                                                                                    
   public :: MPIGather                                                                              
   interface MPIGather                                                                              
      module procedure MPIGather_r,  MPIGather_r_v                                                  
   end interface                                                                                    
                                                                                                    
   public :: MPIAllGather                                                                           
   interface MPIAllGather                                                                           
      module procedure MPIAllGather_r,  MPIAllGather_r_v                                            
      module procedure MPIAllGather2_r, MPIAllGather2_r_v                                           
   end interface                                                                                    
                                                                                                    
   public :: MPIAllGatherV                                                                          
   interface MPIAllGatherV                                                                          
      module procedure MPIAllGatherV_r, MPIAllGatherV_r_v                                           
   end interface                                                                                    
                                                                                                    
   public :: MPISendRecv                                                                            
   interface MPISendRecv                                                                            
      module procedure MPISendRecv_r                                                                
   end interface                                                                                    
                                                                                                    
contains                                                                                            
                                                                                                    
subroutine MPIBCast_r(buffer,sz,type)                                                               
   real(sp), intent(inout) :: buffer                                                                
   integer :: sz                                                                                    
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)                                     
   call MPITimer(0)                                                                                 
end subroutine MPIBCast_r                                                                           
subroutine MPIBCast_r_v(buffer,sz,type)                                                             
   real(sp), intent(inout) :: buffer(*)                                                             
   integer :: sz                                                                                    
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)                                     
   call MPITimer(0)                                                                                 
end subroutine MPIBCast_r_v                                                                         
                                                                                                    
subroutine MPIRedSum_r(sendbuf,recvbuf,sz,type)                                                     
   real(sp), intent(in) :: sendbuf                                                                  
   real(sp), intent(out) :: recvbuf                                                                 
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)                   
   call MPITimer(0)                                                                                 
end subroutine MPIRedSum_r                                                                          
subroutine MPIRedSum_r_v(sendbuf,recvbuf,sz,type)                                                   
   real(sp), intent(in) :: sendbuf(*)                                                               
   real(sp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)                   
   call MPITimer(0)                                                                                 
end subroutine MPIRedSum_r_v                                                                        
                                                                                                    
subroutine MPIAllRedSum_r(sendbuf,recvbuf,sz,type)                                                  
   real(sp), intent(in) :: sendbuf                                                                  
   real(sp), intent(out) :: recvbuf                                                                 
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)                          
   call MPITimer(0)                                                                                 
end subroutine MPIAllRedSum_r                                                                       
subroutine MPIAllRedSum_r_v(sendbuf,recvbuf,sz,type)                                                
   real(sp), intent(in) :: sendbuf(*)                                                               
   real(sp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)                          
   call MPITimer(0)                                                                                 
end subroutine MPIAllRedSum_r_v                                                                     
                                                                                                    
subroutine MPIGather_r(sendbuf,sendsz,recvbuf,recvsz,type)                                          
   real(sp), intent(in) :: sendbuf                                                                  
   real(sp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)           
   call MPITimer(0)                                                                                 
end subroutine MPIGather_r                                                                          
subroutine MPIGather_r_v(sendbuf,sendsz,recvbuf,recvsz,type)                                        
   real(sp), intent(in) :: sendbuf(*)                                                               
   real(sp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)           
   call MPITimer(0)                                                                                 
end subroutine MPIGather_r_v                                                                        
                                                                                                    
subroutine MPIAllGather_r(sendbuf,sendsz,recvbuf,recvsz,type)                                       
   real(sp), intent(in) :: sendbuf                                                                  
   real(sp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather_r                                                                       
subroutine MPIAllGather_r_v(sendbuf,sendsz,recvbuf,recvsz,type)                                     
   real(sp), intent(in) :: sendbuf(*)                                                               
   real(sp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather_r_v                                                                     
subroutine MPIAllGather2_r(sendbuf,sendsz,recvbuf,recvsz,type)                                      
   real(sp), intent(in) :: sendbuf                                                                  
   real(sp), intent(out) :: recvbuf(:,:)                                                            
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather2_r                                                                      
subroutine MPIAllGather2_r_v(sendbuf,sendsz,recvbuf,recvsz,type)                                    
   real(sp), intent(in) :: sendbuf(*)                                                               
   real(sp), intent(out) :: recvbuf(:,:)                                                            
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather2_r_v                                                                    
                                                                                                    
                                                                                                    
subroutine MPIAllGatherV_r(sendbuf,sendsz,recvbuf,recvsz,displs,type)                               
   integer, intent(in) :: sendbuf                                                                   
   real(sp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sendsz, recvsz(*), displs(*)                                              
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)          
   call MPITimer(0)                                                                                 
end subroutine MPIAllGatherV_r                                                                      
subroutine MPIAllGatherV_r_v(sendbuf,sendsz,recvbuf,recvsz,displs,type)                             
   real(sp), intent(in) :: sendbuf(*)                                                               
   real(sp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sendsz, recvsz(*), displs(*)                                              
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)          
   call MPITimer(0)                                                                                 
end subroutine MPIAllGatherV_r_v                                                                    
                                                                                                    
subroutine MPISendRecv_r(sendbuf,sendsz,dest,sendtag,recvbuf,recvsz,source,recvtag,type)            
   real(sp), intent(in) :: sendbuf(*)                                                               
   real(sp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sendsz, recvsz, dest, source                                              
   integer, intent(in) :: sendtag, recvtag, type                                                    
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_SendRecv(sendbuf,sendsz,type,dest,sendtag,recvbuf,recvsz,type,source,recvtag,&          
     MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)                                                         
   call MPITimer(0)                                                                                 
end subroutine MPISendRecv_r                                                                        
                                                                                                    
end module mpi_r                                                                                    
                                                                                                    
module mpi_d                                                                                        
                                                                                                    
   use mpitime,              only : MPITimer                                                        
   use mpi_inc                                                                                      
   use mpikinds                                                                                     
                                                                                                    
   implicit none                                                                                    
                                                                                                    
   PRIVATE                                                                                          
                                                                                                    
   integer :: ierr                                                                                  
                                                                                                    
   public :: MPIBCast                                                                               
   interface MPIBCast                                                                               
      module procedure MPIBCast_d, MPIBCast_d_v                                                     
   end interface                                                                                    
                                                                                                    
   public :: MPIRedSum                                                                              
   interface MPIRedSum                                                                              
      module procedure MPIRedSum_d,  MPIRedSum_d_v                                                  
   end interface                                                                                    
                                                                                                    
   public :: MPIAllRedSum                                                                           
   interface MPIAllRedSum                                                                           
      module procedure MPIAllRedSum_d,  MPIAllRedSum_d_v                                            
   end interface                                                                                    
                                                                                                    
   public :: MPIGather                                                                              
   interface MPIGather                                                                              
      module procedure MPIGather_d,  MPIGather_d_v                                                  
   end interface                                                                                    
                                                                                                    
   public :: MPIAllGather                                                                           
   interface MPIAllGather                                                                           
      module procedure MPIAllGather_d,  MPIAllGather_d_v                                            
      module procedure MPIAllGather2_d, MPIAllGather2_d_v                                           
   end interface                                                                                    
                                                                                                    
   public :: MPIAllGatherV                                                                          
   interface MPIAllGatherV                                                                          
      module procedure MPIAllGatherV_d, MPIAllGatherV_d_v                                           
   end interface                                                                                    
                                                                                                    
   public :: MPISendRecv                                                                            
   interface MPISendRecv                                                                            
      module procedure MPISendRecv_d                                                                
   end interface                                                                                    
                                                                                                    
contains                                                                                            
                                                                                                    
subroutine MPIBCast_d(buffer,sz,type)                                                               
   real(dp), intent(inout) :: buffer                                                                
   integer :: sz                                                                                    
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)                                     
   call MPITimer(0)                                                                                 
end subroutine MPIBCast_d                                                                           
subroutine MPIBCast_d_v(buffer,sz,type)                                                             
   real(dp), intent(inout) :: buffer(*)                                                             
   integer :: sz                                                                                    
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)                                     
   call MPITimer(0)                                                                                 
end subroutine MPIBCast_d_v                                                                         
                                                                                                    
subroutine MPIRedSum_d(sendbuf,recvbuf,sz,type)                                                     
   real(dp), intent(in) :: sendbuf                                                                  
   real(dp), intent(out) :: recvbuf                                                                 
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)                   
   call MPITimer(0)                                                                                 
end subroutine MPIRedSum_d                                                                          
subroutine MPIRedSum_d_v(sendbuf,recvbuf,sz,type)                                                   
   real(dp), intent(in) :: sendbuf(*)                                                               
   real(dp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)                   
   call MPITimer(0)                                                                                 
end subroutine MPIRedSum_d_v                                                                        
                                                                                                    
subroutine MPIAllRedSum_d(sendbuf,recvbuf,sz,type)                                                  
   real(dp), intent(in) :: sendbuf                                                                  
   real(dp), intent(out) :: recvbuf                                                                 
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)                          
   call MPITimer(0)                                                                                 
end subroutine MPIAllRedSum_d                                                                       
subroutine MPIAllRedSum_d_v(sendbuf,recvbuf,sz,type)                                                
   real(dp), intent(in) :: sendbuf(*)                                                               
   real(dp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)                          
   call MPITimer(0)                                                                                 
end subroutine MPIAllRedSum_d_v                                                                     
                                                                                                    
subroutine MPIGather_d(sendbuf,sendsz,recvbuf,recvsz,type)                                          
   real(dp), intent(in) :: sendbuf                                                                  
   real(dp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)           
   call MPITimer(0)                                                                                 
end subroutine MPIGather_d                                                                          
subroutine MPIGather_d_v(sendbuf,sendsz,recvbuf,recvsz,type)                                        
   real(dp), intent(in) :: sendbuf(*)                                                               
   real(dp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)           
   call MPITimer(0)                                                                                 
end subroutine MPIGather_d_v                                                                        
                                                                                                    
subroutine MPIAllGather_d(sendbuf,sendsz,recvbuf,recvsz,type)                                       
   real(dp), intent(in) :: sendbuf                                                                  
   real(dp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather_d                                                                       
subroutine MPIAllGather_d_v(sendbuf,sendsz,recvbuf,recvsz,type)                                     
   real(dp), intent(in) :: sendbuf(*)                                                               
   real(dp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather_d_v                                                                     
subroutine MPIAllGather2_d(sendbuf,sendsz,recvbuf,recvsz,type)                                      
   real(dp), intent(in) :: sendbuf                                                                  
   real(dp), intent(out) :: recvbuf(:,:)                                                            
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather2_d                                                                      
subroutine MPIAllGather2_d_v(sendbuf,sendsz,recvbuf,recvsz,type)                                    
   real(dp), intent(in) :: sendbuf(*)                                                               
   real(dp), intent(out) :: recvbuf(:,:)                                                            
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather2_d_v                                                                    
                                                                                                    
                                                                                                    
subroutine MPIAllGatherV_d(sendbuf,sendsz,recvbuf,recvsz,displs,type)                               
   integer, intent(in) :: sendbuf                                                                   
   real(dp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sendsz, recvsz(*), displs(*)                                              
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)          
   call MPITimer(0)                                                                                 
end subroutine MPIAllGatherV_d                                                                      
subroutine MPIAllGatherV_d_v(sendbuf,sendsz,recvbuf,recvsz,displs,type)                             
   real(dp), intent(in) :: sendbuf(*)                                                               
   real(dp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sendsz, recvsz(*), displs(*)                                              
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)          
   call MPITimer(0)                                                                                 
end subroutine MPIAllGatherV_d_v                                                                    
                                                                                                    
subroutine MPISendRecv_d(sendbuf,sendsz,dest,sendtag,recvbuf,recvsz,source,recvtag,type)            
   real(dp), intent(in) :: sendbuf(*)                                                               
   real(dp), intent(out) :: recvbuf(*)                                                              
   integer, intent(in) :: sendsz, recvsz, dest, source                                              
   integer, intent(in) :: sendtag, recvtag, type                                                    
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_SendRecv(sendbuf,sendsz,type,dest,sendtag,recvbuf,recvsz,type,source,recvtag,&          
     MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)                                                         
   call MPITimer(0)                                                                                 
end subroutine MPISendRecv_d                                                                        
                                                                                                    
end module mpi_d                                                                                    
                                                                                                    
module mpi_s                                                                                        
                                                                                                    
   use mpitime,              only : MPITimer                                                        
   use mpi_inc                                                                                      
   use mpikinds                                                                                     
                                                                                                    
   implicit none                                                                                    
                                                                                                    
   PRIVATE                                                                                          
                                                                                                    
   integer :: ierr                                                                                  
                                                                                                    
   public :: MPIBCast                                                                               
   interface MPIBCast                                                                               
      module procedure MPIBCast_s, MPIBCast_s_v                                                     
   end interface                                                                                    
                                                                                                    
   public :: MPIRedSum                                                                              
   interface MPIRedSum                                                                              
      module procedure MPIRedSum_s,  MPIRedSum_s_v                                                  
   end interface                                                                                    
                                                                                                    
   public :: MPIAllRedSum                                                                           
   interface MPIAllRedSum                                                                           
      module procedure MPIAllRedSum_s,  MPIAllRedSum_s_v                                            
   end interface                                                                                    
                                                                                                    
   public :: MPIGather                                                                              
   interface MPIGather                                                                              
      module procedure MPIGather_s,  MPIGather_s_v                                                  
   end interface                                                                                    
                                                                                                    
   public :: MPIAllGather                                                                           
   interface MPIAllGather                                                                           
      module procedure MPIAllGather_s,  MPIAllGather_s_v                                            
      module procedure MPIAllGather2_s, MPIAllGather2_s_v                                           
   end interface                                                                                    
                                                                                                    
   public :: MPIAllGatherV                                                                          
   interface MPIAllGatherV                                                                          
      module procedure MPIAllGatherV_s, MPIAllGatherV_s_v                                           
   end interface                                                                                    
                                                                                                    
   public :: MPISendRecv                                                                            
   interface MPISendRecv                                                                            
      module procedure MPISendRecv_s                                                                
   end interface                                                                                    
                                                                                                    
contains                                                                                            
                                                                                                    
subroutine MPIBCast_s(buffer,sz,type)                                                               
   character(*), intent(inout) :: buffer                                                            
   integer :: sz                                                                                    
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)                                     
   call MPITimer(0)                                                                                 
end subroutine MPIBCast_s                                                                           
subroutine MPIBCast_s_v(buffer,sz,type)                                                             
   character(*), intent(inout) :: buffer(*)                                                         
   integer :: sz                                                                                    
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)                                     
   call MPITimer(0)                                                                                 
end subroutine MPIBCast_s_v                                                                         
                                                                                                    
subroutine MPIRedSum_s(sendbuf,recvbuf,sz,type)                                                     
   character(*), intent(in) :: sendbuf                                                              
   character(*), intent(out) :: recvbuf                                                             
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)                   
   call MPITimer(0)                                                                                 
end subroutine MPIRedSum_s                                                                          
subroutine MPIRedSum_s_v(sendbuf,recvbuf,sz,type)                                                   
   character(*), intent(in) :: sendbuf(*)                                                           
   character(*), intent(out) :: recvbuf(*)                                                          
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)                   
   call MPITimer(0)                                                                                 
end subroutine MPIRedSum_s_v                                                                        
                                                                                                    
subroutine MPIAllRedSum_s(sendbuf,recvbuf,sz,type)                                                  
   character(*), intent(in) :: sendbuf                                                              
   character(*), intent(out) :: recvbuf                                                             
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)                          
   call MPITimer(0)                                                                                 
end subroutine MPIAllRedSum_s                                                                       
subroutine MPIAllRedSum_s_v(sendbuf,recvbuf,sz,type)                                                
   character(*), intent(in) :: sendbuf(*)                                                           
   character(*), intent(out) :: recvbuf(*)                                                          
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)                          
   call MPITimer(0)                                                                                 
end subroutine MPIAllRedSum_s_v                                                                     
                                                                                                    
subroutine MPIGather_s(sendbuf,sendsz,recvbuf,recvsz,type)                                          
   character(*), intent(in) :: sendbuf                                                              
   character(*), intent(out) :: recvbuf(*)                                                          
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)           
   call MPITimer(0)                                                                                 
end subroutine MPIGather_s                                                                          
subroutine MPIGather_s_v(sendbuf,sendsz,recvbuf,recvsz,type)                                        
   character(*), intent(in) :: sendbuf(*)                                                           
   character(*), intent(out) :: recvbuf(*)                                                          
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)           
   call MPITimer(0)                                                                                 
end subroutine MPIGather_s_v                                                                        
                                                                                                    
subroutine MPIAllGather_s(sendbuf,sendsz,recvbuf,recvsz,type)                                       
   character(*), intent(in) :: sendbuf                                                              
   character(*), intent(out) :: recvbuf(*)                                                          
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather_s                                                                       
subroutine MPIAllGather_s_v(sendbuf,sendsz,recvbuf,recvsz,type)                                     
   character(*), intent(in) :: sendbuf(*)                                                           
   character(*), intent(out) :: recvbuf(*)                                                          
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather_s_v                                                                     
subroutine MPIAllGather2_s(sendbuf,sendsz,recvbuf,recvsz,type)                                      
   character(*), intent(in) :: sendbuf                                                              
   character(*), intent(out) :: recvbuf(:,:)                                                        
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather2_s                                                                      
subroutine MPIAllGather2_s_v(sendbuf,sendsz,recvbuf,recvsz,type)                                    
   character(*), intent(in) :: sendbuf(*)                                                           
   character(*), intent(out) :: recvbuf(:,:)                                                        
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather2_s_v                                                                    
                                                                                                    
                                                                                                    
subroutine MPIAllGatherV_s(sendbuf,sendsz,recvbuf,recvsz,displs,type)                               
   integer, intent(in) :: sendbuf                                                                   
   character(*), intent(out) :: recvbuf(*)                                                          
   integer, intent(in) :: sendsz, recvsz(*), displs(*)                                              
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)          
   call MPITimer(0)                                                                                 
end subroutine MPIAllGatherV_s                                                                      
subroutine MPIAllGatherV_s_v(sendbuf,sendsz,recvbuf,recvsz,displs,type)                             
   character(*), intent(in) :: sendbuf(*)                                                           
   character(*), intent(out) :: recvbuf(*)                                                          
   integer, intent(in) :: sendsz, recvsz(*), displs(*)                                              
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)          
   call MPITimer(0)                                                                                 
end subroutine MPIAllGatherV_s_v                                                                    
                                                                                                    
subroutine MPISendRecv_s(sendbuf,sendsz,dest,sendtag,recvbuf,recvsz,source,recvtag,type)            
   character(*), intent(in) :: sendbuf(*)                                                           
   character(*), intent(out) :: recvbuf(*)                                                          
   integer, intent(in) :: sendsz, recvsz, dest, source                                              
   integer, intent(in) :: sendtag, recvtag, type                                                    
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_SendRecv(sendbuf,sendsz,type,dest,sendtag,recvbuf,recvsz,type,source,recvtag,&          
     MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)                                                         
   call MPITimer(0)                                                                                 
end subroutine MPISendRecv_s                                                                        
                                                                                                    
end module mpi_s                                                                                    
                                                                                                    
module mpi_l                                                                                        
                                                                                                    
   use mpitime,              only : MPITimer                                                        
   use mpi_inc                                                                                      
   use mpikinds                                                                                     
                                                                                                    
   implicit none                                                                                    
                                                                                                    
   PRIVATE                                                                                          
                                                                                                    
   integer :: ierr                                                                                  
                                                                                                    
   public :: MPIBCast                                                                               
   interface MPIBCast                                                                               
      module procedure MPIBCast_l, MPIBCast_l_v                                                     
   end interface                                                                                    
                                                                                                    
   public :: MPIRedSum                                                                              
   interface MPIRedSum                                                                              
      module procedure MPIRedSum_l,  MPIRedSum_l_v                                                  
   end interface                                                                                    
                                                                                                    
   public :: MPIAllRedSum                                                                           
   interface MPIAllRedSum                                                                           
      module procedure MPIAllRedSum_l,  MPIAllRedSum_l_v                                            
   end interface                                                                                    
                                                                                                    
   public :: MPIGather                                                                              
   interface MPIGather                                                                              
      module procedure MPIGather_l,  MPIGather_l_v                                                  
   end interface                                                                                    
                                                                                                    
   public :: MPIAllGather                                                                           
   interface MPIAllGather                                                                           
      module procedure MPIAllGather_l,  MPIAllGather_l_v                                            
      module procedure MPIAllGather2_l, MPIAllGather2_l_v                                           
   end interface                                                                                    
                                                                                                    
   public :: MPIAllGatherV                                                                          
   interface MPIAllGatherV                                                                          
      module procedure MPIAllGatherV_l, MPIAllGatherV_l_v                                           
   end interface                                                                                    
                                                                                                    
   public :: MPISendRecv                                                                            
   interface MPISendRecv                                                                            
      module procedure MPISendRecv_l                                                                
   end interface                                                                                    
                                                                                                    
contains                                                                                            
                                                                                                    
subroutine MPIBCast_l(buffer,sz,type)                                                               
   logical, intent(inout) :: buffer                                                                 
   integer :: sz                                                                                    
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)                                     
   call MPITimer(0)                                                                                 
end subroutine MPIBCast_l                                                                           
subroutine MPIBCast_l_v(buffer,sz,type)                                                             
   logical, intent(inout) :: buffer(*)                                                              
   integer :: sz                                                                                    
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)                                     
   call MPITimer(0)                                                                                 
end subroutine MPIBCast_l_v                                                                         
                                                                                                    
subroutine MPIRedSum_l(sendbuf,recvbuf,sz,type)                                                     
   logical, intent(in) :: sendbuf                                                                   
   logical, intent(out) :: recvbuf                                                                  
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)                   
   call MPITimer(0)                                                                                 
end subroutine MPIRedSum_l                                                                          
subroutine MPIRedSum_l_v(sendbuf,recvbuf,sz,type)                                                   
   logical, intent(in) :: sendbuf(*)                                                                
   logical, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)                   
   call MPITimer(0)                                                                                 
end subroutine MPIRedSum_l_v                                                                        
                                                                                                    
subroutine MPIAllRedSum_l(sendbuf,recvbuf,sz,type)                                                  
   logical, intent(in) :: sendbuf                                                                   
   logical, intent(out) :: recvbuf                                                                  
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)                          
   call MPITimer(0)                                                                                 
end subroutine MPIAllRedSum_l                                                                       
subroutine MPIAllRedSum_l_v(sendbuf,recvbuf,sz,type)                                                
   logical, intent(in) :: sendbuf(*)                                                                
   logical, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)                          
   call MPITimer(0)                                                                                 
end subroutine MPIAllRedSum_l_v                                                                     
                                                                                                    
subroutine MPIGather_l(sendbuf,sendsz,recvbuf,recvsz,type)                                          
   logical, intent(in) :: sendbuf                                                                   
   logical, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)           
   call MPITimer(0)                                                                                 
end subroutine MPIGather_l                                                                          
subroutine MPIGather_l_v(sendbuf,sendsz,recvbuf,recvsz,type)                                        
   logical, intent(in) :: sendbuf(*)                                                                
   logical, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)           
   call MPITimer(0)                                                                                 
end subroutine MPIGather_l_v                                                                        
                                                                                                    
subroutine MPIAllGather_l(sendbuf,sendsz,recvbuf,recvsz,type)                                       
   logical, intent(in) :: sendbuf                                                                   
   logical, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather_l                                                                       
subroutine MPIAllGather_l_v(sendbuf,sendsz,recvbuf,recvsz,type)                                     
   logical, intent(in) :: sendbuf(*)                                                                
   logical, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather_l_v                                                                     
subroutine MPIAllGather2_l(sendbuf,sendsz,recvbuf,recvsz,type)                                      
   logical, intent(in) :: sendbuf                                                                   
   logical, intent(out) :: recvbuf(:,:)                                                             
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather2_l                                                                      
subroutine MPIAllGather2_l_v(sendbuf,sendsz,recvbuf,recvsz,type)                                    
   logical, intent(in) :: sendbuf(*)                                                                
   logical, intent(out) :: recvbuf(:,:)                                                             
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather2_l_v                                                                    
                                                                                                    
                                                                                                    
subroutine MPIAllGatherV_l(sendbuf,sendsz,recvbuf,recvsz,displs,type)                               
   integer, intent(in) :: sendbuf                                                                   
   logical, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sendsz, recvsz(*), displs(*)                                              
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)          
   call MPITimer(0)                                                                                 
end subroutine MPIAllGatherV_l                                                                      
subroutine MPIAllGatherV_l_v(sendbuf,sendsz,recvbuf,recvsz,displs,type)                             
   logical, intent(in) :: sendbuf(*)                                                                
   logical, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sendsz, recvsz(*), displs(*)                                              
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)          
   call MPITimer(0)                                                                                 
end subroutine MPIAllGatherV_l_v                                                                    
                                                                                                    
subroutine MPISendRecv_l(sendbuf,sendsz,dest,sendtag,recvbuf,recvsz,source,recvtag,type)            
   logical, intent(in) :: sendbuf(*)                                                                
   logical, intent(out) :: recvbuf(*)                                                               
   integer, intent(in) :: sendsz, recvsz, dest, source                                              
   integer, intent(in) :: sendtag, recvtag, type                                                    
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_SendRecv(sendbuf,sendsz,type,dest,sendtag,recvbuf,recvsz,type,source,recvtag,&          
     MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)                                                         
   call MPITimer(0)                                                                                 
end subroutine MPISendRecv_l                                                                        
                                                                                                    
end module mpi_l                                                                                    
                                                                                                    
module mpi_c                                                                                        
                                                                                                    
   use mpitime,              only : MPITimer                                                        
   use mpi_inc                                                                                      
   use mpikinds                                                                                     
                                                                                                    
   implicit none                                                                                    
                                                                                                    
   PRIVATE                                                                                          
                                                                                                    
   integer :: ierr                                                                                  
                                                                                                    
   public :: MPIBCast                                                                               
   interface MPIBCast                                                                               
      module procedure MPIBCast_c, MPIBCast_c_v                                                     
   end interface                                                                                    
                                                                                                    
   public :: MPIRedSum                                                                              
   interface MPIRedSum                                                                              
      module procedure MPIRedSum_c,  MPIRedSum_c_v                                                  
   end interface                                                                                    
                                                                                                    
   public :: MPIAllRedSum                                                                           
   interface MPIAllRedSum                                                                           
      module procedure MPIAllRedSum_c,  MPIAllRedSum_c_v                                            
   end interface                                                                                    
                                                                                                    
   public :: MPIGather                                                                              
   interface MPIGather                                                                              
      module procedure MPIGather_c,  MPIGather_c_v                                                  
   end interface                                                                                    
                                                                                                    
   public :: MPIAllGather                                                                           
   interface MPIAllGather                                                                           
      module procedure MPIAllGather_c,  MPIAllGather_c_v                                            
      module procedure MPIAllGather2_c, MPIAllGather2_c_v                                           
   end interface                                                                                    
                                                                                                    
   public :: MPIAllGatherV                                                                          
   interface MPIAllGatherV                                                                          
      module procedure MPIAllGatherV_c, MPIAllGatherV_c_v                                           
   end interface                                                                                    
                                                                                                    
   public :: MPISendRecv                                                                            
   interface MPISendRecv                                                                            
      module procedure MPISendRecv_c                                                                
   end interface                                                                                    
                                                                                                    
contains                                                                                            
                                                                                                    
subroutine MPIBCast_c(buffer,sz,type)                                                               
   complex(sp), intent(inout) :: buffer                                                             
   integer :: sz                                                                                    
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)                                     
   call MPITimer(0)                                                                                 
end subroutine MPIBCast_c                                                                           
subroutine MPIBCast_c_v(buffer,sz,type)                                                             
   complex(sp), intent(inout) :: buffer(*)                                                          
   integer :: sz                                                                                    
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)                                     
   call MPITimer(0)                                                                                 
end subroutine MPIBCast_c_v                                                                         
                                                                                                    
subroutine MPIRedSum_c(sendbuf,recvbuf,sz,type)                                                     
   complex(sp), intent(in) :: sendbuf                                                               
   complex(sp), intent(out) :: recvbuf                                                              
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)                   
   call MPITimer(0)                                                                                 
end subroutine MPIRedSum_c                                                                          
subroutine MPIRedSum_c_v(sendbuf,recvbuf,sz,type)                                                   
   complex(sp), intent(in) :: sendbuf(*)                                                            
   complex(sp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)                   
   call MPITimer(0)                                                                                 
end subroutine MPIRedSum_c_v                                                                        
                                                                                                    
subroutine MPIAllRedSum_c(sendbuf,recvbuf,sz,type)                                                  
   complex(sp), intent(in) :: sendbuf                                                               
   complex(sp), intent(out) :: recvbuf                                                              
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)                          
   call MPITimer(0)                                                                                 
end subroutine MPIAllRedSum_c                                                                       
subroutine MPIAllRedSum_c_v(sendbuf,recvbuf,sz,type)                                                
   complex(sp), intent(in) :: sendbuf(*)                                                            
   complex(sp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)                          
   call MPITimer(0)                                                                                 
end subroutine MPIAllRedSum_c_v                                                                     
                                                                                                    
subroutine MPIGather_c(sendbuf,sendsz,recvbuf,recvsz,type)                                          
   complex(sp), intent(in) :: sendbuf                                                               
   complex(sp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)           
   call MPITimer(0)                                                                                 
end subroutine MPIGather_c                                                                          
subroutine MPIGather_c_v(sendbuf,sendsz,recvbuf,recvsz,type)                                        
   complex(sp), intent(in) :: sendbuf(*)                                                            
   complex(sp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)           
   call MPITimer(0)                                                                                 
end subroutine MPIGather_c_v                                                                        
                                                                                                    
subroutine MPIAllGather_c(sendbuf,sendsz,recvbuf,recvsz,type)                                       
   complex(sp), intent(in) :: sendbuf                                                               
   complex(sp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather_c                                                                       
subroutine MPIAllGather_c_v(sendbuf,sendsz,recvbuf,recvsz,type)                                     
   complex(sp), intent(in) :: sendbuf(*)                                                            
   complex(sp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather_c_v                                                                     
subroutine MPIAllGather2_c(sendbuf,sendsz,recvbuf,recvsz,type)                                      
   complex(sp), intent(in) :: sendbuf                                                               
   complex(sp), intent(out) :: recvbuf(:,:)                                                         
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather2_c                                                                      
subroutine MPIAllGather2_c_v(sendbuf,sendsz,recvbuf,recvsz,type)                                    
   complex(sp), intent(in) :: sendbuf(*)                                                            
   complex(sp), intent(out) :: recvbuf(:,:)                                                         
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather2_c_v                                                                    
                                                                                                    
                                                                                                    
subroutine MPIAllGatherV_c(sendbuf,sendsz,recvbuf,recvsz,displs,type)                               
   integer, intent(in) :: sendbuf                                                                   
   complex(sp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sendsz, recvsz(*), displs(*)                                              
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)          
   call MPITimer(0)                                                                                 
end subroutine MPIAllGatherV_c                                                                      
subroutine MPIAllGatherV_c_v(sendbuf,sendsz,recvbuf,recvsz,displs,type)                             
   complex(sp), intent(in) :: sendbuf(*)                                                            
   complex(sp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sendsz, recvsz(*), displs(*)                                              
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)          
   call MPITimer(0)                                                                                 
end subroutine MPIAllGatherV_c_v                                                                    
                                                                                                    
subroutine MPISendRecv_c(sendbuf,sendsz,dest,sendtag,recvbuf,recvsz,source,recvtag,type)            
   complex(sp), intent(in) :: sendbuf(*)                                                            
   complex(sp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sendsz, recvsz, dest, source                                              
   integer, intent(in) :: sendtag, recvtag, type                                                    
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_SendRecv(sendbuf,sendsz,type,dest,sendtag,recvbuf,recvsz,type,source,recvtag,&          
     MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)                                                         
   call MPITimer(0)                                                                                 
end subroutine MPISendRecv_c                                                                        
                                                                                                    
end module mpi_c                                                                                    
                                                                                                    
module mpi_z                                                                                        
                                                                                                    
   use mpitime,              only : MPITimer                                                        
   use mpi_inc                                                                                      
   use mpikinds                                                                                     
                                                                                                    
   implicit none                                                                                    
                                                                                                    
   PRIVATE                                                                                          
                                                                                                    
   integer :: ierr                                                                                  
                                                                                                    
   public :: MPIBCast                                                                               
   interface MPIBCast                                                                               
      module procedure MPIBCast_z, MPIBCast_z_v                                                     
   end interface                                                                                    
                                                                                                    
   public :: MPIRedSum                                                                              
   interface MPIRedSum                                                                              
      module procedure MPIRedSum_z,  MPIRedSum_z_v                                                  
   end interface                                                                                    
                                                                                                    
   public :: MPIAllRedSum                                                                           
   interface MPIAllRedSum                                                                           
      module procedure MPIAllRedSum_z,  MPIAllRedSum_z_v                                            
   end interface                                                                                    
                                                                                                    
   public :: MPIGather                                                                              
   interface MPIGather                                                                              
      module procedure MPIGather_z,  MPIGather_z_v                                                  
   end interface                                                                                    
                                                                                                    
   public :: MPIAllGather                                                                           
   interface MPIAllGather                                                                           
      module procedure MPIAllGather_z,  MPIAllGather_z_v                                            
      module procedure MPIAllGather2_z, MPIAllGather2_z_v                                           
   end interface                                                                                    
                                                                                                    
   public :: MPIAllGatherV                                                                          
   interface MPIAllGatherV                                                                          
      module procedure MPIAllGatherV_z, MPIAllGatherV_z_v                                           
   end interface                                                                                    
                                                                                                    
   public :: MPISendRecv                                                                            
   interface MPISendRecv                                                                            
      module procedure MPISendRecv_z                                                                
   end interface                                                                                    
                                                                                                    
contains                                                                                            
                                                                                                    
subroutine MPIBCast_z(buffer,sz,type)                                                               
   complex(dp), intent(inout) :: buffer                                                             
   integer :: sz                                                                                    
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)                                     
   call MPITimer(0)                                                                                 
end subroutine MPIBCast_z                                                                           
subroutine MPIBCast_z_v(buffer,sz,type)                                                             
   complex(dp), intent(inout) :: buffer(*)                                                          
   integer :: sz                                                                                    
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)                                     
   call MPITimer(0)                                                                                 
end subroutine MPIBCast_z_v                                                                         
                                                                                                    
subroutine MPIRedSum_z(sendbuf,recvbuf,sz,type)                                                     
   complex(dp), intent(in) :: sendbuf                                                               
   complex(dp), intent(out) :: recvbuf                                                              
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)                   
   call MPITimer(0)                                                                                 
end subroutine MPIRedSum_z                                                                          
subroutine MPIRedSum_z_v(sendbuf,recvbuf,sz,type)                                                   
   complex(dp), intent(in) :: sendbuf(*)                                                            
   complex(dp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)                   
   call MPITimer(0)                                                                                 
end subroutine MPIRedSum_z_v                                                                        
                                                                                                    
subroutine MPIAllRedSum_z(sendbuf,recvbuf,sz,type)                                                  
   complex(dp), intent(in) :: sendbuf                                                               
   complex(dp), intent(out) :: recvbuf                                                              
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)                          
   call MPITimer(0)                                                                                 
end subroutine MPIAllRedSum_z                                                                       
subroutine MPIAllRedSum_z_v(sendbuf,recvbuf,sz,type)                                                
   complex(dp), intent(in) :: sendbuf(*)                                                            
   complex(dp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sz                                                                        
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)                          
   call MPITimer(0)                                                                                 
end subroutine MPIAllRedSum_z_v                                                                     
                                                                                                    
subroutine MPIGather_z(sendbuf,sendsz,recvbuf,recvsz,type)                                          
   complex(dp), intent(in) :: sendbuf                                                               
   complex(dp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)           
   call MPITimer(0)                                                                                 
end subroutine MPIGather_z                                                                          
subroutine MPIGather_z_v(sendbuf,sendsz,recvbuf,recvsz,type)                                        
   complex(dp), intent(in) :: sendbuf(*)                                                            
   complex(dp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)           
   call MPITimer(0)                                                                                 
end subroutine MPIGather_z_v                                                                        
                                                                                                    
subroutine MPIAllGather_z(sendbuf,sendsz,recvbuf,recvsz,type)                                       
   complex(dp), intent(in) :: sendbuf                                                               
   complex(dp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather_z                                                                       
subroutine MPIAllGather_z_v(sendbuf,sendsz,recvbuf,recvsz,type)                                     
   complex(dp), intent(in) :: sendbuf(*)                                                            
   complex(dp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather_z_v                                                                     
subroutine MPIAllGather2_z(sendbuf,sendsz,recvbuf,recvsz,type)                                      
   complex(dp), intent(in) :: sendbuf                                                               
   complex(dp), intent(out) :: recvbuf(:,:)                                                         
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather2_z                                                                      
subroutine MPIAllGather2_z_v(sendbuf,sendsz,recvbuf,recvsz,type)                                    
   complex(dp), intent(in) :: sendbuf(*)                                                            
   complex(dp), intent(out) :: recvbuf(:,:)                                                         
   integer, intent(in) :: sendsz, recvsz                                                            
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)                  
   call MPITimer(0)                                                                                 
end subroutine MPIAllGather2_z_v                                                                    
                                                                                                    
                                                                                                    
subroutine MPIAllGatherV_z(sendbuf,sendsz,recvbuf,recvsz,displs,type)                               
   integer, intent(in) :: sendbuf                                                                   
   complex(dp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sendsz, recvsz(*), displs(*)                                              
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)          
   call MPITimer(0)                                                                                 
end subroutine MPIAllGatherV_z                                                                      
subroutine MPIAllGatherV_z_v(sendbuf,sendsz,recvbuf,recvsz,displs,type)                             
   complex(dp), intent(in) :: sendbuf(*)                                                            
   complex(dp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sendsz, recvsz(*), displs(*)                                              
   integer, intent(in) :: type                                                                      
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)          
   call MPITimer(0)                                                                                 
end subroutine MPIAllGatherV_z_v                                                                    
                                                                                                    
subroutine MPISendRecv_z(sendbuf,sendsz,dest,sendtag,recvbuf,recvsz,source,recvtag,type)            
   complex(dp), intent(in) :: sendbuf(*)                                                            
   complex(dp), intent(out) :: recvbuf(*)                                                           
   integer, intent(in) :: sendsz, recvsz, dest, source                                              
   integer, intent(in) :: sendtag, recvtag, type                                                    
                                                                                                    
   call MPITimer(1)                                                                                 
   call MPI_SendRecv(sendbuf,sendsz,type,dest,sendtag,recvbuf,recvsz,type,source,recvtag,&          
     MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)                                                         
   call MPITimer(0)                                                                                 
end subroutine MPISendRecv_z                                                                        
                                                                                                    
end module mpi_z                                                                                    
