!      Ghanim Ullah, John R Cressman Jr, Ernest Barreto, and Steven J Schiff,
!                        July 29, 2008,
!                        Penn State University,
!                        Email: ghanim.phy@gmail.com
!      Archived with 
!      "The Influence of Sodium and Potassium Dynamics on Excitability, Seizures, 
!      and the Stability of Persistent States: II. Network and Glia Dynamics" (2009) 
!      Journal of Computational Neuroscience, 26:171-183
!      This program averages the activity of the network over larger bin size. 
!      The program takes data from the main program Network.f90, and uses larger bin
!      size for the final result. The output from the main program Network.f90, used here is the
!      activity of the excitatory (inhibitory) network.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
       implicit none   
       integer i,j,n, n1, n2
       PARAMETER (n=400000, n1=100)                                          !n=number of data points in the network activity data file from Network.f90
       double precision,ALLOCATABLE,DIMENSION(:):: data1, data2              !n1=bin size, can be changed according to the desired bin size
       double precision average_activity, average_time                       
       integer status
       open(unit=1,file='activity.dat')                                 !activity.dat stores the output for the results in the paper
       open(33, file='activity_excitatory.dat', status='unknown')      !output file from the main program Network.f90
       open(34, file='time.dat', status='unknown')           !Time from the main program Network.f90
       ALLOCATE(data1(n), STAT=status)
       ALLOCATE(data2(n), STAT=status)
       do i = 1, n                     
          read(33,*)data1(i)                                                 !reads ex_activity into data1 array
          read(34,*)data2(i)                                                 !reads time into data2 array 
       end do
       do  j = 1, n/n1
             average_activity = 0.0d0                   
             average_time = 0.0d0
             do i = 1, n1                                                    !this loop performs the average over larger bin size 
                n2 = (j-1)*n1+i
                average_activity = average_activity + data1(n2)
                average_time = average_time+data2(n2)				
             end do
             write(1,*)average_time/float(n1),average_activity/float(n1)     !output saved in activity.dat   
       end do
       DEALLOCATE(data1, STAT=status)
       DEALLOCATE(data2, STAT=status)
       end
