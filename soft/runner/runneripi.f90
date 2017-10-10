      module runneripi 
!SK: ipistring is the hostaddress: host name or IP address; if you want
!to use UNIX sockets, add 'UNIX' before the host adress
! desc_str will take care of UNIX 
!SK: ipisocket should contain the number of the port that the
!ipi-wrapper is listening to

      implicit none

      character*40 :: desc_str
      character*40 :: ipisocket
      character*40 :: ipistring 
 
      logical luseipi

      end module runneripi 

! changes SK: ipistring - charachter length; ipisocket - type character
