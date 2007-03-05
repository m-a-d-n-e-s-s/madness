! Copyright (c) 2006, R.J. Harrison, UT/ORNL
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
! 
!     * Redistributions of source code must retain the above copyright 
!       notice, this list of conditions and the following disclaimer.
! 
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution.
! 
!     * Neither the names of Oak Ridge National Laboratory and the
!       University of Tennessee, nor the names of its contributors may
!       be used to endorse or promote products derived from this
!       software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

        module test1_mod
        implicit none

!       -----------------------------------------
!       use a structure to pass in more arguments
!       -----------------------------------------
        type hello_args
           integer myid
           integer ndim
           integer ival
           character*80 name
        end type hello_args

        contains

!       ---------------------------------------------
!       use recursive key word
!       to encourage local variables be allocated
!       on the stack
!       may need to use -Mrecursive compiler options
!       ---------------------------------------------


        recursive subroutine hello(arg_in)
        implicit none
        type(hello_args) :: arg_in

        integer ndim
        integer, allocatable, dimension(:) :: A
        integer myid,ival,info
        character*80 name
!       -----------------------------------------------
!       do not use the "save" statement unless all
!       fibers are supposed to share this variable
!       perhaps use common or module variable is better
!       -----------------------------------------------

!       need to make a copy, argument may be shared
!       by all fibers
        myid = arg_in%myid
        ndim = arg_in%ndim
        ival = arg_in%ival
        name = arg_in%name

        allocate( A(ndim), stat=info )
        if (info.ne.0) then
            write(*,*) 'allocate(A), info,myid, name',                    &
     &             info,myid,trim(name)
            stop '** error ** '
        endif

        call set_vec(ndim, A, ival)
        write(*,*) 'starting, myid, name ',myid,trim(name)
        call fiber_yield()


        call check_vec(ndim,A,ival)
        write(*,*) '2nd step, myid ',myid
        call fiber_yield()

        call check_vec(ndim,A,ival)
        write(*,*) '3rd step, myid ',myid
        call fiber_yield()

        call check_vec(ndim,A,ival)
        write(*,*) 'Goodbye from myid ',myid

        deallocate(A,stat=info)
        if (info.ne.0) then
           write(*,*) 'deallocate(A),info,myid,name',info,myid,name
           stop '** error ** '
        endif

        return
        end subroutine hello

        subroutine set_vec(n,A,ival)
        implicit none
        integer n,A(n),ival
        integer i
        do i=1,n
            A(i) = ival
        enddo
        return
        end subroutine set_vec


        subroutine check_vec(n, A,ival)
        integer n,A(n),ival
        integer i
        do i=1,n
          if (A(i).ne.ival) then
             write(*,*) 'bad vector, i,ival,A(i) ',i,ival,A(i)
             stop '** error **'
          endif
        enddo
        return
        end subroutine check_vec

        recursive subroutine noop(idummy)
        implicit none
        integer idummy
        return
        end subroutine noop

        recursive subroutine countdown(vn)
        implicit none
        integer vn
        do while (vn .ne. 0)
           vn = vn - 1
           call fiber_yield()
        enddo
        return
        end subroutine countdown

        end module test1_mod



        program test1
        use test1_mod
        implicit none
!
!       simple fortran program to test fiber package
!
        integer vn
        integer t1,t2,count_rate
        real wall_used, cpu_used, dt1,dt2
        integer n,i,isize,info
        integer fiber_count, fiber_create
        external fiber_count

        integer, parameter :: nfibers = 6
        type(hello_args) :: args(nfibers)

        write(*,*) 'main program started'

        isize = 64*1024+1
        call fiber_set_stacksize(isize)
        do i=1,nfibers
          args(i)%myid = i
          args(i)%ndim = i*1024
          args(i)%ival = -i
          args(i)%name = "worker"
        enddo

        args(1)%name = "Dave"
        args(2)%name = "Mary"
        args(3)%name = "Fred"
        args(4)%name = "Nick"

        do i=1,nfibers
          info =  fiber_create(hello,args(i))

          write(*,*) 'i, info ',i,info
        enddo


        write(*,*) 'after fibers created'
        write(*,*) 'fiber_count() is ',fiber_count()
        call fiber_yield()

        do while (fiber_count().ne.0)
          write(*,*) 'Main thread yielding, fiber_count() ',                &
     &                fiber_count()
          call fiber_yield()
        enddo


        call system_clock(t1,count_rate)
        call cpu_time(cpu_used)

        n = 1000*1000
        do while (.true.)
          call system_clock(t1,count_rate)
          wall_used = real(t1)/real(count_rate)
          call cpu_time( dt1 )
          cpu_used = dt1
          do i=1,n
            info =  fiber_create(noop,0)
            call fiber_yield()
          enddo
          call system_clock(t2,count_rate)
          wall_used = real(t2)/real(count_rate) - wall_used
          call cpu_time(   dt2 )
          cpu_used = dt2 - cpu_used

          write(*,*) 'n, wall, cpu ',n, wall_used, cpu_used
          if (wall_used > 1.0) then
               exit
          else if (wall_used > 0.01) then
            n = 1.04 * n/wall_used
          else
            n = n * 10
          endif

         enddo
         write(*,*) 'time per create+yield+return+destroy '
         write(*,*) 'wall_used, cpu_used ', wall_used/n, cpu_used/n

         vn = n
         info =  fiber_create( countdown, vn )
         call system_clock(t1,count_rate)
         wall_used = real(t1)/real(count_rate)
         call cpu_time(dt1)
         cpu_used = dt1

         do while (vn .ne. 0) 
           call fiber_yield()
         enddo

         call system_clock(t2,count_rate)
         wall_used = real(t2)/real(count_rate) - wall_used
         call cpu_time(dt2)
         cpu_used = dt2 - cpu_used

         write(*,*) 'time per yield+countdown+yield '
         write(*,*) 'wall_used, cpu_used ',wall_used/n,cpu_used/n



        stop 'all done'
        end program test1
