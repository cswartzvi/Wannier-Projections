!*****************************************************************************
Program wanproj
!*****************************************************************************
   !
   use omp_lib
   !
   implicit none
   !
   include 'mpif.h'
   !
   integer, parameter   :: DP  = SELECTED_REAL_KIND(15,99)
   real(DP),parameter   :: ao  = 0.52917721_DP
   !
   complex(DP),allocatable :: read_1d(:),          &  !Initial 1D array for Wavefunction
                              wan_wf(:,:,:),       &  !ix,iy,iz point for Wannier Wavefunction
                              ks_wf(:,:,:)            !ix,iy,iz point for "Kohn-Sham" wavefunction

   real(DP),allocatable ::    space(:,:,:,:),      &  !(ix, iy, iz, point value)  from FFT Grid
                              val(:),              &  !complete overlap
                              val_wan(:),          &  !wannier Functions sum- check
                              val_ks(:)               !Kohn Sham sum - check
   !
   real(DP)                :: valt,                &  !"local" copy of overlap
                              valt_wan,            &  !local copy of wannier sum
                              valt_ks,             &  !local copy of Kohn-Sham sum
                              ans,                 &  !transfer value amoung processors
                              alat(3),             &  !lattice vetors
                              dr(3),               &  !interval distance between space points
                              start_time,          &
                              end_time
                           
   !
   integer                 :: i, j, k,             &      
                              ix, iy, iz,          &  !index for space/wavefunction
                              grid(3), grid_tot,   &  !Grid values and grid totals
                              nbsp,                &  !Number of States
                              ks_state,            &  !Kohn-Sham state to calculate the overlap
                              ks_offset,           &  !added value for the KS fort.* file
                              wan_offset,          &  !added value for the wannier fort.* file
                              ierr
   !
   logical              :: print_xsf
   !
   character(len=100)   :: wan_root, ks_root, tempfile, x1, output, atomfile
   character(len=10)    :: ks_ext, wan_ext
   !
   !----------------------------------------------
   !For parallel execution
   !----------------------------------------------
   integer                 :: myid,       &  !Absolute Index of the processor 0 .... (nproc -1)
                              me,         &  !relative index 
                              nproc,      &  !total number of processors
                              per_proc,   &  !Operations per processor (ideal)
                              remain,     &  !difference in per_proc and nproc
                              np,         &  !processor index
                              status(MPI_STATUS_SIZE)
   !----------------------------------------------
   !
   !Namelist for Standard Input
   NAMELIST /wan_projections/ wan_root, ks_root, wan_offset, grid, nbsp, alat, ks_state, &
                      ks_offset, ks_ext, wan_ext, output, print_xsf, atomfile
   ! 
   !
   !
   !Start Parallel environment Allocate/Initialize variables and read in KS Wavefunction
   Call startup()
   !
   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
   !Loop overall the different wannier Function Files
   !-------------------------------------------------------
   proc_loop: do np=1,per_proc,1
      !
      me = myid + (np-1) * nproc + 1
      if (me > nbsp) exit
      !
      !Determine and open the correct file
      write(x1, '(I0)'), me + wan_offset 
      tempfile = TRIM(wan_root)//TRIM(x1)//TRIM(wan_ext)
      open(unit=1, file=TRIM(tempfile),iostat=ierr, form='unformatted')
      if (ierr /= 0) then
         print *, 'Error: Opening ', tempfile
         stop
      endif
      !
      !Read binary file and create the space and wavefunction files AND convert them to 3D
      read(1) read_1d
      Call array_1D_to_3D (read_1d, wan_wf, grid, grid_tot)
      close(1)
      !
      !If print_xsf if defined print out the xsf file
      if (print_xsf) then
         !
         write(x1, '(I0)'), wan_offset+me
         tempfile = TRIM(wan_root)//TRIM(x1) // '.xsf'
         !
         call write_xsf(dble(wan_wf)**2, grid, tempfile, atomfile)
         !
      endif
      !
      !Loop over all points add up sums
      valt = 0.0_DP
      valt_wan = 0.0_DP
      if(myid == 0) valt_ks  = 0.0_DP
      !
      do iz=1, grid(3)
         do iy=1,grid(2)
            do ix=1,grid(1)
               !
               !SPACE ARRAY
               space(ix,iy,iz,1:3) = (/(ix-1)*dr(1), (iy-1)*dr(2), (iz-1)*dr(3) /)
               !
               !Calculate Overlap
               !valt = valt + ks_wf(ix,iy,iz)*ks_wf(ix,iy,iz) * wan_wf(ix,iy,iz)*wan_wf(ix,iy,iz)
               valt = valt + abs(ks_wf(ix,iy,iz)) * abs(wan_wf(ix,iy,iz))
               !
               !Check wan norm
               valt_wan = valt_wan + wan_wf(ix,iy,iz)*wan_wf(ix,iy,iz)
               !
               !Check ks norm
               if (myid == 0) valt_ks = valt_ks + ks_wf(ix,iy,iz)*ks_wf(ix,iy,iz)
               !
            enddo
         enddo
      enddo
      !
      !Calculate Integrals
      valt = dble(valt/(grid(1)*grid(2)*grid(3)))**2
      valt_wan = valt_wan/(grid(1)*grid(2)*grid(3))
      if (myid ==0) valt_ks = valt_ks/(grid(1)*grid(2)*grid(3))
      !
      !Last loop Clean up
      if (np == per_proc) then
         deallocate( wan_wf, read_1d, space )
         if(allocated(ks_wf)) deallocate(ks_wf)
      endif
      !
      !--------------------------------------------
      ! Communication
      !--------------------------------------------
      Call MPI_Barrier(MPI_COMM_WORLD,ierr)
      Call mpi_send_recv(me, myid, np, nproc, per_proc, remain, val, valt)
      !
      !
      write(6,'(6X,"Wannier Function : ", I4, ", Integrated Chagre : ",F14.5)') me, valt_wan  
      !
   enddo proc_loop 
   !-------------------------------------------------------
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   !
   !
   !-------------------------------------------------------
   Call MPI_Barrier(MPI_COMM_WORLD,ierr)
   !
   if (ierr /= 0 ) then
      write(*,*) ' ERROR: MPI_BARRIER'
   endif
   !
   if(myid == 0)then
      !
      !Open Main output file
      open(unit=3, file=output, IOSTAT=ierr, STATUS='UNKNOWN')
      if (ierr /= 0 ) then
         write (*,*) '  ERROR: There was some issue opeing file: ', output
         stop
      endif
      !
      do i=1,nbsp,1
         !Final package of the overlap
         write(3,'(I8,3X,F14.5)') i, val(i) 
      enddo
      !
      end_time = MPI_Wtime()
      close(3)
      write(6,*)
      write(6,'(5X,"Kohn-Sham Integrated Chagre : ",F14.5)') valt_ks
      write(6,*)
      Write(6,*) '    Total time: ', end_time - start_time
      write(6,*)
      !
   endif
   !
   Call MPI_FINALIZE(ierr)
   !
   !
   !
   !*****************************************************************************
   Contains
   !*****************************************************************************
      !
      !-------------------------------------------------------
      Subroutine startup()
         !
         implicit none
         !
         integer           :: temp
         !
         Call MPI_INIT(ierr)
         Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
         Call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
         if (ierr /= 0 ) then
            print *, '  MPI_INIT error ', myid
         endif
         !
         !Defaults
         atomfile = ' '
         print_xsf = .FALSE.
         wan_offset = 0
         ks_offset = 0
         !
         !Clock
         if (myid == 0) then
            !
            start_time = MPI_Wtime()
            !
            read(5,nml=wan_projections,IOSTAT=ierr)
            !
            call welcome()
            !
         endif
         !
         Call MPI_BCAST(wan_root, 100, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
         Call MPI_BCAST(wan_ext, 100, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
         Call MPI_BCAST(ks_root, 100, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
         Call MPI_BCAST(ks_ext, 100, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
         Call MPI_BCAST(ks_state, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         Call MPI_BCAST(wan_offset, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         Call MPI_BCAST(ks_offset, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         Call MPI_BCAST(alat, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         Call MPI_BCAST(nbsp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         Call MPI_BCAST(grid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         Call MPI_BCAST(alat, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         Call MPI_BCAST(print_xsf, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
         Call MPI_BCAST(atomfile, 100, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
         !
         if (ierr /= 0 ) then
            print *, '  MPI_BCAST error ', myid
         endif
         !
         if (print_xsf .and. atomfile == ' ' ) then 
            if (myid == 0) write(*,*) ' ERORR: atomfile must be defined with print_xsf '
         endif
         !
         !******Number jobs per processor (important for the mpi_send_recv)*****
         per_proc = Ceiling(nbsp/real(nproc))
         remain = per_proc * nproc - nbsp
         !
         !wavefunctions and space grids
         grid_tot = grid(1)*grid(2)*grid(3)
         allocate(read_1d(grid_tot))
         allocate( wan_wf(grid(1), grid(2), grid(3))    )
         allocate( ks_wf(grid(1), grid(2), grid(3))     )
         allocate( space(grid(1), grid(2), grid(3), 3) )
         allocate( val(nbsp)  )
         if (myid == 0) allocate( val_ks(nbsp)   )
         allocate( val_wan(nbsp)  )
         !
         !Initialize
         read_1d(:)     = 0.0_DP
         wan_wf(:,:,:)  = 0.0_DP
         ks_wf(:,:,:)   = 0.0_DP
         space(:,:,:,:) = 0.0_DP
         val(:)         = 0.0_DP
         val_wan(:)     = 0.0_DP
         if (myid == 0) val_ks(:)      = 0.0_DP
         !
         !Set up space with matching indexes 
         do i=1,3,1
            alat(i) = alat(i) *  ao
            dr(i) = alat(i)/DBLE(grid(i) - 1)
         enddo
         !
         !Open the Kohn-Sham Eigenstate File
         if(myid == 0) then
            !
            !Determine and open the correct file
            write(x1, '(I0)'), ks_offset+ks_state
            tempfile = TRIM(ks_root)//TRIM(x1) // TRIM(ks_ext)
            open(unit=2, file=TRIM(tempfile),iostat=ierr, status='unknown', form='UNFORMATTED')
            !
            if (ierr /= 0) then
               print *, 'Error: Opening ', tempfile
               stop
            else 
               write(*,'(/,5X,"Kohn-Sham state read from:",A20)') tempfile
            endif
            !
            read(2) read_1d
            !
            close(2)
            !
         endif
         !
         Call MPI_BCAST(read_1d, grid_tot, MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
         Call array_1D_to_3D (read_1d, ks_wf, grid, grid_tot)
         !
         !If the print_xsf is set then print out the .xsf file for the KS State
         if(myid == 0) then
            if (print_xsf) then
               !
               write(x1, '(I0)'), ks_offset+ks_state
               tempfile = TRIM(ks_root)//TRIM(x1) // '.xsf'
               !
               call write_xsf(dble(ks_wf)**2, grid, tempfile, atomfile)
               !
               write(*,'(5X,"Kohn-Sham xsf :",A20)') tempfile
               !
               close(3)
            endif
            write(*,*)
            !
         endif
         !
         Call MPI_Barrier(MPI_COMM_WORLD,ierr)
         !
         return
         !
         !
      End Subroutine startup
      !-------------------------------------------------------
      !
      !-------------------------------------------------------
      ! Subroutine for all non-root processors to send
      ! thier values to the root processor
      ! The number of cycles per processor is detemined by (in start-up):
      ! 
      !   per_proc = Ceiling(nbsp/real(nproc))
      !   remain = per_proc * nproc - nbsp
      !
      !-------------------------------------------------------
      Subroutine mpi_send_recv(me, myid, np, nproc, per_proc, remain, val_array, val_temp)
         !
         implicit none
         !
         integer,intent(in)       :: me,myid,np,nproc,per_proc,remain
         real(DP),intent(inout)   :: val_array(:) 
         real(DP),intent(inout)   :: val_temp !
         real(DP)                :: ans
         !
         if (myid /= 0) then
            Call MPI_SEND(val_temp, 1, MPI_DOUBLE_PRECISION, 0, me, MPI_COMM_WORLD, ierr)
         endif
         !
         !Receive by Root Process
         if(myid == 0 ) then
            !
            val_array(me) = val_temp
            !
            if (np == per_proc) then !For the last loop, may have unfilled procs
               !
               do i =1,(nproc-1-remain),1 !from the other processors
                  !
                  Call MPI_RECV(ans, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, &
                                 MPI_COMM_WORLD, status, ierr)
                  val_array(status(MPI_TAG)) = ans
                  !
               enddo
               !
            else
               !
               do i =1,(nproc-1),1 !from the other processors
                  Call MPI_RECV(ans, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, &
                                 MPI_COMM_WORLD, status, ierr)
                  val_array(status(MPI_TAG)) = ans
                  !
               enddo
               !
            endif
         endif
         !
         return
         !
      End Subroutine mpi_send_recv
      !-------------------------------------------------------
      !
      !-------------------------------------------------------
      ! Converts a 1D complex array into a 3D complex array
      !-------------------------------------------------------
      Subroutine array_1D_to_3D(array1, array3, grid, grid_tot )
         !
         implicit none
         !
         complex(DP),intent(in)  :: array1(:)
         complex(DP),intent(out) :: array3(:,:,:) 
         integer, intent(in)  :: grid(3), grid_tot       !grid dimensions and totals

         integer              :: i,                   &
                                 loopx,loopy,loopz       !index for 1D -> 3D array conversion
         !
         !
         !Convert From temp (1D) vector to rho (3D) Array 
         !Quantum Espresso outputs the DATA file as 1D vectors, converts to a
         !3D array in chegens.f90
         do i=1,grid_tot
            loopz = INT( ABS(i - 1)/(grid(1)*grid(2)) ) + 1
            loopy = INT( ABS( (i - 1) - (loopz - 1)*grid(1)*grid(2) ) / grid(1) ) + 1
            loopx = i - (loopz-1)*grid(1)*grid(2) - (loopy-1) * grid(1)
            array3(loopx, loopy, loopz) =  array1(i)
         enddo
         !
         return
         !
      End Subroutine array_1D_to_3D
      !-------------------------------------------------------
      !
      !-------------------------------------------------------
      ! Write the KS and WAN Charge Densities (requries atomfile)
      !
      ! Format of atomfile (Note: bohr_or_Ang = 0 for bohr and 1 for Ang ):
      !  
      !  total-number-of-atoms   total-number-of-species  bohr_or_Ang
      !  alat11      alat12      alat13 
      !  alat21      alat22      alat23 
      !  alat31      alat32      alat33 
      !  atomic-number1    pos11    pos12    pos13   
      !  atomic-number1    pos21    pos22    pos23   
      !  atomic-number1    pos31    pos32    pos33   
      !  .
      !  .
      !  .
      !  atomic-numberN    posN1    posN2    posN3   
      !  
      !-------------------------------------------------------
      Subroutine write_xsf(rho, grid, file_name, atomfile)
         !
         implicit none
         !
         real(DP), intent(in)                      :: rho(:,:,:)
         integer,intent(in)                        :: grid(:)
         character(len=*), intent(in)              :: file_name
         character(len=*), intent(in)    :: atomfile
         !
         real(DP), allocatable   :: tau(:,:)
         real(DP)                :: at(3,3)
         !
         integer, allocatable :: ityp(:) 
         integer              :: nat, ntyp, bohr_or_Ang
         integer              :: xunit, aunit
         !
         integer              ::ix, iy, iz, i, j
         !
         open(unit=xunit, file=TRIM(file_name), status='unknown')
         open(unit=aunit, file=TRIM(atomfile), status='unknown')
         !
         !Read the atomfile 
         read(aunit,*) nat, ntyp, bohr_or_Ang
         read(aunit,*) at(1,1:3)
         read(aunit,*) at(2,1:3)
         read(aunit,*) at(3,1:3)
         allocate(ityp(nat), tau(3,nat))
         do i=1,nat
            read(aunit,*) ityp(i), (tau(j, i), j=1,3)
         enddo
         !
         close(aunit)
         if( bohr_or_Ang == 0) then
            tau(:,:) = tau(:,:)*ao
            at(:,:)  = at(:,:)*ao
         endif
         !
         !Print the Files
         WRITE(xunit,*) 'CRYSTAL'
         WRITE(xunit,*) 'PRIMVEC'
         WRITE(xunit,'(2(3f15.9/),3f15.9)') at

         ! write atomic coordinates (and forces)
         WRITE(xunit,*) 'PRIMCOORD'
         WRITE(xunit,*) nat, 1
         DO i = 1, nat
            WRITE (xunit,'(i3,3x,3f15.9,1x,3f12.5)') ityp(i), (tau(j, i), j=1,3)
         END DO

         WRITE(xunit,'(a)') 'BEGIN_BLOCK_DATAGRID_3D'
         WRITE(xunit,'(a)') '3D_PWSCF'
         WRITE(xunit,'(a)') 'DATAGRID_3D_UNKNOWN'

         WRITE(xunit,*) grid(1), grid(2), grid(3)

         ! origin
         WRITE(xunit,'(3f10.6)') 0.0d0, 0.0d0, 0.0d0
         ! lattice vectors
         WRITE(xunit,'(3f10.6)') ((at(i, j), i=1,3), j=1,3)
         ! charge density TWO OPTIONS
         !-------------------------------------------------------------------------
         WRITE(xunit,'(6e13.5)') &
            (((rho(ix, iy, iz), ix=1,grid(1)), iy=1,grid(2)), iz=1,grid(3))
         !-------------------------------------------------------------------------

         WRITE(xunit,'(a)') 'END_DATAGRID_3D'
         WRITE(xunit,'(a)') 'END_BLOCK_DATAGRID_3D'
         !
         deallocate(ityp, tau)
         Close(xunit)
         !
         RETURN
      END SUBROUTINE write_xsf
      !-------------------------------------------------------
      !
      !-------------------------------------------------------
      ! Welcome Print Subroutine
      !-------------------------------------------------------
      Subroutine welcome()
         !
         implicit none
         !
         write(*,*) ' ' 
         write(*,'(2X,"-----------------------------------------------")')
         write(*,'(2X,"             Wannier Projections               ")')
         write(*,'(2X,"-----------------------------------------------")')
         write(*,*)
         write(*,'(3X,"Input Values:")')
         write(*,'(3X,"Kohn Sham State  ", I4)') ks_state
         write(*,'(3X,"Wannier States   ", I4)') nbsp
         write(*,'(3X,"Grid values      ", 1X, I4, I4, I4)') grid(1), grid(2), grid(3)
         write(*,'(3X,"Alat values      ", 3X, 3(F7.4,2X))') alat(1), alat(2), alat(3)
         write(*,*)
         write(*,'(3X,"Output File :      ",A20)') output
         write(*,*)
         if (print_xsf) Then
            write(*,'(3X,"ALL Charge Densities will be converted to xsf Files")')
            write(*,*)
         endif
         write(*,'(3X,"Total Number of Processors: ", I4)') nproc
         write(*,*)
         write(*,'(3X,"Program Start...")')
         write(*,*)
      End Subroutine Welcome
      !-------------------------------------------------------
      !
      !
!*****************************************************************************
End Program wanproj
!*****************************************************************************
