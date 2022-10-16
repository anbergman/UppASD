!> Routine for printing input parameters
!> @copyright
!> GNU Public License.
module PrintInput
   use Parameters
   use Profiling
   implicit none
   public

   real(dblprec), dimension(3) :: fdum
   integer, dimension(3) :: idum
   character, dimension(3) :: cdum

contains


   !> Print input parameters
   !! @todo Restructure whole subroutine
   subroutine prninp()
      !
      use InputData
      use prn_averages
      use prn_microwaves
      use prn_trajectories

      !.. Implicit declarations
      implicit none

      !.. Local scalars
      integer :: i, file_id
      character(len=20) :: filn

      !.. Executable statements
      file_id=ofileno

      ! Open outputfile and print starting brackets
      write (filn,'(''inp.'',a,''.json'')') trim(simid)
      open(file_id, file=filn)
      write (file_id,'(a)') '{'

      ! Print header
      call print_header()
      ! Print physical constants 
      !!! call print_constants()

      ! Print general system info
      call print_input1()

      ! Print simulation parameters
      call print_simulation()

      ! Print initial magnetization info
      call print_initmagnetization()

      ! Conditional printing of unusual features
      if(roteul > 0 )    call print_rotation()
      if(demag .ne. 'N') call print_demag()

      ! Print initial phase info
      call print_initialphase()

      ! Print measurement phase info
      call print_measurementphase()

      ! Print correlation info
      call print_correlation_info()

      ! print comment and close file
      call json_key(file_id,'comment')
      write(file_id,'(a)') '"Input data from UppASD simulation"'
      write (file_id,'(a)') "  } "
      close(file_id)

      ! Print unported outputs to old text format
      write (filn,'(''inp.'',a,''.out'')') trim(simid)
      open(file_id, file=filn)
      call print_mwf_info()
      close(file_id)


   CONTAINS


      subroutine print_header()
         call json_key(file_id,'simid')
         call json_string(file_id,simid)
      end subroutine print_header

      !!! subroutine print_constants()
      !!!    use Constants, only : k_bolt, mub
      !!!    write (file_id,'(a)') "# Physical constants"
      !!!    write (file_id,'(a, /20x ,es16.8)') "Boltzman constant:", k_bolt
      !!!    write (file_id,'(a, /20x ,es16.8)') "Bohr magneton:", mub
      !!!    !write (file_id,'(a,es16.8)') "mRy             ", mry
      !!!    write (file_id,'(a)')  " "
      !!! end subroutine print_constants

      subroutine print_input1()
         !
         ! ncell
         call json_key(file_id,'ncell')
         idum(1)=N1;idum(2)=N2;idum(3)=N3
         call json_int(file_id,idum,3)

         ! bc
         call json_key(file_id,'bc')
         cdum(1)=BC1;cdum(2)=BC2;cdum(3)=BC3
         call json_char(file_id,cdum,3)

         ! cell
         call json_key(file_id,'cell')
         write(file_id,'(a)') '['
         !write (file_id,'(a, 20x,3es15.8/20x ,3es15.8,/20x ,3es15.8)') ' "cell" : ',C1, C2, C3
         call json_float(file_id,C1,3,indflag=.true.)
         call json_float(file_id,C2,3,indflag=.true.)
         call json_float(file_id,C3,3,indflag=.true.,contflag=.true.)
         write(file_id,'(25x,a)') '],'

         ! positions
         !!!write (file_id,'(a)') '  "positions" : '
         !!!do i=1,NA
         !!!   write (file_id,'(20x,i0,2x,i0,2x,3es15.8)') anumb_inp(i), atype_inp(i), bas(1:3,i)
         !!!end do

         ! na
         call json_key(file_id,'na')
         call json_int(file_id,(/NA/),1)

         ! nt
         call json_key(file_id,'nt')
         call json_int(file_id,(/NT/),1)

         ! sym
         call json_key(file_id,'sym')
         call json_int(file_id,(/Sym/),1)

         ! natom not used in input
         !write (file_id,'(a, 20x ,i0)') '   "natom" :', Natom

      end subroutine print_input1


      subroutine print_simulation()
         !write (file_id,'(a)') "# Simulation parameters"
         ! llg
         call json_key(file_id,'llg')
         call json_int(file_id,(/llg/),1)

         ! sdealgh
         call json_key(file_id,'sdealgh')
         call json_int(file_id,(/sdealgh/),1)

         ! mensemble
         call json_key(file_id,'mensemble')
         call json_int(file_id,(/mensemble/),1)

      end subroutine print_simulation


      subroutine print_initmagnetization()
         !write (file_id,'(a)') "# Initial magnetization"
         ! initmag
         call json_key(file_id,'initmag')
         call json_int(file_id,(/initmag/),1)
         !!! if(initmag==1) then
         !!!    write (file_id,'(a,i16)') "mseed           ", mseed
         !!! endif
         !!! if(initmag==2) then
         !!!    write (file_id,'(a,i16)') "mseed           ", mseed
         !!!    write (file_id,'(a,es16.8)') "Theta0          ", theta0
         !!!    write (file_id,'(a,es16.8)') "Phi0            ", phi0
         !!! endif
         if(initmag==3) then
            call json_key(file_id,'moments')
            write(file_id,'(a)') '['
            do i=1,NA-1
               call json_float(file_id,aemom_inp(:,i,1,1),3,indflag=.true.)
            enddo
            call json_float(file_id,aemom_inp(:,i,1,1),3,contflag=.true.,indflag=.true.)
            write(file_id,'(25x,a)') '],'
         endif
         if(initmag==4) then
            call json_key(file_id,'restartfile')
            call json_string(file_id,restartfile)
         endif
      end subroutine print_initmagnetization


      subroutine print_rotation()
         call json_key(file_id,"roteul")
         call json_int(file_id,(/roteul/),1)
         call json_key(file_id,"rotang")
         call json_float(file_id,rotang,3)
      end subroutine print_rotation

      !!! subroutine print_rotation()
      !!!    write (file_id,'(a)') "**************** Rotation parameters            ****************"
      !!!    write (file_id,'(a,i16)') "Rotation        ", roteul
      !!!    write (file_id,'(a,3es16.8)') "Rotation Angle  ", rotang(1), rotang(2), rotang(3)
      !!!    write (file_id,*) " "
      !!! end subroutine print_rotation


      subroutine print_demag()
        
         call json_key(file_id,'demag')
         call json_char(file_id,(/demag/),1)
         call json_key(file_id,'demagfield')
         call json_char(file_id,(/demag1,demag2,demag3/),3)
         call json_key(file_id,'demagvol')
         call json_float(file_id,(/demagvol/),1)
      end subroutine print_demag

      !!! subroutine print_demag()
      !!!    write (file_id,'(a)') "**************** Demagnetization field          ****************"
      !!!    write (file_id,'(a,a,a,a,a,a)') "Demag                          ", demag1, '               ', &
      !!!       demag2, '               ', demag3
      !!!    write (file_id,'(a,a)') "Demag                          ", demag
      !!!    write (file_id,'(a,es16.8)') "Volume          ", demagvol
      !!!    write (file_id,*) " "
      !!! end subroutine print_demag


      subroutine print_mwf_info()
         write (file_id,'(a)') "**************** Micro-wave field               ****************"
         if (mwf=='Y'.or.mwf=='P'.or.mwf=='I') then
            write (file_id,'(a,a)') "MWF                            ", mwf
            write (file_id,'(a,3es16.8)') "Direction       ", mwfdir(1), mwfdir(2), mwfdir(3)
            write (file_id,'(a,es16.8)') "Amplitude       ", mwfampl
            write (file_id,'(a,es16.8)') "Frequency       ", mwffreq
            if (mwf=='P') write (file_id,'(a,i8)') "Pulse steps       ", mwf_pulse_time
         else if (mwf=='S'.or.mwf=='W') then
            write (file_id,'(a,a)') "MWF Site dependent             ", mwf
            write (file_id,'(a,3es16.8)') "Direction       ", mwfdir(1), mwfdir(2), mwfdir(3)
            write (file_id,'(a,es16.8)') "Amplitude       ", mwfampl
            write (file_id,'(a,es16.8)') "Frequency       ", mwffreq
            if (mwf=='W') write (file_id,'(a,i8)') "Pulse steps       ", mwf_pulse_time
         endif
         if (mwf_gauss=='Y'.or.mwf_gauss=='P') then
            write (file_id,'(a,a)') "MWF frequency broadened        ", mwf_gauss
            write (file_id,'(a,3es16.8)') "Direction       ", mwf_gauss_dir(1), mwf_gauss_dir(2), mwf_gauss_dir(3)
            write (file_id,'(a,es16.8)') "Amplitude       ", mwf_gauss_ampl
            write (file_id,'(a,es16.8)') "Frequency       ", mwf_gauss_freq
            if (mwf_gauss=='P') write (file_id,'(a,i8)') "Pulse steps       ", mwf_gauss_pulse_time
         else if (mwf_gauss=='S'.or.mwf_gauss=='W') then
            write (file_id,'(a,a)') "MWF freqeuncy broadened Site dependent ", mwf_gauss
            write (file_id,'(a,3es16.8)') "Direction       ", mwf_gauss_dir(1), mwf_gauss_dir(2), mwf_gauss_dir(3)
            write (file_id,'(a,es16.8)') "Amplitude       ", mwf_gauss_ampl
            write (file_id,'(a,es16.8)') "Frequency       ", mwf_gauss_freq
            if (mwf_gauss=='W') write (file_id,'(a,i8)') "Pulse steps       ", mwf_gauss_pulse_time
         endif
         write (file_id,*) " "
      end subroutine print_mwf_info


      subroutine print_initialphase()
         
         call json_key(file_id,'ip_mode')
         call json_string(file_id,ipmode)

         call json_key(file_id,'ip_hfield')
         call json_float(file_id,iphfield,3)

         if(mode .eq. 'S' .or. mode .eq. 'R' .or. mode .eq. 'L') then
            call json_key(file_id,'ip_nphase')
            call json_int(file_id,(/ipnphase/),1)

            if (ipnphase>0) then
               call json_key(file_id,'ip_nstep')
               call json_int(file_id,ipnstep,ipnphase)

               call json_key(file_id,'ip_timestep')
               call json_float(file_id,ipdelta_t,ipnphase)

               call json_key(file_id,'ip_damping')
               call json_float(file_id,iplambda1,ipnphase)
            end if

         else if (mode .eq. 'M' .or. mode .eq. 'H' .or. mode .eq. 'I') then
            call json_key(file_id,'ip_mcanneal')
            call json_int(file_id,(/ipmcnphase/),1)

            if (ipmcnphase>0) then
               call json_key(file_id,'ip_mcnstep')
               call json_int(file_id,ipmcnstep,ipmcnphase)

               call json_key(file_id,'ip_temp')
               call json_float(file_id,ipTemp,ipmcnphase)
            end if
         end if

      end subroutine print_initialphase

      !!! subroutine print_initialphase()
      !!!    write (file_id,'(a)') "**************** Initial phase                  ****************"
      !!!    write (file_id,'(a,a)') "MC, SD, EM                       ", ipmode
      !!!    write (file_id,'(a,3es16.8)') "Field           ", iphfield(1), iphfield(2), iphfield(3)
      !!!    if(ipmode.eq.'M' .or. ipmode.eq.'H') then
      !!!       write (file_id,'(a,i16)') "nphase          ", ipmcnphase
      !!!       write (file_id,'(a)') "MC steps, T :"
      !!!       do i=1,ipmcnphase
      !!!          write (file_id,'(i16,es16.8)') ipmcnstep(i), ipTemp(i)
      !!!       enddo
      !!!    elseif(ipmode.eq.'S') then
      !!!       write (file_id,'(a,i16)') "nphase          ", ipnphase
      !!!       write(file_id,'(a)') "Nstep         dt            T        lambda1      lambda2 "
      !!!       if (ipnphase>0) then
      !!!          do i=1,ipnphase
      !!!             write(file_id,'(i10,5es12.5)') ipnstep(i),ipdelta_t(i), iptemp(i), iplambda1(i), iplambda2(i)
      !!!          enddo
      !!!       endif
      !!!    endif
      !!!    write (file_id,'(a)') " "
      !!! end subroutine print_initialphase


      subroutine print_measurementphase()

         call json_key(file_id,'mode')
         call json_char(file_id,(/mode/),1)

         call json_key(file_id,'temp')
         call json_float(file_id,(/temp/),1)

         call json_key(file_id,'hfield')
         call json_float(file_id,hfield,3)

         if(mode .eq. 'S' .or. mode .eq. 'R' .or. mode .eq. 'L') then
            call json_key(file_id,'nstep')
            call json_int(file_id,(/nstep/),1)

            call json_key(file_id,'timestep')
            call json_float(file_id,(/delta_t/),1)

            call json_key(file_id,'damping')
            call json_float(file_id,(/mplambda1/),1)
         else if (mode .eq. 'M' .or. mode .eq. 'H' .or. mode .eq. 'I') then
            call json_key(file_id,'mcnstep')
            call json_int(file_id,(/mcnstep/),1)
         end if

         call json_key(file_id,'do_avrg')
         call json_char(file_id,(/do_avrg/),1)
         if(do_avrg .eq. 'Y') then
            call json_key(file_id,'avrg_step')
            call json_int(file_id,(/avrg_step/),1)
         end if

         call json_key(file_id,'do_tottraj')
         call json_char(file_id,(/do_tottraj/),1)
         if(do_tottraj .eq. 'Y') then
            call json_key(file_id,'tottraj_step')
            call json_int(file_id,(/tottraj_step/),1)
         end if

         call json_key(file_id,'do_cumu')
         call json_char(file_id,(/do_cumu/),1)
         if(do_cumu .eq. 'Y') then
            call json_key(file_id,'cumu_step')
            call json_int(file_id,(/cumu_step/),1)
         end if

         call json_key(file_id,'plotenergy')
         call json_int(file_id,(/plotenergy/),1)
         !if(plotenergy>0) then
         !   call json_key(file_id,'ene_step')
         !   call json_int(file_id,(/ene_step/),1)
         !end if

         call json_key(file_id,'do_spintemp')
         call json_char(file_id,(/do_spintemp/),1)
         if(do_spintemp .eq. 'Y') then
            call json_key(file_id,'spintemp_step')
            call json_int(file_id,(/spintemp_step/),1)
         end if

         call json_key(file_id,'ntraj')
         call json_int(file_id,(/ntraj/),1)
         if(ntraj>0) then
            call json_key(file_id,'traj_atom')
            call json_int(file_id,traj_atom,ntraj)
            call json_key(file_id,'traj_step')
            call json_int(file_id,traj_step,ntraj)
         end if
      end subroutine print_measurementphase

      subroutine print_correlation_info()

         use Correlation
         use Correlation_utils
         use Omegas

         call json_key(file_id,'do_sc')
         call json_string(file_id,do_sc)

         if (do_sc .ne. 'N') then

            if (do_sc .eq. 'Y' .or. do_sc .eq. 'C') then
               call json_key(file_id,'sc_sep')
               call json_int(file_id,(/sc%sc_sep/),1)
            end if

            if (do_sc .eq. 'Y' .or. do_sc .eq. 'Q') then
               call json_key(file_id,'sc_step')
               call json_int(file_id,(/sc%sc_step/),1)

               call json_key(file_id,'sc_nstep')
               call json_int(file_id,(/sc%sc_nstep/),1)
            end if

            call json_key(file_id,'do_sc_local_axis')
            call json_string(file_id,sc%do_sc_local_axis)

            call json_key(file_id,'sc_local_axis_mix')
            call json_float(file_id,(/sc%sc_local_axis_mix/),1)

            call json_key(file_id,'sc_window_fun')
            call json_int(file_id,(/sc_window_fun/),1)

            call json_key(file_id,'sc_average')
            call json_string(file_id,sc%sc_average)

            call json_key(file_id,'qpoints')
            call json_string(file_id,qpoints)

            call json_key(file_id,'qfile')
            call json_string(file_id,qfile)

            call json_key(file_id,'do_conv')
            call json_string(file_id,do_conv)

            if (do_conv .ne. 'N') then
               call json_key(file_id,'sigma_q')
               call json_float(file_id,(/sigma_q/),1)

               call json_key(file_id,'sigma_w')
               call json_float(file_id,(/sigma_w/),1)

               call json_key(file_id,'lorentz_q')
               call json_float(file_id,(/LQfactor/),1)

               call json_key(file_id,'lorenrz_w')
               call json_float(file_id,(/LWfactor/),1)

            end if


         end if

      end subroutine print_correlation_info

      !!! subroutine print_measurementphase()
      !!!    use Temperature, only : grad

      !!!    write (file_id,'(a)') "**************** Measurement phase              ****************"
      !!!    write (file_id,'(a,a)') "MC, SD, MEP                         ", mode
      !!!    write (file_id,'(a,a)') "Gradient				 ", grad
      !!!    if(ham_inp%do_ewald.ne.'N')write (file_id,'(a,a)') "Ewald Summation                ", ham_inp%do_ewald
      !!!    if(ham_inp%do_ewald.ne.'N') write(file_id,'(a,es16.8)') "Ewald Parameter                ", ham_inp%Ewald_alpha
      !!!    if(grad.eq.'N') then
      !!!       write (file_id,'(a,es16.8)') "T               ", Temp
      !!!    end if
      !!!    write (file_id,'(a,3es16.8)') "Field           ", hfield(1), hfield(2), hfield(3)
      !!!    if(mode.eq.'M') then
      !!!       write (file_id,'(a,i16)') "Nstep           ", mcnstep
      !!!       write (file_id,'(a,2i16)') "Step, Buffer    ", mcavrg_step, mcavrg_buff
      !!!    elseif(mode.eq.'S') then
      !!!       write (file_id,'(a,2es16.8)') "Lambda          ", mplambda1, mplambda2
      !!!       write (file_id,'(a,i16)') "Nstep           ", Nstep
      !!!       write (file_id,'(a,es16.8)') "delta_t         ", delta_t
      !!!    endif
      !!!    write (file_id,'(a)') "   "
      !!!    write (file_id,'(a,a,2i16)') "Tottraj                        ", do_tottraj, tottraj_step, tottraj_buff
      !!!    write (file_id,'(a,i16)')    "No. trajectories", ntraj
      !!!    do i=1,ntraj
      !!!       write (file_id,'(a,3i16)')   "Atom Trajectory ", traj_atom(i), traj_step(i), traj_buff(i)
      !!!    enddo
      !!! end subroutine print_measurementphase


   end subroutine prninp

         subroutine json_key(fileno,key)
            implicit none
            integer, intent(in) :: fileno
            character(*), intent(in) :: key
            !
            write(fileno,'(a20,a)',advance='no') '"'//trim(key)//'"','  :  '
         end subroutine json_key

         subroutine json_string(fileno,val)
            implicit none
            integer, intent(in) :: fileno
            character(*), intent(in) :: val
            !
            write(fileno,'(1x,a,a,a)') ' "',val,'" ,'
         end subroutine json_string

         subroutine json_char(fileno,val,nel,contflag,indflag)
            implicit none
            integer, intent(in) :: fileno
            integer, intent(in) :: nel
            character, dimension(nel), intent(in) :: val
            logical, intent(in), optional :: contflag
            logical, intent(in), optional :: indflag
            !
            integer :: iel
            logical :: cont 
            logical :: ind 

            !
            cont = .false.
            if (present(contflag)) cont=contflag
            ind = .false.
            if (present(indflag)) ind=indflag
            !
            if(nel==1) then
               write(fileno,'(a,a,a)') ' "',val(1),'" ,'
            else
               if(ind) then
                  write(fileno,'(26x,a)',advance='no') ' [ '
               else
                  write(fileno,'(1x,a)',advance='no') ' [ '
               end if
               do iel=1,nel-1
                  write(fileno,'(a,a,a)',advance='no') '"',val(iel),'" , '
               end do
               if (cont) then 
                  write(fileno,'(a,a,a)') '"',val(iel),'" ] '
               else
                  write(fileno,'(a,a,a)') '"',val(iel),'" ] , '
               end if
            end if
         end subroutine json_char

         subroutine json_float(fileno,val,nel,contflag,indflag)
            implicit none
            integer, intent(in) :: fileno
            integer, intent(in) :: nel
            real(dblprec), dimension(nel), intent(in) :: val
            logical, intent(in), optional :: contflag
            logical, intent(in), optional :: indflag
            !
            integer :: iel
            logical :: cont
            logical :: ind
            !
            cont = .false.
            if(present(contflag)) cont=contflag
            ind = .false.
            if(present(indflag)) ind=indflag

            if(nel==1) then
               write(fileno,'(g14.6,a)') val(1),' ,'
            else
               if(ind) then
                  write(fileno,'(26x,a)',advance='no') ' [ '
               else 
                  write(fileno,'(1x,a)',advance='no') ' [ '
               end if
               do iel=1,nel-1
                  write(fileno,'(g14.6,a)',advance='no') val(iel),' , '
               end do
               if(cont) then
                  write(fileno,'(g14.6,a)') val(iel),' ]  '
               else
                  write(fileno,'(g14.6,a)') val(iel),' ] , '
               end if
            end if
         end subroutine json_float

         subroutine json_int(fileno,val,nel,contflag,indflag)
            implicit none
            integer, intent(in) :: fileno
            integer, intent(in) :: nel
            integer, dimension(nel), intent(in) :: val
            logical, intent(in), optional :: contflag
            logical, intent(in), optional :: indflag
            !
            integer :: iel
            logical :: cont
            logical :: ind

            cont = .false.
            if(present(contflag)) cont=contflag
            ind = .false.
            if(present(indflag)) ind=indflag
            !
            if(nel==1) then
               write(fileno,'(i8,a)') val(1),' ,'
            else
               if(ind) then
                  write(fileno,'(26x,a)',advance='no') ' [ '
               else
                  write(fileno,'(1x,a)',advance='no') ' [ '
               end if
               do iel=1,nel-1
                  write(fileno,'(i8,a)',advance='no') val(iel),' , '
               end do
               if(cont) then
                   write(fileno,'(i8,a)') val(iel),' ]  '
               else
                   write(fileno,'(i8,a)') val(iel),' ] , '
                end if
            end if
         end subroutine json_int

!!!          subroutine json_vector(fileno,val,nel)
!!!             implicit none
!!!             integer, intent(in) :: fileno
!!!             integer, intent(in) :: nel
!!!             real(dblprec), dimension(nel), intent(in) :: val
!!!             !
!!!             integer :: iel
!!!             !
!!!             write(fileno,'(a)',advance='no') ' [ '
!!!             do iel=1,nel-1
!!!                write(fileno,'(f12.6,a)',advance='no') val(iel),' , '
!!!             end do
!!!             write(fileno,'(f12.6,a)') val(iel),' ] '
!!!          end subroutine json_vector

end module PrintInput
