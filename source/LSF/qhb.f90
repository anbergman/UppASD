!> @brief
!> Rescales the temperature used in the Langevin heat bath according to quantum statistics
!> @author
!> Lars Bergqvist, Anders Bergman
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module QHB
   !
   use Parameters
   !
   implicit none


   ! Quantum Heat Bath flag
   character(LEN=1) :: do_qhb                  !< Quantum Heat Bath (N/Q)
   character(LEN=2) :: qhb_mode                !< Temperature control of QHB (TM/TR/MT)
   real(dblprec)    :: tcurie

   public

contains

   !> Rescaling of QHB using quasiharmonic approximation
   subroutine qhb_rescale(temperature,temprescale,temprescalegrad,domode,mode,mag)
      use Constants, only : k_bolt_ev
      use AMS, only : magdos, tcmfa, tcrpa, msat, magdos_freq

      implicit none

      real(dblprec),intent(in) :: temperature
      real(dblprec),intent(out) :: temprescale
      real(dblprec),intent(out) :: temprescalegrad
      character(len=1),intent(in) :: domode
      character(len=2),intent(in) :: mode
      real(dblprec),intent(in) :: mag


      real(dblprec) :: emin,emax,chb,qhb,deltae,beta,fac,cfac,tcrit,ecut,dT,Tg
      real(dblprec),dimension(magdos_freq) :: mdos,edos,qint
      real(dblprec),dimension(-2:2) :: Tresc
      integer :: i

      mdos=0.0_dblprec ; edos=0.0_dblprec ; dT=1.0_dblprec
      beta=0.3650_dblprec         !beta critical exponent (3D)
      emin=minval(magdos(:,1))
      emax=maxval(magdos(:,1))
      deltae=(emax-emin)/(magdos_freq-1)
      fac=1.0_dblprec
      if (mode=='TM') then
         tcrit=tcmfa
      elseif(mode=='TR') then
         tcrit=tcrpa
      elseif(mode=='TC') then
         tcrit=tcurie
      endif

      if (do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='P') then
         do i=-2,2
            Tg=temperature+i*dT
            if (Tg <= 0.0_dblprec) Tg=dbl_tolerance
            if (mode=='MT') then
               fac=mag/msat
            else
               if (Tg < tcrit) then
                  fac=(1.0_dblprec-Tg/tcrit)**beta
               else
                  fac=-1.0_dblprec
               endif
            endif
            if (fac>0) then
               chb=k_bolt_ev*Tg*1000
               if(do_qhb=='Q') then
                  mdos=magdos(:,2)/fac
                  edos=magdos(:,1)*fac
                  qint=edos(:)/(exp(edos(:)/chb)-1.0_dblprec)
                  qint(1)=chb
                  qhb=sum(qint(:)*mdos(:))*deltae*fac
                  Tresc(i)=qhb/chb
               elseif(do_qhb=='R') then
                  ! Enforce classical statistics at Tdebye with additional renormalization
                  qint=magdos(:,1)/(exp(magdos(:,1)/emax)-1.0_dblprec)
                  qint(1)=emax
                  cfac=(sum(qint(:)*magdos(:,2))*deltae)/emax
                  mdos=magdos(:,2)/fac/cfac
                  !--------
                  edos=magdos(:,1)*fac*cfac
                  qint=edos(:)/(exp(edos(:)/chb)-1.0_dblprec)
                  qint(1)=chb
                  qhb=sum(qint(:)*mdos(:))*deltae*fac*cfac
                  Tresc(i)=qhb/chb
                  if (Tresc(i) > 1.0_dblprec) Tresc(i)=1.0_dblprec
               else
                  ! Enforce classical statistics at Tdebye with additional renormalization
                  ecut=k_bolt_ev*tcrit*1000.0_dblprec
                  qint=magdos(:,1)/(exp(magdos(:,1)/ecut)-1.0_dblprec)
                  qint(1)=ecut
                  cfac=(sum(qint(:)*magdos(:,2))*deltae)/ecut
!               mdos=magdos(:,2)/fac/cfac
                  mdos=magdos(:,2)/cfac
               !--------
!               edos=magdos(:,1)*fac
                  edos=magdos(:,1)
                  qint=edos(:)/(exp(edos(:)/chb)-1.0_dblprec)
                  qint(1)=chb
                  qhb=sum(qint(:)*mdos(:))*deltae
!               qhb=sum(qint(:)*mdos(:))*deltae*fac
                  Tresc(i)=qhb/chb
                  if (Tresc(i) > 1.0_dblprec) Tresc(i)=1.0_dblprec
               endif
            else
               Tresc(i)=1.0_dblprec
            endif
         enddo
         temprescale=Tresc(0)
         temprescalegrad=max(0.0_dblprec,((1.0_dblprec/12.0_dblprec)*(Tresc(-2)-Tresc(2))+(2.0_dblprec/3.0_dblprec)*(Tresc(1)-Tresc(-1)))/dT)
      else    ! use m-DOS from SQW or BLS, direct integration
         do i=-2,2
            Tg=temperature+i*dT
            if (Tg <= 0.0_dblprec) Tg=dbl_tolerance
            chb=k_bolt_ev*Tg*1000
            qint(:)=magdos(:,1)/(exp(magdos(:,1)/chb)-1.0_dblprec)
            qint(1)=chb
            qhb=sum(qint*magdos(:,2))*deltae
            Tresc(i)=qhb/chb
         enddo
         temprescale=Tresc(0)
         temprescalegrad=max(0.0_dblprec,((1.0_dblprec/12.0_dblprec)*(Tresc(-2)-Tresc(2))+(2.0_dblprec/3.0_dblprec)*(Tresc(1)-Tresc(-1)))/dT)
      endif

      if (mode=='MT') then
      else
         write(*,3001) 'Temperature rescaling from Quantum Heat Bath:',temprescale, 'T_sim: ', temprescale*temperature, 'dx/dt: ',temprescalegrad
      endif
      3001 format (2x,a,f10.4,2x,a,f7.2,2x,a,f8.4)
   end subroutine qhb_rescale

   subroutine init_qhb()
      !
      implicit none
      !
      !Quantum Heat Bath
      do_qhb = 'N'
      qhb_mode= 'TC'
      tcurie=0.01
   end subroutine init_qhb


end module qhb
