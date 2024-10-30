! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_interior_restoring
!
!> \brief MPAS ocean restoring
!> \author Todd Ringler
!> \date   06/08/2015
!> \details
!>  This module contains routines for computing the tracer tendency due to restoring
!
!-----------------------------------------------------------------------

module ocn_tracer_interior_restoring

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_forcing
   use mpas_stream_manager
   use ocn_config
   use ocn_framework_forcing
   use ocn_constants

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_tracer_interior_restoring_compute, &
             ocn_get_interiorRestoringData,                   &
             ocn_tracer_interior_restoring_forcing_write_restart,  &
             ocn_tracer_interior_restoring_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_interior_restoring_compute
!
!> \brief   computes a tracer tendency due to interior restoring
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This routine computes a tracer tendency due to interior restoring
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_interior_restoring_compute(nTracers, nCellsSolve, minLevelCell, maxLevelCell, layerThickness, &
                  tracers, tracersInteriorRestoringRate, tracersInteriorRestoringValue, tracer_tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! one dimensional arrays
      integer, dimension(:), intent(in) :: &
         minLevelCell, maxLevelCell

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers, &
         tracersInteriorRestoringRate, &
         tracersInteriorRestoringValue

      ! scalars
      integer, intent(in) :: nTracers, nCellsSolve

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
        tracer_tend

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, iLevel, iTracer

      err = 0

!tcx
!      write(*,*) 'otirc tcx1a ',minval(tracersInteriorRestoringValue(1,:,:)),maxval(tracersInteriorRestoringValue(1,:,:))
!      write(*,*) 'otirc tcx1b ',minval(tracersInteriorRestoringRate(1,:,:)),maxval(tracersInteriorRestoringRate(1,:,:))
!      write(*,*) 'otirc tcx2a ',minval(tracersInteriorRestoringValue(2,:,:)),maxval(tracersInteriorRestoringValue(2,:,:))
!      write(*,*) 'otirc tcx2b ',minval(tracersInteriorRestoringRate(2,:,:)),maxval(tracersInteriorRestoringRate(2,:,:))
!      write(*,*) 'otirc tcx3a ',minval(tracer_tend(1,:,:)),maxval(tracer_tend(1,:,:))
!      write(*,*) 'otirc tcx3b ',minval(tracer_tend(2,:,:)),maxval(tracer_tend(2,:,:))

      !$omp parallel
      !$omp do schedule(runtime) private(iLevel, iTracer)
      do iCell=1,nCellsSolve
        do iLevel=minLevelCell(iCell),maxLevelCell(iCell)
          do iTracer=1,nTracers
             tracer_tend(iTracer, iLevel, iCell) = tracer_tend(iTracer, iLevel, iCell)  - layerThickness(iLevel,iCell) &
                                                 * ( tracers(iTracer, iLevel, iCell) &
                                                 - tracersInteriorRestoringValue(iTracer, iLevel, iCell) )  &
                                                 * tracersInteriorRestoringRate(iTracer, iLevel, iCell)
          enddo
        enddo
      enddo
      !$omp end do
      !$omp end parallel

!      write(*,*) 'otirc tcx4a ',minval(tracer_tend(1,:,:)),maxval(tracer_tend(1,:,:))
!      write(*,*) 'otirc tcx4b ',minval(tracer_tend(2,:,:)),maxval(tracer_tend(2,:,:))

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_interior_restoring_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_interior_restoring_init
!
!> \brief   Initializes ocean interior restoring
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This routine initializes fields required for tracer interior restoring
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_interior_restoring_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_interior_restoring_init!}}}

!***********************************************************************

!***********************************************************************
!
!  routine get_interiorRestoringData
!
!> \brief   retrieve data needed to compute interior T,S restoring using monthly climatology
!> \author  Mathew Maltrud
!> \date    08/07/23
!> \details
!>  This routine calls mpas_forcing routines to acquire needed interior forcing data and interpolates
!>    between time levels.
!
!-----------------------------------------------------------------------

    subroutine ocn_get_interiorRestoringData( streamManager, &
        domain, &
        simulationClock, &
        firstTimeStep) !{{{

        type (MPAS_streamManager_type), intent(inout) :: streamManager

        type (domain_type) :: domain
        type (MPAS_timeInterval_type) :: timeStepInteriorRestoring
        type (MPAS_clock_type) :: simulationClock

        type (MPAS_Time_Type) :: currTime
        character(len=strKind) :: timeStamp

        logical, intent(in) :: firstTimeStep
        real(kind=RKIND) :: dt
        type (block_type), pointer :: block

        type (dm_info) :: dminfo

        type (mpas_pool_type), pointer :: forcingPool
        type (mpas_pool_type), pointer :: meshPool
        type (mpas_pool_type), pointer :: tracersInteriorRestoringFieldsPool
        type (mpas_pool_type), pointer :: tracersInteriorMonthlyRestoringPool

        real (kind=RKIND), dimension(:,:,:), pointer :: &
           activeTracersInteriorRestoringValue

        real (kind=RKIND), dimension(:,:), pointer :: &
           interiorTemperatureMonthlyClimatologyValue,  &
           interiorSalinityMonthlyClimatologyValue

        integer, dimension(:), pointer :: minLevelCell, maxLevelCell

        integer, pointer :: nCells, nCellsSolve,  &
           indexTemperatureInteriorRestoringValue, &
           indexSalinityInteriorRestoringValue
        integer :: iCell, iLevel

        character(len=strKIND) :: &
           forcingIntervalMonthly,  &
           forcingReferenceTimeMonthly

        integer :: ierr

        ! initialize monthly forcing to be read from file

!maltrud debug
! write(*,*)' tcxw getData ',firstTimeStep

        if (firstTimeStep) then

           currTime = mpas_get_clock_time( simulationClock, MPAS_NOW, ierr)
           call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
           timeStamp = '0000'//trim(timeStamp(5:))

           forcingIntervalMonthly = "0000-01-00_00:00:00"
           forcingReferenceTimeMonthly = "0000-01-15_00:00:00"

           call MPAS_forcing_init_group( forcingGroupHead,  &
                "interiorRestoringMonthlyClimatology", &
                domain, &
                timeStamp, &
                '0000-01-01_00:00:00', &
                '0001-00-00_00:00:00', &
                config_do_restart)

           call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'interiorRestoringMonthlyClimatology', &  !  forcing group name
                'interiorTemperatureMonthlyClimatologyValue', &  !  array name
                'interior_restoring_monthly_data', &  !  stream name
                'tracersInteriorMonthlyRestoring',  & !  pool name
                'interiorTemperatureMonthlyClimatologyValue',  &   !  array name
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

           call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'interiorRestoringMonthlyClimatology', &  !  forcing group name
                'interiorSalinityMonthlyClimatologyValue', &  !  array name
                'interior_restoring_monthly_data', &  !  stream name
                'tracersInteriorMonthlyRestoring',  & !  pool name
                'interiorSalinityMonthlyClimatologyValue',  &   !  array name
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

           call MPAS_forcing_init_field_data( forcingGroupHead, &
                'interiorRestoringMonthlyClimatology', &
                domain % streamManager, &
                config_do_restart, &
                .false.)

           return
        endif  !  first timestep

!maltrud debug
! write(*,*)'tcxw getAgain ',firstTimeStep

        call MPAS_set_timeInterval(timeStepInteriorRestoring, timeString=config_surface_salinity_monthly_restoring_compute_interval)
        call MPAS_get_timeInterval(timeStepInteriorRestoring,dt=dt)

        call MPAS_forcing_get_forcing(forcingGroupHead, &
             'interiorRestoringMonthlyClimatology', streamManager, dt)

      block => domain % blocklist
      do while (associated(block))

        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(forcingPool, 'tracersInteriorRestoringFields',tracersInteriorRestoringFieldsPool)
        call mpas_pool_get_subpool(block % structs, 'tracersInteriorMonthlyRestoring', tracersInteriorMonthlyRestoringPool)

        call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
        call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(tracersInteriorRestoringFieldsPool, 'index_temperatureInteriorRestoringValue',  &
                                                                          indexTemperatureInteriorRestoringValue)
        call mpas_pool_get_dimension(tracersInteriorRestoringFieldsPool, 'index_salinityInteriorRestoringValue',  &
                                                                          indexSalinityInteriorRestoringValue)

        call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

        call MPAS_pool_get_array(tracersInteriorRestoringFieldsPool, 'activeTracersInteriorRestoringValue',   &
                                                                      activeTracersInteriorRestoringValue)

        call MPAS_pool_get_array(tracersInteriorMonthlyRestoringPool, 'interiorTemperatureMonthlyClimatologyValue',   &
                                                                       interiorTemperatureMonthlyClimatologyValue)
        call MPAS_pool_get_array(tracersInteriorMonthlyRestoringPool, 'interiorSalinityMonthlyClimatologyValue',   &
                                                                       interiorSalinityMonthlyClimatologyValue)

        do iCell = 1, nCells
          do iLevel=minLevelCell(iCell),maxLevelCell(iCell)
            activeTracersInteriorRestoringValue(indexTemperatureInteriorRestoringValue, iLevel, iCell) = &
              interiorTemperatureMonthlyClimatologyValue(iLevel, iCell)
            activeTracersInteriorRestoringValue(indexSalinityInteriorRestoringValue, iLevel, iCell) = &
              interiorSalinityMonthlyClimatologyValue(iLevel, iCell)
          enddo
        enddo

        block => block % next
      end do

!tcx
!      write(*,*) 'ogird tcx0 ',indexTemperatureInteriorRestoringValue,indexSalinityInteriorRestoringValue
!      write(*,*) 'ogird tcx1a ',minval(activeTracersInteriorRestoringValue(1,:,:)),maxval(activeTracersInteriorRestoringValue(1,:,:))
!      write(*,*) 'ogird tcx1b ',minval(activeTracersInteriorRestoringRate(1,:,:)),maxval(activeTracersInteriorRestoringRate(1,:,:))
!      write(*,*) 'ogird tcx2a ',minval(activeTracersInteriorRestoringValue(2,:,:)),maxval(activeTracersInteriorRestoringValue(2,:,:))
!      write(*,*) 'ogird tcx2b ',minval(activeTracersInteriorRestoringRate(2,:,:)),maxval(activeTracersInteriorRestoringRate(2,:,:))
!      write(*,*) 'ogird tcx3a ',minval(tracersInteriorRestoringValue(1,:,:)),maxval(tracersInteriorRestoringValue(1,:,:))
!      write(*,*) 'ogird tcx3b ',minval(tracersInteriorRestoringRate(1,:,:)),maxval(tracersInteriorRestoringRate(1,:,:))
!      write(*,*) 'ogird tcx4a ',minval(tracersInteriorRestoringValue(2,:,:)),maxval(tracersInteriorRestoringValue(2,:,:))
!      write(*,*) 'ogird tcx4b ',minval(tracersInteriorRestoringRate(2,:,:)),maxval(tracersInteriorRestoringRate(2,:,:))

    end subroutine ocn_get_interiorRestoringData

!***********************************************************************
!
!  routine ocn_tracer_interior_restoring_forcing_write_restart
!
!> \brief   writes restart timestamp for interior restoring data to be read in on future restart
!> \author  Mathew Maltrud
!> \date    08/08/2023

!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_interior_restoring_forcing_write_restart(domain)!{{{

      type(domain_type) :: domain

      call MPAS_forcing_write_restart_times(forcingGroupHead)

    end subroutine ocn_tracer_interior_restoring_forcing_write_restart!}}}

!***********************************************************************

end module ocn_tracer_interior_restoring

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
