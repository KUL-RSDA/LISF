!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp36_snowsimple_update
! \label{noahmp36_snowsimple_update}
!
! !REVISION HISTORY:
!  22 Nov 2022: Isis Brangers; Adjusted update from da_snow
!  This version divides updates over layers and compared to da_snow 
!  is simplified by removing non update related model routines
! !INTERFACE:
subroutine noahmp36_snowsimple_update(n, t, dsneqv, dsnowh)

  use LIS_coreMod
  use NoahMP36_lsmMod
  use module_sf_noahlsm_36 
  use NOAHMP_ROUTINES_36
  
  implicit none
! 
! !DESCRIPTION: 
!  This subroutine updates relevant snow prognostics based
!  on the update to the total SWE (dsneqv) and total
!  snow depth (dsnowh). The updated variables include
!  number of snow layers, snice, snliq, snow temperature 
!  and snow thickness. 
! 
! !ARGUMENTS:
  integer, intent(in)  :: n
  integer, intent(in)  :: t
  real                 :: dsneqv !mm
  real                 :: dsnowh !m
!EOP

  real, allocatable, dimension(:) :: zsoil
  real, allocatable, dimension(:) :: snice
  real, allocatable, dimension(:) :: snliq
  real, allocatable, dimension(:) :: stc
  real, allocatable, dimension(:) :: dzsnso
  real, allocatable, dimension(:) :: zsnso
  real, allocatable, dimension(:) :: sice

  integer :: snl_idx,i,j,iz
  integer :: iloc, jloc
  real    :: sneqv,snowh
  real    :: sneqv1,snowh1
  real    :: ponding1,ponding2
  integer :: newnode
  integer :: isnow, nsoil, nsnow,soiltyp

! local
  real    :: BDSNOW

  isnow = noahmp36_struc(n)%noahmp36(t)%isnow
  nsoil = noahmp36_struc(n)%nsoil
  nsnow = noahmp36_struc(n)%nsnow

  allocate(snice(-nsnow+1:0))
  allocate(snliq(-nsnow+1:0))
  allocate(stc(-nsnow+1:nsoil))
  allocate(dzsnso(-nsnow+1:nsoil))
  allocate(zsnso(-nsnow+1:nsoil))
  allocate(sice(nsoil))
  
  sneqv = noahmp36_struc(n)%noahmp36(t)%sneqv
  snowh = noahmp36_struc(n)%noahmp36(t)%snowh

  zsnso(-nsnow+1:nsoil) = noahmp36_struc(n)%noahmp36(t)%zss(1:nsnow+nsoil) 

! snow/soil layer thickness (m)

  do iz = isnow+1, nsoil
     if(iz == isnow+1) then
        dzsnso(iz) = - zsnso(iz)
     else
        dzsnso(iz) = zsnso(iz-1) - zsnso(iz)
     end if
  end do 
  
  ! set ZSOIL 
  allocate(zsoil(nsoil))
  ! zsoil is negative.
  zsoil(1) = -NOAHMP36_struc(n)%sldpth(1)
  do i = 2, nsoil
     zsoil(i) = zsoil(i-1) - NOAHMP36_struc(n)%sldpth(i)
  enddo


  ! state variables 
  snice(-nsnow+1:0) = &
       NOAHMP36_struc(n)%noahmp36(t)%snowice(1:nsnow)
  snliq(-nsnow+1:0) = &
       NOAHMP36_struc(n)%noahmp36(t)%snowliq(1:nsnow) 
  stc(-nsnow+1:nsoil) = &
       NOAHMP36_struc(n)%noahmp36(t)%sstc(1:nsnow+&
       nsoil) 

! from snowfall routine
  ! no snow layer - update adds snow
  IF(ISNOW == 0.and.(dsneqv.gt.0.and.dsnowh.gt.0))  THEN
     SNOWH = SNOWH + dsnowh
     SNEQV = SNEQV + dsneqv
  END IF
  
  NEWNODE = 0 

  ! no snow layer - after update snowh is sufficient (>0.025) to create first layer
  IF(ISNOW == 0 .AND. SNOWH >= 0.025.and.&
       (dsneqv.gt.0.and.dsnowh.gt.0)) then 
     ISNOW    = -1
     NEWNODE  =  1
     DZSNSO(0)= SNOWH
     STC(0)   = MIN(273.16, NOAHMP36_struc(n)%noahmp36(t)%sfctmp)   ! temporary setup
     SNICE(0) = SNEQV
     SNLIQ(0) = 0.
  END IF
  
  ! snow with layers
  if(isnow.eq.0.and.(dsneqv.lt.0.or.dsnowh.lt.0)) then !no snow, 

  else 
     snowh1 = snowh + dsnowh
     sneqv1 = sneqv + dsneqv
     ! snow layer(s) are present, update adds or removes snow 
     if(isnow.lt.0.and.snowh1.ge.0.and.sneqv1.ge.0.and.newnode==0) then          
        dzsnso(-nsnow+1:(isnow-1)) = 0 
        snice(-nsnow+1:(isnow-1))=0
        ! divide perturbation over layers relative to their size
        do iz=isnow+1,0
           dzsnso(iz) = dzsnso(iz)+dsnowh* dzsnso(iz)/snowh
           snice(iz) = snice(iz)+dsneqv* snice(iz)/sneqv
        enddo 
        SNOWH = SNOWH + dsnowh
        SNEQV = SNEQV + dsneqv
     ! update decreases snow to negative values
     elseif(snowh1.lt.0.or.sneqv1.lt.0) then 
        isnow = 0
        snowh = 0
        sneqv = 0
     endif
    
     sice(:) = max(0.0, NOAHMP36_struc(n)%noahmp36(t)%smc(:)&
          - NOAHMP36_struc(n)%noahmp36(t)%sh2o(:))   
     
     if(isnow < 0) &
          call  combine (nsnow  ,&
          nsoil  ,iloc   ,jloc   ,         & !in
          isnow  ,&
          noahmp36_struc(n)%noahmp36(t)%sh2o   ,&
          stc    ,snice  ,snliq  , & !inout
          dzsnso ,sice   ,snowh  ,sneqv  ,         & !inout
          ponding1       ,ponding2)                  !out
     
     
     if(isnow < 0) &        
          call divide (nsnow  , nsoil  ,   & !in
          isnow  , stc    ,&
          snice  , snliq  , dzsnso )   !inout
     

!set empty snow layers to zero
     do iz = -nsnow+1, isnow
        snice(iz) = 0.
        snliq(iz) = 0.
        stc(iz)   = 0.
        dzsnso(iz)= 0.
        zsnso(iz) = 0.
     enddo
     

! sum up snow mass for layered snow
     IF(ISNOW < 0) THEN  ! MB: only do for multi-layer
        SNEQV = 0.
        SNOWH = 0.      ! Yeosang Yoon
        DO IZ = ISNOW+1,0
           SNEQV = SNEQV + SNICE(IZ) + SNLIQ(IZ)
           SNOWH = SNOWH + DZSNSO(IZ)             ! Yeosang Yoon
        ENDDO
     END IF


! Yeosag Yoon, no snow layer case, limit snow density to 1000
      IF (ISNOW == 0 .AND. SNEQV > 0. .AND. SNOWH > 0.) THEN
           BDSNOW = SNEQV/SNOWH
           IF (BDSNOW >= DENH2O) THEN
               SNOWH  = SNOWH*(BDSNOW/1000.) ! change unit, SNEQV=[mm] SNOWH=[m]
           END IF
      END IF


! Reset ZSNSO and layer thinkness DZSNSO
     do iz = isnow+1, nsoil
        if(iz==isnow+1) then
           zsnso(iz) = -dzsnso(iz)
        else
           zsnso(iz) = zsnso(iz-1)-dzsnso(iz)
        endif
     end do

    
     if (sneqv.le.0 .or. snowh.le.0) then
        sneqv=0
        snowh=0
     endif

     noahmp36_struc(n)%noahmp36(t)%isnow = isnow
     noahmp36_struc(n)%noahmp36(t)%sneqv = sneqv
     noahmp36_struc(n)%noahmp36(t)%snowh = snowh 
     
     NOAHMP36_struc(n)%noahmp36(t)%zss(1:nsnow+&
          nsoil) =  zsnso(-nsnow+1:nsoil)
     NOAHMP36_struc(n)%noahmp36(t)%snowice(1:nsnow) = & 
          snice(-nsnow+1:0) 
     NOAHMP36_struc(n)%noahmp36(t)%snowliq(1:nsnow)  = &        
          snliq(-nsnow+1:0) 
     NOAHMP36_struc(n)%noahmp36(t)%sstc(1:nsnow+&
          nsoil) =  stc(-nsnow+1:nsoil) 
  endif

  deallocate(snice)
  deallocate(snliq)
  deallocate(stc)
  deallocate(dzsnso)
  deallocate(zsnso)
  deallocate(sice)
     
end subroutine noahmp36_snowsimple_update
