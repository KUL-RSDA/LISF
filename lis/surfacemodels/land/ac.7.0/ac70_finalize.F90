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
!
! !ROUTINE: Ac70_finalize
! \label{Ac70_finalize}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   9/4/14: Shugong Wang; initial implementation for Ac70 with LIS-7
!
! !INTERFACE:
subroutine Ac70_finalize(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
    use Ac70_lsmMod
!
! !DESCRIPTION:
!
!  This routine cleans up the allocated memory structures in Ac70
!
!EOP
    implicit none   
    
    integer :: t, n 

    do n=1, LIS_rc%nnest
        ! free memory allocated for each tile
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            deallocate(AC70_struc(n)%ac70(t)%shdfac_monthly)
            deallocate(AC70_struc(n)%ac70(t)%smceq)
            deallocate(AC70_struc(n)%ac70(t)%sstc)
            deallocate(AC70_struc(n)%ac70(t)%sh2o)
            deallocate(AC70_struc(n)%ac70(t)%smc)
            deallocate(AC70_struc(n)%ac70(t)%zss)
            deallocate(AC70_struc(n)%ac70(t)%snowice)
            deallocate(AC70_struc(n)%ac70(t)%snowliq)
        end do  ! tile loop
 
        ! free memory for ac70, the data at tile level
        deallocate(AC70_struc(n)%ac70)

        ! free momory for constant parameter 
        deallocate(AC70_struc(n)%sldpth)
        deallocate(AC70_struc(n)%Thickness)

        ! free momory for initial state variable
        deallocate(AC70_struc(n)%init_stc)
        deallocate(AC70_struc(n)%init_sh2o)
        deallocate(AC70_struc(n)%init_smc)
        !deallocate(AC70_struc(n)%init_zss)
        !deallocate(AC70_struc(n)%init_snowice)
        !deallocate(AC70_struc(n)%init_snowliq)
    end do ! nest loop
  
    deallocate(AC70_struc)
 
end subroutine Ac70_finalize

