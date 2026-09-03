!------------------------------------------------------------------------------
! NASA/GSFC, Software Systems Support Office, Code 610.3
!------------------------------------------------------------------------------
!
! PROGRAM:  DerivedVars_mod
!
! AUTHOR:
! Eric Kemp, NASA SSSO/SSAI
!
! DESCRIPTION:
! Contains utility routines for calculating derived variables.
!
! REVISION:
! 9 May 2012 - First version
!
!------------------------------------------------------------------------------

module DerivedVars_mod

   ! Import modules
   use FieldWPS_mod
   use WrfRH_mod

   ! Change defaults
   implicit none
   private

   ! Public routines
   public :: calcEdgePressures
   public :: calcLayerPressures
   public :: checkDuplicatesXlvls
   public :: bubbleSortXlvl
   public :: bubbleSort
   public :: copySlabsToArray3d
   public :: copySlabToArray2d
   public :: calcEdgeGeopotentialHeights
   public :: calcLayerGeopotentialHeights

   ! Local constants
   real, parameter :: G_0 = 9.80665 ! Gravity at mean sea level
   real, parameter :: R_d = 287.    ! Dry air gas constant

contains

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calcEdgePressures
   !
   ! DESCRIPTION:  Calculate pressures at edge of GEOS model levels.
   !
   !---------------------------------------------------------------------------

   subroutine calcEdgePressures(longitudeDimension,latitudeDimension, &
        verticalDimension,modelTopPressure,layerPressureThicknesses, &
        edgePressures)

      ! Arguments
      integer,intent(in) :: longitudeDimension
      integer,intent(in) :: latitudeDimension
      integer,intent(in) :: verticalDimension
      real,intent(in) :: modelTopPressure
      real,intent(in) :: layerPressureThicknesses(longitudeDimension, &
           latitudeDimension,verticalDimension)
      real,allocatable,intent(out) :: edgePressures(:,:,:)

      ! Local variables
      integer :: i,j,k

      ! Allocate memory
      if (allocated(edgePressures)) then
         print*,'ERROR, edgePressures already allocated!'
         stop 1
      end if
      ! There are two edges to each layer, so there will be one extra
      ! edge when looping through all the layers
      allocate(edgePressures(longitudeDimension,latitudeDimension, &
           verticalDimension+1))

      ! Calculate edge pressures
      do k = 1, verticalDimension+1
         do j = 1, latitudeDimension
            do i = 1, longitudeDimension
               if (k == 1) then
                  edgePressures(i,j,k) = modelTopPressure
               else
                  edgePressures(i,j,k) = &
                       edgePressures(i,j,k-1) + &
                       layerPressureThicknesses(i,j,k-1)
               end if
            end do
         end do
      end do

      return
   end subroutine calcEdgePressures

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calcLayerPressures
   !
   ! DESCRIPTION:  Calculate pressures at middle of each GEOS layers.
   !
   !---------------------------------------------------------------------------

   subroutine calcLayerPressures(longitudeDimension,latitudeDimension, &
        verticalDimension,edgePressures,layerPressures)

      ! Arguments
      integer,intent(in) :: longitudeDimension
      integer,intent(in) :: latitudeDimension
      integer,intent(in) :: verticalDimension
      real,intent(in) :: edgePressures(longitudeDimension,latitudeDimension,&
           verticalDimension+1)
      real,allocatable,intent(out) :: layerPressures(:,:,:)

      ! Local variables
      integer :: i,j,k

      ! Allocate memory
      if (allocated(layerPressures)) then
         print*,'ERROR, layerPressures already allocated!'
         stop 1
      end if
      allocate(layerPressures(longitudeDimension,latitudeDimension, &
           verticalDimension))
      
      ! Calculate mid-layer pressures
      do k = 1, verticalDimension
         do j = 1, latitudeDimension
            do i = 1, longitudeDimension
               layerPressures(i,j,k) = (edgePressures(i,j,k) + &
                    edgePressures(i,j,k+1)) * 0.5
            end do
         end do
      end do
      
      return
   end subroutine calcLayerPressures

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  checkDuplicatesXlvls
   !
   ! DESCRIPTION:  Loops through FieldWPS variables and checks for duplicate
   ! levels.
   !
   !---------------------------------------------------------------------------

   subroutine checkDuplicatesXlvls(numberOfSlabs,slabs)
      
      ! Arguments
      integer,intent(in) :: numberOfSlabs
      type(FieldWPS), intent(in) :: slabs(numberOfSlabs)

      ! Local variables
      integer :: i,j

      do i = 1, numberOfSlabs-1
         do j = i+1,numberOfSlabs
            if (slabs(i)%xlvl .eq. slabs(j)%xlvl) then
               print*,'ERROR, found slabs with identical levels!'
               print*,'xlvl = ',slabs(i)%xlvl
               stop 1
            end if
         end do
      end do
      return
   end subroutine checkDuplicatesXlvls

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  bubbleSortXlvl
   !
   ! DESCRIPTION:  Sort through array of FieldWPS variables and order by
   ! level.  Based on bubble sort pseudocode on Wikipedia.
   !
   !---------------------------------------------------------------------------

   subroutine bubbleSortXlvl(numberOfSlabs,slabs)

      ! Arguments
      integer,intent(in) :: numberOfSlabs
      type(FieldWPS),intent(inout) :: slabs(numberOfSlabs)

      ! Local variables
      type(FieldWPS) :: field
      logical :: swapped
      integer :: i,n

      n = numberOfSlabs
      do
         swapped=.false.
         do i = 2, n
            if (slabs(i-1)%xlvl > slabs(i)%xlvl) then
               field = copyFieldWPS(slabs(i-1))
               call destroyFieldWPS(slabs(i-1))
               slabs(i-1) = copyFieldWPS(slabs(i))
               call destroyFieldWPS(slabs(i))
               slabs(i) = copyFieldWPS(field)
               call destroyFieldWPS(field)
               swapped=.true.
            end if
         end do
         n = n - 1
         if (.not. swapped) exit
      end do
      return
   end subroutine bubbleSortXlvl

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  bubbleSort
   !
   ! DESCRIPTION:  Sort through real 1d array and order by value.
   ! Based on bubble sort pseudocode on Wikipedia.
   !
   !---------------------------------------------------------------------------

   subroutine bubbleSort(arraySize,array)

      ! Arguments
      integer,intent(in) :: arraySize
      real,intent(inout) :: array(arraySize)

      ! Local variables
      real :: val
      logical :: swapped
      integer :: i,n

      n = arraySize
      do
         swapped=.false.
         do i = 2, n
            if (array(i-1) > array(i)) then
               val = array(i-1)
               array(i-1) = array(i)
               array(i) = val
               swapped=.true.
            end if
         end do
         n = n - 1
         if (.not. swapped) exit
      end do

      return
   end subroutine bubbleSort

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  copySlabsToArray3d
   !
   ! DESCRIPTION:  Copy FieldWPS slabs to a 3D array
   !
   !---------------------------------------------------------------------------

   subroutine copySlabsToArray3d(numberOfSlabs,slabs,array3d)

      ! Arguments
      integer,intent(in) :: numberOfSlabs
      type(FieldWPS),intent(in) :: slabs(numberOfSlabs)
      real,allocatable,intent(out) :: array3d(:,:,:)

      ! Local variables
      integer :: iDim,jDim
      integer :: i,j,k

      iDim = slabs(1)%nx
      jDim = slabs(1)%ny

      allocate(array3d(iDim,jDim,numberOfSlabs))

      do k = 1, numberOfSlabs
         if (iDim .ne. slabs(k)%nx .or. &
             jDim .ne. slabs(k)%ny) then
            print*,'ERROR, dimension mismatch for slabs!'
            print*,'Expected ',iDim,' ',jDim
            print*,'Found ',slabs(k)%nx,' ',slabs(k)%ny
            stop 1
         end if
         do j = 1, jDim
            do i = 1, iDim
               array3d(i,j,k) = slabs(k)%slab(i,j)
            end do
         end do
      end do
      return
   end subroutine copySlabsToArray3d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  copySlabToArray2d
   !
   ! DESCRIPTION:  Copy FieldWPS slab to a 2D array
   !
   !---------------------------------------------------------------------------

   subroutine copySlabToArray2d(slab,array2d)

      ! Arguments
      type(FieldWPS),intent(in) :: slab
      real,allocatable,intent(out) :: array2d(:,:)

      ! Local variables
      integer :: iDim,jDim
      integer :: i,j

      iDim = slab%nx
      jDim = slab%ny

      allocate(array2d(iDim,jDim))

      do j = 1, jDim
         do i = 1, iDim
            array2d(i,j) = slab%slab(i,j)
         end do
      end do

      return
   end subroutine copySlabToArray2d

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calcEdgeGeopotentialHeights
   !
   ! DESCRIPTION:  Calculate geopotential heights at edges of GEOS levels.
   !
   !---------------------------------------------------------------------------
   
   subroutine calcEdgeGeopotentialHeights(longitudeDimension, &
        latitudeDimension,verticalDimension,edgePressures, &
        layerTemperatures,layerSpecificHumidities,soilHeights, &
        edgeGeopotentialHeights)
      
      ! Arguments
      integer,intent(in) :: longitudeDimension
      integer,intent(in) :: latitudeDimension
      integer,intent(in) :: verticalDimension
      real,intent(in) :: edgePressures(longitudeDimension,latitudeDimension, &
           verticalDimension+1)
      real,intent(in) :: layerTemperatures(longitudeDimension, &
           latitudeDimension,verticalDimension)
      real,intent(in) :: layerSpecificHumidities(longitudeDimension, &
           latitudeDimension,verticalDimension)
      real,intent(in) :: soilHeights(longitudeDimension,latitudeDimension)
      real,allocatable,intent(out) :: edgeGeopotentialHeights(:,:,:)

      ! Local variables
      integer :: i,j,k
      real :: t, qv, tv, p1, p2, h1, h2
      real :: inverse_G_0
      
      ! Allocate memory
      allocate(edgeGeopotentialHeights(longitudeDimension,latitudeDimension, &
           verticalDimension+1))
      
      ! Start from ground level
      inverse_G_0 = 1./G_0
      do k = verticalDimension+1, 1, -1
         do j = 1, latitudeDimension
            do i = 1, longitudeDimension               
               if (k == verticalDimension+1) then
                  edgeGeopotentialHeights(i,j,k) = soilHeights(i,j)
               else
                  t = layerTemperatures(i,j,k)
                  qv = layerSpecificHumidities(i,j,k)
                  Tv = calcVirtualTemperature(t, qv)
                  p2 = edgePressures(i,j,k  )
                  p1 = edgePressures(i,j,k+1)
                  h1 = edgeGeopotentialHeights(i,j,k+1)
                  h2 = h1 + R_d*Tv*inverse_G_0*(log(p1) - log(p2))
                  edgeGeopotentialHeights(i,j,k) = h2
               end if                       
            end do
         end do
      end do
      return
   end subroutine calcEdgeGeopotentialHeights

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calcLayerGeopotentialHeights
   !
   ! DESCRIPTION:  Calculate geopotential heights at middle of layers.
   !
   !---------------------------------------------------------------------------

   subroutine calcLayerGeopotentialHeights(longitudeDimension, &
        latitudeDimension,verticalDimension,edgePressures,layerPressures,&
        layerTemperatures,layerSpecificHumidities,edgeGeopotentialHeights, &
        layerGeopotentialHeights)

      ! Arguments
      integer,intent(in) :: longitudeDimension
      integer,intent(in) :: latitudeDimension
      integer,intent(in) :: verticalDimension
      real,intent(in) :: edgePressures(longitudeDimension,latitudeDimension, &
           verticalDimension+1)
      real,intent(in) :: layerPressures(longitudeDimension,latitudeDimension, &
           verticalDimension)
      real,intent(in) :: layerTemperatures(longitudeDimension, &
           latitudeDimension,verticalDimension)
      real,intent(in) :: layerSpecificHumidities(longitudeDimension, &
           latitudeDimension,verticalDimension)
      real,intent(in) :: edgeGeopotentialHeights(longitudeDimension, &
           latitudeDimension,verticalDimension+1)
      real,allocatable,intent(out) :: layerGeopotentialHeights(:,:,:)

      ! Local variables
      integer :: i,j,k
      real :: t, qv, tv, p1, p2, h1, h2
      real :: inverse_G_0

      allocate(layerGeopotentialHeights(longitudeDimension,latitudeDimension,&
           verticalDimension))

      inverse_G_0 = 1./G_0
      do k = verticalDimension, 1, -1
         do j = 1, latitudeDimension
            do i = 1, longitudeDimension
               t  = layerTemperatures(i,j,k)
               qv = layerSpecificHumidities(i,j,k)
               tv = calcVirtualTemperature(t, qv)
               p1 = edgePressures (i,j,k+1)
               p2 = layerPressures(i,j,k  )
               h1 = edgeGeopotentialHeights(i,j,k+1)
               h2 = h1 + R_d*Tv*inverse_G_0*(log(p1) - log(p2))
               layerGeopotentialHeights(i,j,k) = h2
            end do
         end do
      end do

      return
   end subroutine calcLayerGeopotentialHeights

   !---------------------------------------------------------------------------
   !
   ! ROUTINE:  calcVirtualTemperature
   !
   ! DESCRIPTION:  Calculate virtual temperature.
   !
   !---------------------------------------------------------------------------

   function calcVirtualTemperature(temperature,specificHumidity) &
        result(virtualTemperature)

      ! Arguments
      real,intent(in) :: temperature
      real,intent(in) :: specificHumidity
      
      ! Return argument
      real :: virtualTemperature

      ! Local variables
      real :: mixingRatio
      
      ! Calculate mixing ratio from specific humidity
      mixingRatio = calcMixingRatio(specificHumidity)

      ! Now calculate virtual temperature
      virtualTemperature = temperature*(1. + 0.61*mixingRatio)

      return
   end function calcVirtualTemperature

end module DerivedVars_mod
