module CoreUtils

contains

    !---------------------------------------------------------------------------
    !> @brief   Calculates core indexes
    !> @detail  Calculates core indexes based on the number of snps, requested
    !>          core length and whether to offset cores.  Will try to make
    !>          each core the same length so in practice each core will be
    !>          longer than requested.  Offset controls whether the first core
    !>          is half length (so offsetting all cores).
    !> @date    October 30, 2017
    !> @param[in] nSnp          Number of snps
    !> @param[in] Jump          Requested core length (legacy name!)
    !> @param[in] offset        Whether to offset
    !> @return  Array of core indexes.  CoreIndex(:,1) are start positions,
    !>          CoreIndex(:,2) are end positions.
    !---------------------------------------------------------------------------
    function calculateCores(nSnp, Jump, offset) result(CoreIndex)
        implicit none

        integer, intent(in) :: nSnp, Jump
        logical, intent(in) :: offset
        integer, dimension(:,:), pointer :: CoreIndex

        double precision :: corelength
        integer :: i, nCores, left

        if (jump /= 0) then
            nCores = nSnp / Jump
        else
            nCores = 10
        end if

        ! Catch the case where Jump (i.e. wanted core length) is greater than nsnps
        if (nCores == 0) then
            nCores = 1
        end if

        corelength = nSnp / nCores
        left = nSnp - nCores * corelength

        if (.not. Offset) then
            allocate(CoreIndex(nCores, 2))
            CoreIndex(1, 1) = 1
            if (left /= 0) then
                CoreIndex(1, 2) = 1 + corelength
            else
                CoreIndex(1, 2) = corelength
            end if
        else
            nCores = nCores + 1
            allocate(CoreIndex(nCores, 2))
            CoreIndex(1, 1) = 1
            if (left /= 0) then
                CoreIndex(1, 2) = 1 + floor(dble(corelength) / 2.0)
            else
                CoreIndex(1, 2) = floor(dble(corelength) / 2.0)
            end if
        end if

        do i = 2, nCores
            CoreIndex(i,1) = CoreIndex(i - 1, 2) + 1
            if (i <= left) then
                CoreIndex(i, 2) = CoreIndex(i - 1, 2) + corelength + 1
            else
                CoreIndex(i, 2) = CoreIndex(i - 1, 2) + corelength
            end if
        end do

        if (Offset) then
            CoreIndex(nCores,2) = nSnp
        end if
    end function CalculateCores

    !---------------------------------------------------------------------------
    !> @brief   Calculates core lengths from an array of haplotype libraries
    !> @date    October 30, 2017
    !> @param[in] libraries     Array of Haplotype libraries
    !> @return  Array of core indexes.  CoreIndex(:,1) are start positions,
    !>          CoreIndex(:,2) are end positions.
    !---------------------------------------------------------------------------
    function getCoresFromLibraries(libraries) result(CoreIndex)
        use HaplotypeLibraryModule

        type(HaplotypeLibrary), dimension(:), intent(in) :: libraries
        integer, dimension(:,:), pointer :: CoreIndex

        integer :: i, start

        allocate(CoreIndex(size(libraries),2))

        start = 1
        do i = 1, size(libraries)
            CoreIndex(i,1) = start
            start = start + libraries(i)%nSnps
            CoreIndex(i,2) = start - 1
        end do
    end function getCoresFromLibraries

    !---------------------------------------------------------------------------
    !> @brief   Reads in core lengths from a file
    !> @detail  Reads in core lengths from a file.  First line of file should be
    !>          number of cores.  Next two lines are irrelevant for this
    !>          function.  Next lines should be start and end positions of each
    !>          core.
    !> @date    October 30, 2017
    !> @param[in] file          File to read core lengths from
    !> @return  Array of core indexes.  CoreIndex(:,1) are start positions,
    !>          CoreIndex(:,2) are end positions.
    !---------------------------------------------------------------------------
    function readInCores(file) result(CoreIndex)
        implicit none

        character(len=300), intent(in) :: file
        integer, dimension(:,:), pointer :: CoreIndex

        integer :: nCores, dumI, i

        open (unit = 25, file = trim(file), status = "unknown")
        read (25, *) nCores
        read (25, *) dumI
        read (25, *) dumI
        allocate(CoreIndex(nCores,2))
        do i = 1, nCores
            read (25, *) dumI , CoreIndex(i,1), CoreIndex(i,2)
        end do
        close(25)
    end function readInCores

    !---------------------------------------------------------------------------
    !> @brief   Calculates tail indexes
    !> @detail  Calculates tail indexes based on core indexes and a required
    !>          tail length.
    !> @date    October 30, 2017
    !> @param[in] CoreIndex     Array of core indexes
    !> @param[in] tailLength    The requested tail length
    !> @param[in] nSnp          Number of snps
    !> @return  Array of tail indexes.  TailIndex(:,1) are start positions,
    !>          TailIndex(:,2) are end positions.
    !---------------------------------------------------------------------------
    function calculateTails(CoreIndex, tailLength, nSnp) result(TailIndex)
        implicit none

        integer, dimension(:,:), intent(in) :: CoreIndex
        integer, intent(in) :: tailLength, nSnp
        integer, dimension(:,:), pointer :: TailIndex

        integer :: i, nCores

        nCores = size(CoreIndex,1)
        allocate(TailIndex(nCores,2))

        do i = 1, nCores
            TailIndex(i,1) = max(1,CoreIndex(i,1) - taillength)
            TailIndex(i,2) = min(nSnp,CoreIndex(i,2) + taillength)
        end do
    end function calculateTails

    !---------------------------------------------------------------------------
    !> @brief   DEPRECEATED Calculates tail indexes
    !> @detail  DEPRECEATED Calculates tail indexes based on core indexes and
    !>          a required core and tail length.  calcuateTails should be used
    !>          since it's more general purpose.
    !> @date    October 30, 2017
    !> @param[in] CoreIndex     Array of core indexes
    !> @param[in] nSnp          Number of snps
    !> @param[in] Jump          Core length requested in call to calculateCores
    !> @param[in] tailLength    The requested core and tail length
    !> @return  Array of tail indexes.  TailIndex(:,1) are start positions,
    !>          TailIndex(:,2) are end positions.
    !---------------------------------------------------------------------------
    function OldCalculateTails(CoreIndex, nSnp, Jump, CoreAndTailLength) result(TailIndex)
        implicit none

        integer, dimension(:,:), intent(in) :: CoreIndex
        integer, intent(in) :: nSnp, Jump, CoreAndTailLength
        integer, dimension(:,:), pointer :: TailIndex

        integer :: ltail, rtail, nCores
        integer :: i

        nCores = size(CoreIndex,1)
        allocate(TailIndex(nCores,2))

        ltail = floor(dble(CoreAndTailLength - Jump) / 2.0)
        rtail = ceiling(dble(CoreAndTailLength - Jump) / 2.0)


        do i = 1, nCores
            TailIndex(i,1) = max(1,CoreIndex(i,1) - ltail)
            TailIndex(i,2) = min(nSnp,CoreIndex(i,2) + rtail)
        end do
    end function OldCalculateTails
end module CoreUtils

