!======MODULE CODE FOR maplib======


subroutine Map(AtomMap1, Weight, StartInd1, NAtom2, Pos1, Pos2, &
    & NAtomMap1, NStartInd1, Dim, NAtom1)
    implicit none
    integer, intent(in) :: NAtomMap1, NStartInd1, Dim, NAtom1, NAtom2
    integer, dimension(0:NAtomMap1-1), intent(in) :: AtomMap1
    real(8), dimension(0:NAtomMap1-1), intent(in) :: Weight
    integer, dimension(0:NStartInd1-1), intent(in) :: StartInd1
    real(8), dimension(0:NAtom1-1, 0:Dim-1), intent(in) :: Pos1
    real(8), dimension(0:NAtom2-1, 0:Dim-1), intent(out) :: Pos2
    integer :: i, j, ind1, ind2
    Pos2 = 0.d0
    do i = 0, NStartInd1 - 2
        ind1 = StartInd1(i)
        ind2 = StartInd1(i+1) - 1
        if (ind1 < 0 .or. ind1 > NAtomMap1 - 1 .or. &
          & ind2 < 0 .or. ind2 > NAtomMap1 - 1) cycle
        do j = ind1, ind2
            Pos2(i,:) = Pos2(i,:) + Weight(j) * Pos1(AtomMap1(j), :)
        enddo
    enddo
end subroutine


function modulehash()
    character(len=40) :: modulehash
    modulehash = 'e17683fb3a35fecbadd500f5f1cc8a59d438aaa4'
end function
