SUBROUTINE SKBI(Pos, BoxL, NAtom, AtomTypes, AtomType_B, AtomType_W, L0, RandIter, G_BB, G_WW, G_BW)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NAtom, AtomType_B, AtomType_W, RandIter
    INTEGER, INTENT(IN), DIMENSION(0:NAtom-1) :: AtomTypes
    REAL(8), INTENT(IN) :: L0
    REAL(8), INTENT(IN), DIMENSION(0:NAtom-1, 0:2) :: Pos
    REAL(8), INTENT(IN), DIMENSION(0:2) :: BoxL
    REAL(8), INTENT(OUT) :: G_BB, G_WW, G_BW
    
    INTEGER :: i, iter
    REAL(8) :: muB, muW, muBB, muWW, muBW, SVol
    REAL(8), DIMENSION(0:RandIter-1) :: NB, NW 
    REAL(8), DIMENSION(0:2) :: rij, invBoxL, Origin
    
    NB = 0.d0
    NW = 0.d0
    invBoxL = MERGE(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    SVol = L0 * L0 * L0
    
    DO iter = 0, RandIter - 1
        ! find origin of cubical subvolume
        Origin = (/ 0.5*L0 + (BoxL(0)-L0) * RAND(0), &
                    0.5*L0 + (BoxL(1)-L0) * RAND(0), &
                    0.5*L0 + (BoxL(2)-L0) * RAND(0)  /)
        
        ! calculate particle numbers within cube
        DO i = 0, NAtom-1
            rij = Pos(i,:) - Origin(:)
            rij = rij - BoxL * DNINT(invBoxL * rij) ! reimage
            IF ( (ABS(rij(0)) <= 0.5*L0) .AND. (ABS(rij(1)) <= 0.5*L0) .AND. (ABS(rij(2)) <= 0.5*L0) ) THEN
                IF (AtomTypes(i) == AtomType_B) NB(iter) = NB(iter) + 1
                IF (AtomTypes(i) == AtomType_W) NW(iter) = NW(iter) + 1
            END IF
        END DO
    END DO
    
    ! calculate KBIs from particle fluctuations
    muB = SUM(NB) / RandIter
    muW = SUM(NW) / RandIter
    muBB = SUM(NB*NB) / RandIter
    muWW = SUM(NW*NW) / RandIter
    muBW = SUM(NB*NW) / RandIter
    G_BB = SVol * ( (muBB - muB*muB)/(muB*muB) - 1.d0/muB )
    G_WW = SVol * ( (muWW - muW*muW)/(muW*muW) - 1.d0/muW )
    G_BW = SVol * ( (muBW - muB*muW)/(muB*muW) )
END SUBROUTINE


SUBROUTINE ZDENSITY(Pos, BoxL, NAtom, rho, NBins)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NAtom, NBins
    REAL(8), INTENT(IN), DIMENSION(0:NAtom-1, 0:2) :: Pos
    REAL(8), INTENT(IN), DIMENSION(0:2) :: BoxL
    REAL(8), INTENT(OUT), DIMENSION(0:NBins-1) :: rho
    
    INTEGER :: i, iter
    REAL(8) :: di, z0, dz
    REAL(8), DIMENSION(0:2) :: invBoxL
    
    dz = (BoxL(2) - 0.d0) / NBins
    invBoxL = MERGE(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    DO iter = 0, NBins-1
        z0 = 0.d0 + dz * (iter + 0.5)
        DO i = 0, NAtom-1
            di = Pos(i,2) - z0
            di = di - BoxL(2) * DNINT(invBoxL(2) * di) ! reimage along Z axis
            IF ( ABS(di) <= 0.5 * dz) rho(iter) = rho(iter) + 1
        END DO
    END DO            
END SUBROUTINE


SUBROUTINE GRIDINSERT(Pos, BoxL, NAtom, InsPos, NIns, Rcav, Ncav)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NAtom, NIns
    REAL(8), INTENT(IN) :: Rcav
    REAL(8), INTENT(IN), DIMENSION(0:2) :: BoxL
    REAL(8), INTENT(IN), DIMENSION(0:NAtom-1, 0:2) :: Pos
    REAL(8), INTENT(IN), DIMENSION(0:NIns-1, 0:2) :: InsPos
    INTEGER, INTENT(OUT) :: Ncav
    
    INTEGER :: i, j
    LOGICAL :: isCav
    REAL(8) :: rsq,  Rcavsq
    REAL(8), DIMENSION(0:2) :: invBoxL, Posi, Posj, rij
    
    !VDW Radii (Bondii et. al.)
    Rcavsq = Rcav * Rcav
    invBoxL = MERGE(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    Ncav = 0.d0
    
    DO i = 0, NIns - 1
        isCav = .TRUE.
        Posi = InsPos(i,:)
        DO j = 0, NAtom - 1
            ! obtain cavity location from supplied grid
            Posj = Pos(j,:)
            rij = Posi - Posj
            rij = rij - BoxL * ANINT(rij * invBoxL) ! minimage
            rsq = SUM(rij * rij)
            IF (rsq < Rcavsq) THEN
                isCav = .FALSE.
                EXIT
            END IF
        END DO
        
        IF (isCav) Ncav = Ncav + 1
    END DO
END SUBROUTINE    
    

SUBROUTINE Smeared_RDF_frame(nr, NBins, Bin_centers, Bin_delta, &
AtomTypes, Pos, NAtom, BoxL, CentAtomtype, NeighAtomType, &
NPoints, SpherePos_cent, SpherePos_neigh, NWorkers)

    IMPLICIT NONE

    include '/opt/mpich2/gnu/include/mpif.h'

    INTEGER, INTENT(IN) :: NAtom, NBins, NPoints, CentAtomType, NeighAtomType, NWorkers
    REAL(8), INTENT(IN) :: Bin_delta
    REAL(8), INTENT(IN), DIMENSION(0:NBins-1) :: Bin_centers
    REAL(8), INTENT(IN), DIMENSION(0:2) :: BoxL
    REAL(8), INTENT(IN), DIMENSION(0:NAtom-1,0:2) :: Pos
    REAL(8), INTENT(IN), DIMENSION(0:NPoints-1,0:2) :: SpherePos_cent, SpherePos_neigh
    INTEGER, INTENT(IN), DIMENSION(0:NAtom-1) :: AtomTypes

    REAL(8), INTENT(OUT), DIMENSION(0:NBins-1) :: nr

    EXTERNAL :: Bin_worker

    INTEGER :: i,j
    LOGICAL :: Applyi, Applyj
    REAL(8) :: Bin_min, Bin_max, Bin_minsq, Bin_maxsq, invBin_delta
    REAL(8), DIMENSION(0:2) :: Posi, Posj, rij, invBoxL

    ! prep worker variables
    REAL(8), DIMENSION(0:NPoints-1, 0:2) :: iPoints, jPoints
    REAL(8), DIMENSION(0:NBins-1) :: worker_nr, master_nr
    INTEGER :: k, StartInd, StopInd, BlockLen, NPointssq

    ! prep mpi variables
    INTEGER :: MASTER, ierr
    PARAMETER (MASTER = 0)
    INTEGER :: worker_id, status(MPI_STATUS_SIZE)

    Bin_min  = MINVAL(Bin_centers) - 0.5 * Bin_delta
    Bin_max =  MAXVAL(Bin_centers) + 0.5 * Bin_delta
    Bin_minsq = Bin_min * Bin_min
    Bin_maxsq = Bin_max * Bin_max
    invBoxL = MERGE(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    invBin_delta = 1.d0/Bin_delta

    nr = 0.0d0
    NPointssq = NPoints*NPoints
    BlockLen = NPointssq / NWorkers
    master_nr = 0.d0
    worker_nr = 0.d0
       
    ! init MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, worker_id, ierr)

    ! outer loop for histogramming
    DO i = 0, NAtom-2
        Posi = Pos(i,:)
        DO j = i+1, NAtom-1
            Posj = Pos(j,:)
            Applyi = ( (CentAtomType == AtomTypes(i)) .AND. (NeighAtomType == AtomTypes(j)) )
            Applyj = ( (CentAtomType == AtomTypes(j)) .AND. (NeighAtomType == AtomTypes(i)) )
            IF (.NOT. (Applyi .OR. Applyj) ) CYCLE

            ! inner loop for Monte carlo sampling over a smeared pair of
            ! molecular centers. This loop is parallelized. 
            ! Prepare data to send to workers
            rij = Posj - Posi
            master_nr = 0.d0
            worker_nr = 0.d0
            IF (Applyi) THEN
                iPoints = SpherePos_cent
                jPoints = SpherePos_neigh
            ELSE IF (Applyj) THEN
                iPoints = SpherePos_neigh
                jPoints = SpherePos_cent
            END IF

            ! parallel loop
            DO k = 0, NWorkers - 1
                ! determine chunk to send to worker
                ! remember fortran loops are inclusive of 
                ! both start and end indices
                StartInd = k * BlockLen
                StopInd = StartInd + BlockLen - 1
                !IF (k == NWorkers - 1) StopInd = NPointssq - 1
                ! call worker (let kth worker handle kth chunk for clarity)
                
                IF (worker_id == k) THEN
                    CALL Bin_worker(worker_nr, rij, BoxL, invBoxL, &
                                    Bin_min, Bin_minsq, Bin_maxsq, invBin_Delta, NBins, &
                                    iPoints, jPoints, StartInd, StopInd, NPoints)
                END IF
                ! sum up bin counts on master
                call MPI_REDUCE(worker_nr, master_nr, 1, MPI_REAL8, MPI_SUM, MASTER, MPI_COMM_WORLD, ierr)
            ENDDO
            ! average and add to main variable nr only on master
            IF (worker_id == MASTER) nr = nr + master_nr / NPointssq
        ENDDO
    ENDDO
    CALL MPI_FINALIZE(ierr)
END SUBROUTINE

   
SUBROUTINE Bin_worker(nr, rij0, BoxL, invBoxL, & 
Bin_min, Bin_minsq, Bin_maxsq, invBin_delta, NBins, & 
iPoints, jPoints, StartInd, StopInd, NPoints)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NBins, NPoints, StartInd, StopInd
    REAL(8), INTENT(IN) :: Bin_min, Bin_minsq, Bin_maxsq, invBin_delta
    REAL(8), DIMENSION(0:2), INTENT(IN) :: rij0, BoxL, invBoxL
    REAL(8), DIMENSION(0:NPoints-1, 0:2), INTENT(IN) :: iPoints, jPoints

    REAL(8), DIMENSION(0:NBins-1), INTENT(OUT) :: nr

    INTEGER :: ii, jj, i_start, i_stop, j_start, j_stop, bin_idx
    REAL(8) :: r, rsq
    REAL(8), DIMENSION(0:2) :: rij

    nr = 0.d0

    ! get start and stop indices
    i_start = INT(StartInd / NPoints)
    j_start = INT(StopInd / NPoints)
    i_stop = StartInd - i_start * NPoints
    j_stop = StopInd - i_stop * NPoints

    ! doubly smeared MC samples from iPoints and jPoints
    DO ii = i_start, i_stop
        DO jj = j_start, j_stop
            ! smear
            rij = rij0 + (jPoints(jj,:) - iPoints(ii,:))
            ! minimage
            rij = rij - BoxL * DNINT(invBoxL * rij)
            ! check if good
            rsq = SUM(rij * rij)
            IF (rsq < Bin_minsq .OR. rsq > Bin_maxsq) CYCLE
            ! bin
            r = SQRT(rsq)
            bin_idx = AINT( (r-Bin_min) * invBin_delta)
            nr(bin_idx) = nr(bin_idx) + 1.0
        ENDDO
    ENDDO

END SUBROUTINE
