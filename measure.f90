
SUBROUTINE RDF(g, Nbins, Bin_centers, Bin_delta, &
AtomTypes, Pos, NAtom, BoxL)
	
	IMPLICIT NONE
	
	INTEGER, INTENT(IN) :: NAtom, Nbins
	REAL(8), INTENT(IN) :: Bin_delta
	REAL(8), INTENT(IN), DIMENSION(0:Nbins-1) :: Bin_centers
	REAL(8), INTENT(IN) :: BoxL
	REAL(8), INTENT(IN), DIMENSION(0:NAtom-1,0:2) :: Pos
	INTEGER, INTENT(IN), DIMENSION(0:NAtom-1) :: AtomTypes
	
	REAL(8), INTENT(INOUT), DIMENSION(0:Nbins-1, 0:2) :: g
	!f2py intent(in,out) :: g
	
	INTEGER :: i,j, Bin_assignment, rdf_type = -1
	REAL(8) :: r, rsq, Bin_min, Bin_max, Bin_minsq, Bin_maxsq, invBoxL, invBin_delta
	REAL(8), DIMENSION(0:2) :: Posi, rij
	
	! Precompute bin extrema for fast comparison
	Bin_min  = MINVAL(Bin_centers) - Bin_delta/2.0
	Bin_max = MAXVAL(Bin_centers) + Bin_delta/2.0
	Bin_minsq = Bin_min * Bin_min
	Bin_maxsq = Bin_max * Bin_max
	invBoxL = 1.d0/BoxL
	invBin_delta = 1.d0/Bin_delta

	DO i = 0, NAtom-2			
		! precompute central atom position for speed
		Posi = Pos(i,:)
		DO j = i+1, NAtom-1
			IF (AtomTypes(i) == 1 .AND. AtomTypes(j) == 1) THEN
				rdf_type = 0 ! benzene-benzene
			ELSE IF (AtomTypes(i) == 2 .AND. AtomTypes(j) == 2) THEN
				rdf_type = 1 ! water-water
			ELSE IF (AtomTypes(i) == 1 .AND. AtomTypes(j) == 2) THEN
				rdf_type = 2 ! benzene-water
			ELSE IF (AtomTypes(i) == 2 .AND. AtomTypes(j) == 1) THEN
				rdf_type = 2 ! water-benzene
			END IF
			
			! Compute pairwise distance
			rij = Pos(j,:) - Posi
			rij(:) = rij(:) - BoxL * DNINT(invBoxL * rij(:)) ! minimum-imaging
			rsq = SUM(rij*rij)
			
			! Bin the data
			IF (rsq < Bin_minsq .OR. rsq > Bin_maxsq) THEN
				CYCLE
			END IF
			
			r = sqrt(rsq)
			Bin_assignment = AINT((r-Bin_min) * invBin_delta)
			g(Bin_assignment, rdf_type) = g(Bin_assignment, rdf_type) + 1.0	
		ENDDO
	ENDDO
END SUBROUTINE		



SUBROUTINE NEARESTNEIGH(fsn_BW, ld_BW, ld_WB,& 
Pos, NAtom, NB, NW, AtomTypes, NCuts, LDUpperCuts_BW, LDUpperCuts_WB, &
FirstShellCut_BW, FirstShellCut_WB, coeff_BW, coeff_WB, BoxL)

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: NAtom, NCuts, NB, NW
	REAL(8), INTENT(IN) :: FirstShellCut_BW, FirstShellCut_WB
	REAL(8), INTENT(IN), DIMENSION(0:NAtom-1,0:2) :: Pos
	INTEGER, INTENT(IN), DIMENSION(0:NAtom-1) :: AtomTypes
	
	REAL(8), INTENT(IN) :: BoxL
	REAL(8), INTENT(IN), DIMENSION(0:NCuts-1,0:3) :: coeff_BW, coeff_WB
	REAL(8), INTENT(IN), DIMENSION(0:NCuts-1) :: LDUpperCuts_BW, LDUpperCuts_WB
	
	REAL(8), INTENT(INOUT), DIMENSION(0:NB - 1, 0:NCuts - 1) :: ld_BW
	REAL(8), INTENT(INOUT), DIMENSION(0:NW - 1, 0:NCuts - 1) :: ld_WB
	REAL(8), INTENT(INOUT), DIMENSION(0:NB - 1) :: fsn_BW
	REAL(8), INTENT(INOUT), DIMENSION(0:NW - 1) :: fsn_WB
	!f2py intent(in,out) :: fsn_BW, ld_BW, fsn_WB, ld_WB

	INTEGER :: i,j,k
	REAL(8) :: rsq, phi, FirstShellCutsq_BW, FirstShellCutsq_WB, invBoxL
	REAL(8), DIMENSION(0:2) :: Posi, rij 
	REAL(8) :: c0, c2, c4, c6, R1sq, R2sq
	REAL(8), DIMENSION(0:NCuts-1) :: LDUpperCutsq_BW, LDLowerCutsq_BW, LDUpperCutsq_WB, LDLowerCutsq_WB
	
	! Precompute for speed
	FirstShellCutsq_BW = FirstShellCut_BW * FirstShellCut_BW
	FirstShellCutsq_WB = FirstShellCut_WB * FirstShellCut_WB
	invBoxL = 1.0/BoxL
	
	DO i = 0, NCuts-1
		LDUpperCutsq_BW(i) = LDUpperCuts_BW(i) * LDUpperCuts_BW(i)
		LDLowerCutsq_BW(i) = 0.8 * 0.8 * LDUpperCutsq_BW(i)
		LDUpperCutsq_WB(i) = LDUpperCuts_WB(i) * LDUpperCuts_WB(i)
		LDLowerCutsq_WB(i) = 0.8 * 0.8 * LDUpperCutsq_WB(i)	 
	ENDDO
	
	DO i = 0, NAtom-2
		! Precompute central atom position for speed
		Posi = Pos(i,:)	
		DO j = i+1, NAtom -1
			! compute pair-distances
			rij = Pos(j,:) - Posi
			rij(:) = rij(:) - BoxL * DNINT(invBoxL * rij(:)) ! minimum-imaging
			rsq = SUM(rij*rij)				
			
			!! no need to do (j,i) typechecking since all the benzene sites are in the beginning
			IF (AtomTypes(i) == 1 .AND. AtomTypes(j) == 2) THEN ! benzene-water
				
	            ! compute first shell neighbors
				IF (rsq <= FirstShellCutsq_BW) THEN
					fsn_BW(i) = fsn_BW(i) + 1.00
				END IF			
				IF (rsq <= FirstShellCutsq_WB) THEN
					fsn_WB(j-NB) = fsn_WB(j-NB) + 1.00
				END IF
			
			    ! loop over each cutoff
				DO k = 0, NCuts - 1
					
					! BW local density
					R1sq = LDLowerCutsq_BW(k)
					R2sq = LDUpperCutsq_BW(k)
					c0 = coeff_BW(k,0)
					c2 = coeff_BW(k,1)
					c4 = coeff_BW(k,2)
					c6 = coeff_BW(k,3)
					IF (rsq >= R2sq) THEN
						phi = 0
					ELSE IF (rsq <= R1sq) THEN
						phi = 1
					ELSE
						phi = c0 + rsq*(c2 + rsq*(c4 + rsq*c6))
					ENDIF
					ld_BW(i,k) = ld_BW(i,k) + phi
					ld_BW(j,k) = ld_BW(j,k) + phi
				
				    ! WB local density
				    R1sq = LDLowerCutsq_WB(k)
					R2sq = LDUpperCutsq_WB(k)
					c0 = coeff_WB(k,0)
					c2 = coeff_WB(k,1)
					c4 = coeff_WB(k,2)
					c6 = coeff_WB(k,3)
					IF (rsq >= R2sq) THEN
						phi = 0
					ELSE IF (rsq <= R1sq) THEN
						phi = 1
					ELSE
						phi = c0 + rsq*(c2 + rsq*(c4 + rsq*c6))
					ENDIF
					ld_WB(i-NB,k) = ld_WB(i-NB,k) + phi
					ld_WB(j-NB,k) = ld_WB(j-NB,k) + phi
				ENDDO
			ENDIF
		ENDDO
	ENDDO
END SUBROUTINE	
