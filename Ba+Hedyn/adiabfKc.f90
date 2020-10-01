subroutine DIMpSOadiabfc(eigenSO,ndimeigenSO,splim, ifail, linit)
!!!

!!! >>>>> NOT DEBUGGED
!!! >>>>> NOT COMPLETE
!!! >>>>> CHECK for orthogonalization in the case of degeneracy between 2 pairs
  
!!! adiabatic following to fix the phase and linear combination of Kramer's degenerate pairs
!!! Choose initially one coefficient to be real positive for the 1st eigenvector of each pair
!!!  and force the phase of the 2nd one of the pair to be such that
!!! if the first of the pair has components (a1,a2,a3,b1,b2,b3)
!!! the 2nd is (-b1*,-b2*,-b3*,a1*,a2*,a3*)
!!! before performing the adiabatic following.
!!! (this guarantees orthogonality of the two transformed vectors of a Kramer's pair
!!!  after the adiabatic following operation)
!!! Adiabatic following:
!!! transform the Kramer's pair vectors (|phi1>, \phi2>) coming from the diagonalization
!!! (after forcing the phase for the 2nd as described above)
!!! as:
!!! |psi1> = <phi1|psi'1> |phi1> + <phi2|psi'1> |phi2> 
!!! |psi2> = <phi1|psi'2> |phi1> + <phi2|psi'2> |phi2>
!!! where (|psi'1>, \psi'2>) are the (final) vectors from the previous step

!!! in addition, follows adiabatically the pairs if a crossing is encountered
!!!  (they come out of diagonalization in increasing energy order)
!!! by assigning them the rank of the eigenvectors in the previous call
!!!  (eigvcpRe,eigvcpIm) with overlap >= splim

!!! ifail flag set to 0 if no problem, set to 1 if some vectors are not assigned

!!! linit = .true. to reinitialize the adiabatic following
!!! lreorder = .true. if re-ordering has taken place in a previous step


!!! >>>>>   eigenE are the eigenenergies coming out of the diagonalization
!!! >>>>> eigenSO are the same ones on input, but they are exchanged if needed
!!! >>>>> in the course of reordering the vectors
!!! >>>>> and eigenE is updated only at the end
  
  implicit none

  logical               :: linit, lreorder   !!! 
  integer               :: ndimeigenSO, ifail
  real (kind=8), dimension(ndimeigenSO) :: eigenSO   ! eigenenergies
  real (kind=8)         :: splim

  integer               :: il, ic, iv, ivp, ivd
  real (kind=8)         :: splim2

  include 'nbrat12.f90'
  include 'units.f90'

  integer               :: ifail
  

!!! real, imaginary part of eigenvectors
  double precision, dimension (ndimorb,ndimorb):: eigvcRe, eigvcIm, eigvcoldRe, eigvcoldIm
  double precision, dimension (ndimorb):: eigenE

!!! (complex) eigenvector (output),
!!! intermediate (complex) eigenvectors (for adiabatic following)
!!!  eigenvector from previous call (input)
  double complex, dimension(ndimorb,ndimorb) :: zeigvc, zeigvctemp, zeigvcold
  real (kind=8), dimension(ndimorb,ndimorb) :: scprodsq
  double precision, dimension(ndimorb) :: worknorm  ! to renormalize eigvec after adiab. follow.
!! eigenvector phase and modulus, and phase of previous eigenvectors
  double precision, dimension (ndimorb,ndimorb) :: evph, evmod
  real (kind=8) :: Remax, evphcorr, evmodmaxc
  integer, dimension (ndimorb) :: ifph   !! which component will be real positive for each eigenvector
  integer :: iat

  integer, dimension(ndimorb) :: indx, indxr, indxp, indxrp
  integer :: nassign

  common/DIM/eigenE,eigvcRe,eigvcIm,eigvcoldRe,eigvcoldIm

  save ifph, zeigvcold
  save indx, indxr, splim2

  ifail = 0

!print *,' check orthogonality before adiabatic following '
!print *
!do ic = 1, ndimorb-1, 2
!   print *,' <',ic,'|',ic,'> = ',dot_product(dcmplx(eigvcRe(:,ic),eigvcIm(:,ic)),dcmplx(eigvcRe(:,ic),eigvcIm(:,ic)))
!   print *,' <',ic+1,'|',ic+1,'> = ',dot_product(dcmplx(eigvcRe(:,ic+1),eigvcIm(:,ic+1)),dcmplx(eigvcRe(:,ic+1),eigvcIm(:,ic+1)))
!   print *,' <',ic,'|',ic+1,'> = ',dot_product(dcmplx(eigvcRe(:,ic),eigvcIm(:,ic)),dcmplx(eigvcRe(:,ic+1),eigvcIm(:,ic+1)))
!enddo

! write(63,'(1x,g15.8,/,6(12(1x,g15.8),/))')xyzdist(3,1),((eigvcRe(il,idimp),eigvcIm(il,idimp),il=1,ndimorb),idimp=1,ndimorb)

!!! eigenvectors are determined with an arbitrary phase
!!! modulus of the components of the eigenvectors
  evmod = abs(dcmplx(eigvcRe,eigvcIm))

  if (linit) then   !!! initialization
     splim2 = splim*2.d0
!!! ! array containing the index of the component with largest modulus for 1st eigenvector of each Kramer's pair
     ifph = maxloc(evmod,1)
!!! ! when 2 components have their modulus equal to max, maxloc gives the lowest index:
!!! ! here choose the one with the largest absolute value of the real part
     do ic=1, ndimorb-1, 2
        Remax = dabs(eigvcRe(ifph(ic),ic))
        evmodmaxc = evmod(ifph(ic),ic)
        do il = 1, ndimorb
           if (il.ne.ifph(ic)) then
              if (dabs(evmod(il,ic)-evmodmaxc).lt.1.d-10) then  ! if 2 components have nearly equal modulus:
                 if (dabs(eigvcRe(il,ic)-Remax).gt.1.d-10) then ! choose the one with the largest real part in absolute value
                    ifph(ic) = il
                    Remax = dabs(eigvcRe(il,ic))
                 endif
              endif
           endif
        enddo
     enddo
         
!   print *,' component with largest modulus for each eigenvector '
!   print *,(ifph(ic),evmod(ifph(ic),ic),ic=1,ndimorb)

     evph = 0.d0   ! phases of the components of the eigenvectors
     do ic = 1, ndimorb-1, 2
!   print *,' ic, ifph = ',ic,ifph(ic)
        evphcorr = atan2(eigvcIm(ifph(ic),ic),eigvcRe(ifph(ic),ic))
        evph(ifph(ic),ic) = evphcorr
!   print *,' evphcorr = ',evphcorr
        do il = 1, ndimorb
           evph(il,ic) = atan2(eigvcIm(il,ic),eigvcRe(il,ic))-evphcorr
        enddo
     enddo
     eigvcRe(:,1:ndimorb:2) = evmod(:,1:ndimorb:2) * dcos(evph(:,1:ndimorb:2))
     eigvcIm(:,1:ndimorb:2) = evmod(:,1:ndimorb:2) * dsin(evph(:,1:ndimorb:2))
! COMPONENTS of the 2ND VECTOR FOR EACH KRAMER'S PAIR
     do ic = 1, ndimorb-1, 2
        eigvcRe(1:ndimorb0,ic+1) = -eigvcRe(ndimorb0+1:ndimorb,ic)
        eigvcIm(1:ndimorb0,ic+1) = eigvcIm(ndimorb0+1:ndimorb,ic)
        eigvcRe(ndimorb0+1:ndimorb,ic+1) = eigvcRe(1:ndimorb0,ic)
        eigvcIm(ndimorb0+1:ndimorb,ic+1) = -eigvcIm(1:ndimorb0,ic)
     enddo
     zeigvc = dcmplx(eigvcRe,eigvcIm)
     zeigvcold = zeigvc
     indx = (/(iv, iv=1, ndimorb)/)   ! initial ordering of the eigenvectors
     indxr = indx
     lreorder = .false.   ! set reordering flag
     linit = .false.

     RETURN

  endif
!!! end of initialization !!! -----------------------------------------------
     
!!! enter here for adiabatic following of the eigenvectors from the previous step
     
  indxp = indx   !!! previous assignment
  indxrp = indxr
  indx = 0.d0
  indxr = 0.d0
  nassign = 0
     
     
!!! COMPONENTS of the 2ND VECTOR FOR EACH KRAMER'S PAIR
  ! done at the end, in case degeneracy 4x4 mixes the pairs
!!! to fix the relative phase of the 2 degenerate eigenvectors
!  do ic = 1, ndimorb-1, 2
!     eigvcRe(1:ndimorb0,ic+1) = -eigvcRe(ndimorb0+1:ndimorb,ic)
!     eigvcIm(1:ndimorb0,ic+1) = eigvcIm(ndimorb0+1:ndimorb,ic)
!     eigvcRe(ndimorb0+1:ndimorb,ic+1) = eigvcRe(1:ndimorb0,ic)
!     eigvcIm(ndimorb0+1:ndimorb,ic+1) = -eigvcIm(1:ndimorb0,ic)
!  enddo

  zeigvctemp = dcmplx(eigvcRe,eigvcIm)

!!!  Check ordering of the pairs:
!!! STEP 1 : first assign the obvious: levels with no crossing since the beginning
  do iv = 1, ndimorb
     scprodsq(iv,iv) = abs(dot_product(zeigvcold(:,iv),zeigvctemp(:,iv)))**2
  enddo
  do iv = 1, ndimorb-1, 2
     if (scprodsq(iv,iv)+scprodsq(iv+1,iv+1).gt.pslim2) then
        indx(iv) = iv
        indx(iv+1) = iv+1
        indxr(iv) = iv
        indxr(iv+1) = iv+1
        nassign = nassign+2
        indxp(iv) = indx(iv)
        indxrp(iv) = indxr(iv)
        indxp(iv+1) = indx(iv+1)
        indxrp(iv+1) = indxr(iv+1)
        eigvcpRe(:,iv) = eigvcRe(:,iv)
        eigvcpIm(:,iv) = eigvcIm(:,iv)
        eigvcpRe(:,iv+1) = eigvcRe(:,iv+1)
        eigvcpIm(:,iv+1) = eigvcIm(:,iv+1)
     endif
  enddo
     
  if (nassign.neq.ndimorb) then   !!!
!!! STEP 2:
!!! try same assignment as in previous call for levels which have not been assigned in step 1
     tempE = eigenE
     
     do ivd = 1, ndimorb-1, 2
        if (indx(ivd).eq.0) then  ! level not assigned in step 1
           ivp = indxp(ivd)   ! assignment in previous adiabKc call
           ivp1 = indxp(ivd+1)
           ! calculate squared scalar products with corresponding pair in previouss tep
           scprodsq(ivp,ivd) = abs(dot_product(zeigvcold(:,ivp),zeigvctemp(:,ivd)))**2
           scprodsq(ivp1,ivd) = abs(dot_product(zeigvcold(:,ivp1),zeigvctemp(:,ivd)))**2
           scprodsq(ivp,ivd1) = abs(dot_product(zeigvcold(:,ivp),zeigvctemp(:,ivd1)))**2
           scprodsq(ivp1,ivd1) = abs(dot_product(zeigvcold(:,ivp1),zeigvctemp(:,ivd1)))**2
           if (scprodsq(ivp,ivd)+scprodsq(ivp1,ivd)+scprodsq(ivp,ivd1)+scprodsq(ivp1,ivd1).gt.splim2) then  ! ivd matches iv from previous call
              nassign = nassign + 2
              indx(ivd) = ivp
              indxr(ivp) = ivd
              indx(ivd+1) = ivp1
              indxr(ivp1) = ivd+1
              eigenSO(ivp) = tempE(ivd)
              eigenSO(ivp1) = tempE(ivd+1)
                  
              eigvcRe(1:ndimorb,ivp) = dreal(zeigvctemp(1:ndimorb,ivd))
              eigvcIm(1:ndimorb,ivp) = dimag(zeigvctemp(1:ndimorb,ivd))
              eigvcRe(1:ndimorb,ivp1) = dreal(zeigvctemp(1:ndimorb,ivd+1))
              eigvcIm(1:ndimorb,ivp1) = dimag(zeigvctemp(1:ndimorb,ivd+1))
           endif
        endif
     enddo
  endif
     
  if (nassign.neq.ndimorb)  then     !  need for further re-ordering
!!! STEP 3
!!! some vectors remain unassigned
!!! calculate all the scalar products of the unassigned eigenvectors
     do ivd = 1, ndimorb-1, 2
        if (indx(ivd).ne.0) cycle    !!! already assigned
        ivd1 = ivd+1
        do ivp = 1, ndimorb-1, 2
           if (indxr(ivp).ne.0) cycle   !!! already used for assignment
           ivp1 = ivp+1
           scprodsq(ivp,ivd) = abs(dot_product(zeigvcold(:,ivp)*zeigvc(:,ivd))**2
           scprodsq(ivp1,ivd) = abs(dot_product(zeigvcold(:,ivp1)*zeigvc(:,ivd))**2
           scprodsq(ivp,ivd1) = abs(dot_product(zeigvcold(:,ivp)*zeigvc(:,ivd1))**2
           scprodsq(ivp1,ivd1) = abs(dot_product(zeigvcold(:,ivp1)*zeigvc(:,ivd1))**2
           if (scprodsq(ivp,ivd)+scprodsq(ivp1,ivd)+scprodsq(ivp,ivd1)+scprodsq(ivp1,ivd1).gt.splim2) then  ! ivd matches iv from previous call
              nassign = nassign + 2
              indx(ivd) = ivp
              indxr(ivp) = ivd
              indx(ivd+1) = ivp1
              indxr(ivp1) = ivd+1
              eigenSO(ivp) = tempE(ivd)
              eigenSO(ivp1) = tempE(ivd+1)
                  
              eigvcRe(1:ndimorb,ivp) = dreal(zeigvctemp(1:ndimorb,ivd))
              eigvcIm(1:ndimorb,ivp) = dimag(zeigvctemp(1:ndimorb,ivd))
              eigvcRe(1:ndimorb,ivp1) = dreal(zeigvctemp(1:ndimorb,ivd+1))
              eigvcIm(1:ndimorb,ivp1) = dimag(zeigvctemp(1:ndimorb,ivd+1))
           endif
        enddo
     enddo
  endif

  if (nassign.eq.ndimorb)  then

!!! adiabatic following for Kramer's pairs which are now in correct order:
          
!!! COMPONENTS of the 2ND VECTOR FOR EACH KRAMER'S PAIR
!!! to fix the relative phase of the 2 degenerate eigenvectors
     do ic = 1, ndimorb-1, 2
        eigvcRe(1:ndimorb0,ic+1) = -eigvcRe(ndimorb0+1:ndimorb,ic)
        eigvcIm(1:ndimorb0,ic+1) = eigvcIm(ndimorb0+1:ndimorb,ic)
        eigvcRe(ndimorb0+1:ndimorb,ic+1) = eigvcRe(1:ndimorb0,ic)
        eigvcIm(ndimorb0+1:ndimorb,ic+1) = -eigvcIm(1:ndimorb0,ic)
     enddo

     zeigvctemp = dcmplx(eigvcRe,eigvcIm)

     do ic = 1, ndimorb-1,2
        zeigvc(:,ic) = dot_product(zeigvctemp(:,ic),zeigvcold(:,ic))*zeigvctemp(:,ic) &
           +dot_product(zeigvctemp(:,ic+1),zeigvcold(:,ic))*zeigvctemp(:,ic+1)
        worknorm(ic) = dot_product(zeigvc(:,ic),zeigvc(:,ic))
        zeigvc(:,ic+1) = dot_product(zeigvctemp(:,ic),zeigvcold(:,ic+1))*zeigvctemp(:,ic) &
           +dot_product(zeigvctemp(:,ic+1),zeigvcold(:,ic+1))*zeigvctemp(:,ic+1)
        worknorm(ic+1) = dot_product(zeigvc(:,ic+1),zeigvc(:,ic+1))
     enddo

     worknorm = 1.d0/dsqrt(worknorm)
!!! renormalize (they are "almost" normalized, ie, to 1.d-7 instead of 1.d-14)
   
!print *,' check orthogonality after adiabatic following '
!print *
!do ic = 1, ndimorb-1, 2
!   print *,' <',ic,'|',ic,'> = ',dot_product(zeigvc(1:ndimorb,ic),zeigvc(1:ndimorb,ic))
!   print *,' <',ic+1,'|',ic+1,'> = ',dot_product(zeigvc(1:ndimorb,ic+1),zeigvc(1:ndimorb,ic+1))
!   print *,' <',ic,'|',ic+1,'> = ',dot_product(zeigvc(1:ndimorb,ic),zeigvc(1:ndimorb,ic+1))
!enddo

     do ic = 1, ndimorb
        zeigvc(:,ic) = zeigvc(:,ic)*worknorm(ic)
     enddo

!print *,' check orthogonality after adiabatic following AND renormalization '
!print *
!do ic = 1, ndimorb-1, 2
!   print *,' <',ic,'|',ic,'> = ',dot_product(zeigvc(1:ndimorb,ic),zeigvc(1:ndimorb,ic))
!   print *,' <',ic+1,'|',ic+1,'> = ',dot_product(zeigvc(1:ndimorb,ic+1),zeigvc(1:ndimorb,ic+1))
!   print *,' <',ic,'|',ic+1,'> = ',dot_product(zeigvc(1:ndimorb,ic),zeigvc(1:ndimorb,ic+1))
!enddo

!print *
!print *,' check orthogonality after adiabatic following '
!do idimp = 1, ndimorb
!   print *,' < i  |',idimp,'> = '
!   write(6,'(6(1x,2(1x,d15.8)))') (dot_product(zeigvc(1:ndimorb,ic),zeigvc(1:ndimorb,idimp)),ic=1,ndimorb)
!enddo

!print *,' check d_{jj} (diagonal non-adiabatic coupling) after adiabatic following '
!print *
!print *,' print all non-adiabatic couplings '
!do idimp = 1, ndimorb
!   print *,' < i old |',idimp,'> = '
!   write(6,'(6(1x,2(1x,d15.8)))') (dot_product(zeigvcold(1:ndimorb,ic),zeigvc(1:ndimorb,idimp)),ic=1,ndimorb)
!enddo

!!! save for next call, save for commons
     eigenE = eigenSO
     eigvcRe = real(zeigvc)
     eigvcIm = imag(zeigvc)
     zeigvcold = zeigvc
     eigvcoldRe = real(zeigvcold)
     eigvcoldIm = imag(zeigvcold)

     RETURN

  endif

!!! Here not all the eigenvector pairs have been assigned:
!!! either there is a degeneracy between 2 pairs
!!! or pslim is not judiciously chosen

!!! check for degeneracy between 2 pairs (for Ba+(p) it is not possible for more)
  nbdeg = 2   !!! array containing the number of degeneracies for each level
                        ! (each level is doubly degenerate, Kramer's pairs)
  ideg = 0    !!!
  ideg(1:ndimorb,1) = /(iv,iv=1,ndimorb)/   !!! each vector is degenerate with itself
  ideg(1:ndimorb-1:2,2) = /(2*iv,iv=1,ndimorb0)/ !!! and with its partner from
  ideg(2:ndimorb:2,2) = /(2*iv-1,iv=1,ndimorb0)/ !!! the same Kramer's pair
  do ivd = 1, ndimorb-1, 2
     if (indx(ivd).ne.0) cycle   !!! vector pair already assigned
     do ivp = ivd+2, ndimorb-1, 2
        if (indx(ivp).ne.0) cycle   !!! vector pair already used for assignment
        if (dabs(eigenE(ivd)-eigenE(ivp)).lt.1.d-11) then
           ideg(ivd,nbdeg(ivd)+1) = ivp
           ideg(ivd,nbdeg(ivd)+2) = ivp+1
           ideg(ivd+1,nbdeg(ivd)+1) = ivp
           ideg(ivd+1,nbdeg(ivd)+2) = ivp+1
           nbdeg(ivd) = nbdeg(ivd) + 2
           ideg(ivp,nbdeg(ivp)+1) = ivd
           ideg(ivp,nbdeg(ivp)+2) = ivd+1
           ideg(ivp+1,nbdeg(ivp)+1) = ivd
           ideg(ivp+1,nbdeg(ivp)+2) = ivd+1
           nbdeg(ivp) = nbdeg(ivp) + 2
        endif
     enddo
  enddo

  if (sum(nbdeg).eq.6) then   !!! error return (STOP for the time being)
     print *
     print *,' >>>>> ERROR STOP in adiabfKc <<<<< '
     print *,'   Not all eigenvectors could be assigned, '
     print *,'   but there are no degeneracies: '
     print *,'   try  change pslim or use smaller integration steps '
     print *,' pslim = ',pslim
     STOP    !!! or RETURN but should set a flag to terminate trajectory
  endif
 
!!! treatment of degeneracy 4x4 (2 pairs x 2 pairs)
  ltreated = .false.
  do ivd = 1, ndimorb-1, 2
     if (indx(ivd).ne.0) then   !!! this vector has been assigned, hence it is not part of the degeneracy
        ! fix the relative phase of the 2nd eigenvector
        ivp = indx(ivd)
        ivp1 = indx(ivd+1)
        eigvcRe(1:ndimorb0,ivp1) = -eigvcRe(ndimorb0+1:ndimorb,ivp)
        eigvcIm(1:ndimorb0,ivp1) = eigvcIm(ndimorb0+1:ndimorb,ivp)
        eigvcRe(ndimorb0+1:ndimorb,ivp1) = eigvcRe(1:ndimorb0,ivp)
        eigvcIm(ndimorb0+1:ndimorb,ivp1) = -eigvcIm(1:ndimorb0,ivp)

        zeigvctemp(:,ivp) = dcmplx(eigvcRe(:,ivp),eigvcIm(:,ivp))
        zeigvctemp(:,ivp1) = dcmplx(eigvcRe(:,ivp1),eigvcIm(:,ivp1))
!!! adiabatic following for each vector of the pair:
        zeigvc(:,ivp) = dot_product(zeigvctemp(:,ivp),zeigvcold(:,ivp))*zeigvctemp(:,ivp) &
           +dot_product(zeigvctemp(:,ivp1),zeigvcold(:,ivp))*zeigvctemp(:,ivp1)
        worknorm(ivp) = dot_product(zeigvc(:,ivp),zeigvc(:,ivp))
        zeigvc(:,ivp1) = dot_product(zeigvctemp(:,ivp),zeigvcold(:,ivp1))*zeigvctemp(:,ivp) &
             +dot_product(zeigvctemp(:,ivp1),zeigvcold(:,ivp1))*zeigvctemp(:,ivp1)
        worknorm(ivp1) = dot_product(zeigvc(:,ivp1),zeigvc(:,ivp1))
!!! renormalize (they are "almost" normalized, ie, to 1.d-7 instead of 1.d-14)   
        zeigvc(:,ivp) = zeigvc(:,ivp)/dsqrt(worknorm(ivp))
        zeigvc(:,ivp1) = zeigvc(:,ivp1)/dsqrt(worknorm(ivp1))

     elseif (nbdeg(ivd).gt.2.and.not.Ltreated(ivd)) then   !!! case of degenerate eigenvectors
!!! search the vectors at t-dt which correspond to the degenerate vectors at t
!!! (since there are 3 pairs of levels for a P state,
!!!  one is assigned and the other 2 are in the degeneracy)        
        do iv = 1, nbdeg(ivd)
           zeigvctemp(:,ideg(ivd,iv)) = dcmplx(eigvcRe(:,ideg(ivd,iv)),eigvcIm(:,ideg(ivd,iv)))
        enddo
        ivp = 0
        do iv = 1, ndimorb
           if (indxr(ivp).ne.0) cycle
           ivp = ivp + 1
           idegr(ivd,ivp) = iv
        enddo
        
        do iv = 1, nbdeg(ivd)
           ivr = ideg(ivd,iv)   !!! "real" value of the index for this vector
           zeigvc(:,ivr) = 0.d0
           do ivp = 1, nbdeg(ivd)
              ivpr = idegr(ivd,ivp)    !!! "real" value of the index for this vector
              zeigvc(:,ivr) = zeigvc(:,ivr) + dot_product(zeigvctemp(:,ivpr),zeigvcold(:,ivr))*zeigvctemp(:,ivpr)
           enddo
           worknorm(ivr) = dot_product(zeigvc(:,ivr),zeigvc(:,ivr))
           zeigvc(:,ivr) = zeigvc(:,ivr)/dsqrt(worknorm(ivr))
        enddo

!!! >>>>> NH: TO BE DO (only Franck Sinatra's addicts can understand)
!!! >>>>>>  at this stage, check if orthogonalization (Schmidt) is required <<<<<<
!!! >>>   cf classdegvalp_so from David  <<<
        
!!! mark all the other eigenvectors involved in this degeneracy as treated
        do ivp = 1, nbdeg(ivd)
           ltreated(ideg(ivd,ivp)) = .true.  
        enddo

     endif

  enddo
   ! save for next call, save for commons
  eigenE = eigenSO
  eigvcRe = real(zeigvc)
  eigvcIm = imag(zeigvc)
  zeigvcold = zeigvc
  eigvcoldRe = real(zeigvcold)
  eigvcoldIm = imag(zeigvcold)
   
  !write(53,'(37(1x,g15.8))') xyzdist(3,1), ((evmod(il,idimp),il=1,ndimorb),idimp=1,ndimorb)
  !write(54,'(37(1x,g15.8))') xyzdist(3,1), ((evph(il,idimp),il=1,ndimorb),idimp=1,ndimorb)

  return
end subroutine DIMpSOadiabfc


