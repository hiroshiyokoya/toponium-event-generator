      PROGRAM MAIN81
C...Example program to illustrate how to read in LHEF files, e.g., from 
C...a matrix-element generator, into Pythia 6 for further processing 
C...(resonance decays, parton showering, underlying event, and hadronization)
C...     For LHEF, see: hep-ph/0609017
C...     For LHEF+SLHA/BSM, see also: arXiv:0712.3311

C--------------- PREAMBLE: COMMON BLOCK DECLARATIONS ETC -------------
C...All real arithmetic done in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C...The PYTHIA event record:
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C...Pythia parameters
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
C...EXTERNAL statement links PYDATA on most machines.
      CHARACTER*3 chlun
      EXTERNAL PYDATA

C-------------------------- PYTHIA SETUP -----------------------------

C...1) Open LHEF file on unit LUN, and tell Pythia where to find it.
      LUN=88
      OPEN(LUN,FILE='events.lhe')
      WRITE(CHLUN,'(I3)') LUN
      CALL PYGIVE('MSTP(161)='//CHLUN)
      CALL PYGIVE('MSTP(162)='//CHLUN)

C...2) Initialize Pythia for user proces
      CALL PYINIT('USER',' ',' ',0D0)
      
C------------------------- GENERATE EVENTS ---------------------------

C...Initial values for number of events and cumulative charged multiplicity
      IEV=0
      DNSUM=0D0
      DN2SUM=0D0

C...Get next event from file and process it
 100  CALL PYEVNT

C...If event generation failed, quit loop
      IF(MSTI(51).EQ.1) THEN
        GOTO 999
      ENDIF

C...Else count up number of generated events
      IEV=IEV+1

C...Print first few events, both LHEF input and Pythia output.
      IF(IEV.LE.3) THEN
        CALL PYLIST(7)
        CALL PYLIST(2)
      ENDIF

C.../PYJETS/ now contains a fully generated event.
C...Insert user analysis here (or save event to output) 
C...(example: count charged multiplicity)
      CALL PYEDIT(3)
      DNSUM=DNSUM+N
      DN2SUM=DN2SUM+N**2

C...Loop back to look for next event
      GOTO 100

C...Jump point when end-of-file reached (or other problem encountered)
C...Print final statistics.
 999  CALL PYSTAT(1)

      DNAVG=DNSUM/IEV
      DNRMS=SQRT(DN2SUM/IEV)
      SIGMA2=MAX(0D0,DNRMS**2-DNAVG**2)
      WRITE(*,'(1x,A,1x,F6.1,5x,A,1x,F6.1)') 
     &    'Average Charged Multiplicity was <Nch> =', 
     &    DNAVG, '+/-',SQRT(SIGMA2/IEV)
      PRINT*

      END
