*DECK SPLPMU
      SUBROUTINE SPLPMU (MRELAS, NVARS, LMX, LBM, NREDC, INFO, IENTER,
     +   ILEAVE, IOPT, NPP, JSTRT, IBASIS, IMAT, IBRC, IPR, IWR, IND,
     +   IBB, ANORM, EPS, UU, GG, RPRNRM, ERDNRM, DULNRM, THETA, COSTSC,
     +   XLAMDA, RHSNRM, AMAT, BASMAT, CSC, WR, RPRIM, WW, BU, BL, RHS,
     +   ERD, ERP, RZ, RG, COLNRM, COSTS, PRIMAL, DUALS, SINGLR, REDBAS,
     +   ZEROLV, STPEDG)
C
C    This is a hacked version of the slatec program
C    All ASSIGN statements have been replaced
C
C***BEGIN PROLOGUE  SPLPMU
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SPLP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SPLPMU-S, DPLPMU-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
C     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
C
C     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
C     /REAL (12 BLANKS)/DOUBLE PRECISION/,
C     /SASUM/DASUM/,/SCOPY/DCOPY/,/SDOT/DDOT/,
C     /.E0/.D0/
C
C     THIS SUBPROGRAM IS FROM THE SPLP( ) PACKAGE.  IT PERFORMS THE
C     TASKS OF UPDATING THE PRIMAL SOLUTION, EDGE WEIGHTS, REDUCED
C     COSTS, AND MATRIX DECOMPOSITION.
C     IT IS THE MAIN PART OF THE PROCEDURE (MAKE MOVE AND UPDATE).
C
C     REVISED 821122-1100
C     REVISED YYMMDD
C
C***SEE ALSO  SPLP
C***ROUTINES CALLED  IPLOC, LA05BS, LA05CS, PNNZRS, PRWPGE, SASUM,
C                    SCOPY, SDOT, SPLPDM, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   811215  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890605  Removed unreferenced labels.  (WRB)
C   890606  Removed unused COMMON block LA05DS.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  SPLPMU
      INTEGER IBASIS(*),IMAT(*),IBRC(LBM,2),IPR(*),IWR(*),IND(*),IBB(*)
      REAL             AIJ,ALPHA,ANORM,COSTSC,ERDNRM,DULNRM,EPS,GAMMA,
     * GG,GQ,ONE,RPRNRM,RZJ,SCALR,THETA,TWO,UU,WP,XLAMDA,RHSNRM,
     * ZERO,AMAT(*),BASMAT(*),CSC(*),WR(*),RPRIM(*),WW(*),BU(*),BL(*),
     * RHS(*),ERD(*),ERP(*),RZ(*),RG(*),COSTS(*),PRIMAL(*),DUALS(*),
     * COLNRM(*),RCOST,SASUM,SDOT
      LOGICAL SINGLR,REDBAS,PAGEPL,TRANS,ZEROLV,STPEDG
     integer NPR001
     integer NPR003
C
C***FIRST EXECUTABLE STATEMENT  SPLPMU
      ZERO=0.E0
      ONE=1.E0
      TWO=2.E0
      LPG=LMX-(NVARS+4)
C
C     UPDATE THE PRIMAL SOLUTION WITH A MULTIPLE OF THE SEARCH
C     DIRECTION.
      I=1
      N20002=MRELAS
      GO TO 20003
20002 I=I+1
20003 IF ((N20002-I).LT.0) GO TO 20004
      RPRIM(I)=RPRIM(I)-THETA*WW(I)
      GO TO 20002
C
C     IF EJECTED VARIABLE IS LEAVING AT AN UPPER BOUND,  THEN
C     TRANSLATE RIGHT HAND SIDE.
20004 IF (.NOT.(ILEAVE.LT.0)) GO TO 20006
      IBAS=IBASIS(ABS(ILEAVE))
      SCALR=RPRIM(ABS(ILEAVE))
      NPR001=20009
      GO TO 30001
20009 IBB(IBAS)=ABS(IBB(IBAS))+1
C
C     IF ENTERING VARIABLE IS RESTRICTED TO ITS UPPER BOUND, TRANSLATE
C     RIGHT HAND SIDE.  IF THE VARIABLE DECREASED FROM ITS UPPER
C     BOUND, A SIGN CHANGE IS REQUIRED IN THE TRANSLATION.
20006 IF (.NOT.(IENTER.EQ.ILEAVE)) GO TO 20010
      IBAS=IBASIS(IENTER)
      SCALR=THETA
      IF (MOD(IBB(IBAS),2).EQ.0) SCALR=-SCALR
      NPR001=20013 
      GO TO 30001
20013 IBB(IBAS)=IBB(IBAS)+1
      GO TO 20011
20010 IBAS=IBASIS(IENTER)
C
C     IF ENTERING VARIABLE IS DECREASING FROM ITS UPPER BOUND,
C     COMPLEMENT ITS PRIMAL VALUE.
      IF (.NOT.(IND(IBAS).EQ.3.AND.MOD(IBB(IBAS),2).EQ.0)) GO TO 20014
      SCALR=-(BU(IBAS)-BL(IBAS))
      IF (IBAS.LE.NVARS) SCALR=SCALR/CSC(IBAS)
      NPR001=20017
      GO TO 30001
20017 THETA=-SCALR-THETA
      IBB(IBAS)=IBB(IBAS)+1
20014 CONTINUE
      RPRIM(ABS(ILEAVE))=THETA
      IBB(IBAS)=-ABS(IBB(IBAS))
      I=IBASIS(ABS(ILEAVE))
      IBB(I)=ABS(IBB(I))
      IF(PRIMAL(ABS(ILEAVE)+NVARS).GT.ZERO) IBB(I)=IBB(I)+1
C
C     INTERCHANGE COLUMN POINTERS TO NOTE EXCHANGE OF COLUMNS.
20011 IBAS=IBASIS(IENTER)
      IBASIS(IENTER)=IBASIS(ABS(ILEAVE))
      IBASIS(ABS(ILEAVE))=IBAS
C
C     IF VARIABLE WAS EXCHANGED AT A ZERO LEVEL, MARK IT SO THAT
C     IT CAN'T BE BROUGHT BACK IN.  THIS IS TO HELP PREVENT CYCLING.
      IF(ZEROLV) IBASIS(IENTER)=-ABS(IBASIS(IENTER))
      RPRNRM=MAX(RPRNRM,SASUM(MRELAS,RPRIM,1))
      K=1
      N20018=MRELAS
      GO TO 20019
20018 K=K+1
20019 IF ((N20018-K).LT.0) GO TO 20020
C
C     SEE IF VARIABLES THAT WERE CLASSIFIED AS INFEASIBLE HAVE NOW
C     BECOME FEASIBLE.  THIS MAY REQUIRED TRANSLATING UPPER BOUNDED
C     VARIABLES.
      IF (.NOT.(PRIMAL(K+NVARS).NE.ZERO .AND.
     *          ABS(RPRIM(K)).LE.RPRNRM*ERP(K))) GO TO 20022
      IF (.NOT.(PRIMAL(K+NVARS).GT.ZERO)) GO TO 20025
      IBAS=IBASIS(K)
      SCALR=-(BU(IBAS)-BL(IBAS))
      IF(IBAS.LE.NVARS)SCALR=SCALR/CSC(IBAS)
      NPR001=20028
      GO TO 30001
20028 RPRIM(K)=-SCALR
      RPRNRM=RPRNRM-SCALR
20025 PRIMAL(K+NVARS)=ZERO
20022 CONTINUE
      GO TO 20018
C
C     UPDATE REDUCED COSTS, EDGE WEIGHTS, AND MATRIX DECOMPOSITION.
20020 IF (.NOT.(IENTER.NE.ILEAVE)) GO TO 20029
C
C     THE INCOMING VARIABLE IS ALWAYS CLASSIFIED AS FEASIBLE.
      PRIMAL(ABS(ILEAVE)+NVARS)=ZERO
C
      WP=WW(ABS(ILEAVE))
      GQ=SDOT(MRELAS,WW,1,WW,1)+ONE
C
C     COMPUTE INVERSE (TRANSPOSE) TIMES SEARCH DIRECTION.
      TRANS=.TRUE.
      CALL LA05BS(BASMAT,IBRC,LBM,MRELAS,IPR,IWR,WR,GG,WW,TRANS)
C
C     UPDATE THE MATRIX DECOMPOSITION.  COL. ABS(ILEAVE) IS LEAVING.
C     THE ARRAY DUALS(*) CONTAINS INTERMEDIATE RESULTS FOR THE
C     INCOMING COLUMN.
      CALL LA05CS(BASMAT,IBRC,LBM,MRELAS,IPR,IWR,DUALS,GG,UU,
     * ABS(ILEAVE))
      REDBAS=.FALSE.
      IF (.NOT.(GG.LT.ZERO)) GO TO 20032
C
C     REDECOMPOSE BASIS MATRIX WHEN AN ERROR RETURN FROM
C     LA05CS( ) IS NOTED.  THIS WILL PROBABLY BE DUE TO
C     SPACE BEING EXHAUSTED, GG=-7.
      CALL SPLPDM(
     *MRELAS,NVARS,LMX,LBM,NREDC,INFO,IOPT,
     *IBASIS,IMAT,IBRC,IPR,IWR,IND,IBB,
     *ANORM,EPS,UU,GG,
     *AMAT,BASMAT,CSC,WR,
     *SINGLR,REDBAS)
      IF (.NOT.(SINGLR)) GO TO 20035
      NERR=26
      CALL XERMSG ('SLATEC', 'SPLPMU',
     +   'IN SPLP, MOVED TO A SINGULAR POINT.  THIS SHOULD NOT HAPPEN.',
     +   NERR, IOPT)
      INFO=-NERR
      RETURN
20035 CONTINUE
      GO TO 30002
20038 CONTINUE
20032 CONTINUE
C
C     IF STEEPEST EDGE PRICING IS USED, UPDATE REDUCED COSTS
C     AND EDGE WEIGHTS.
      IF (.NOT.(STPEDG)) GO TO 20039
C
C     COMPUTE COL. ABS(ILEAVE) OF THE NEW INVERSE (TRANSPOSE) MATRIX
C     HERE ABS(ILEAVE) POINTS TO THE EJECTED COLUMN.
C     USE ERD(*) FOR TEMP. STORAGE.
      CALL SCOPY(MRELAS,ZERO,0,ERD,1)
      ERD(ABS(ILEAVE))=ONE
      TRANS=.TRUE.
      CALL LA05BS(BASMAT,IBRC,LBM,MRELAS,IPR,IWR,WR,GG,ERD,TRANS)
C
C     COMPUTE UPDATED DUAL VARIABLES IN DUALS(*).
      NPR003=20042
      GO TO 30003
C
C     COMPUTE THE DOT PRODUCT OF COL. J OF THE NEW INVERSE (TRANSPOSE)
C     WITH EACH NON-BASIC COLUMN.  ALSO COMPUTE THE DOT PRODUCT OF THE
C     INVERSE (TRANSPOSE) OF NON-UPDATED MATRIX (TIMES) THE
C     SEARCH DIRECTION WITH EACH NON-BASIC COLUMN.
C     RECOMPUTE REDUCED COSTS.
20042 PAGEPL=.TRUE.
      CALL SCOPY(NVARS+MRELAS,ZERO,0,RZ,1)
      NNEGRC=0
      J=JSTRT
20043 IF (.NOT.(IBB(J).LE.0)) GO TO 20045
      PAGEPL=.TRUE.
      RG(J)=ONE
      GO TO 20046
C
C     NONBASIC INDEPENDENT VARIABLES (COLUMN IN SPARSE MATRIX STORAGE)
20045 IF (.NOT.(J.LE.NVARS)) GO TO 20048
      RZJ=COSTS(J)*COSTSC
      ALPHA=ZERO
      GAMMA=ZERO
C
C     COMPUTE THE DOT PRODUCT OF THE SPARSE MATRIX NONBASIC COLUMNS
C     WITH THREE VECTORS INVOLVED IN THE UPDATING STEP.
      IF (.NOT.(J.EQ.1)) GO TO 20051
      ILOW=NVARS+5
      GO TO 20052
20051 ILOW=IMAT(J+3)+1
20052 IF (.NOT.(PAGEPL)) GO TO 20054
      IL1=IPLOC(ILOW,AMAT,IMAT)
      IF (.NOT.(IL1.GE.LMX-1)) GO TO 20057
      ILOW=ILOW+2
      IL1=IPLOC(ILOW,AMAT,IMAT)
20057 CONTINUE
      IPAGE=ABS(IMAT(LMX-1))
      GO TO 20055
20054 IL1=IHI+1
20055 IHI=IMAT(J+4)-(ILOW-IL1)
20060 IU1=MIN(LMX-2,IHI)
      IF (.NOT.(IL1.GT.IU1)) GO TO 20062
      GO TO 20061
20062 CONTINUE
      DO 10 I=IL1,IU1
      RZJ=RZJ-AMAT(I)*DUALS(IMAT(I))
      ALPHA=ALPHA+AMAT(I)*ERD(IMAT(I))
      GAMMA=GAMMA+AMAT(I)*WW(IMAT(I))
10    CONTINUE
      IF (.NOT.(IHI.LE.LMX-2)) GO TO 20065
      GO TO 20061
20065 CONTINUE
      IPAGE=IPAGE+1
      KEY=1
      CALL PRWPGE(KEY,IPAGE,LPG,AMAT,IMAT)
      IL1=NVARS+5
      IHI=IHI-LPG
      GO TO 20060
20061 PAGEPL=IHI.EQ.(LMX-2)
      RZ(J)=RZJ*CSC(J)
      ALPHA=ALPHA*CSC(J)
      GAMMA=GAMMA*CSC(J)
      RG(J)=MAX(RG(J)-TWO*ALPHA*GAMMA+ALPHA**2*GQ,ONE+ALPHA**2)
C
C     NONBASIC DEPENDENT VARIABLES (COLUMNS DEFINED IMPLICITLY)
      GO TO 20049
20048 PAGEPL=.TRUE.
      SCALR=-ONE
      IF(IND(J).EQ.2) SCALR=ONE
      I=J-NVARS
      ALPHA=SCALR*ERD(I)
      RZ(J)=-SCALR*DUALS(I)
      GAMMA=SCALR*WW(I)
      RG(J)=MAX(RG(J)-TWO*ALPHA*GAMMA+ALPHA**2*GQ,ONE+ALPHA**2)
20049 CONTINUE
20046 CONTINUE
C
      RCOST=RZ(J)
      IF (MOD(IBB(J),2).EQ.0) RCOST=-RCOST
      IF (.NOT.(IND(J).EQ.3)) GO TO 20068
      IF(BU(J).EQ.BL(J)) RCOST=ZERO
20068 CONTINUE
      IF (IND(J).EQ.4) RCOST=-ABS(RCOST)
      CNORM=ONE
      IF (J.LE.NVARS) CNORM=COLNRM(J)
      IF (RCOST+ERDNRM*DULNRM*CNORM.LT.ZERO) NNEGRC=NNEGRC+1
      J=MOD(J,MRELAS+NVARS)+1
      IF (.NOT.(NNEGRC.GE.NPP .OR. J.EQ.JSTRT)) GO TO 20071
      GO TO 20044
20071 CONTINUE
      GO TO 20043
20044 JSTRT=J
C
C     UPDATE THE EDGE WEIGHT FOR THE EJECTED VARIABLE.
      RG(ABS(IBASIS(IENTER)))= GQ/WP**2
C
C     IF MINIMUM REDUCED COST (DANTZIG) PRICING IS USED,
C     CALCULATE THE NEW REDUCED COSTS.
      GO TO 20040
C
C     COMPUTE THE UPDATED DUALS IN DUALS(*).
20039 NPR003=20074
      GO TO 30003
20074 CALL SCOPY(NVARS+MRELAS,ZERO,0,RZ,1)
      NNEGRC=0
      J=JSTRT
      PAGEPL=.TRUE.
C
20075 IF (.NOT.(IBB(J).LE.0)) GO TO 20077
      PAGEPL=.TRUE.
      GO TO 20078
C
C     NONBASIC INDEPENDENT VARIABLE (COLUMN IN SPARSE MATRIX STORAGE)
20077 IF (.NOT.(J.LE.NVARS)) GO TO 20080
      RZ(J)=COSTS(J)*COSTSC
      IF (.NOT.(J.EQ.1)) GO TO 20083
      ILOW=NVARS+5
      GO TO 20084
20083 ILOW=IMAT(J+3)+1
20084 CONTINUE
      IF (.NOT.(PAGEPL)) GO TO 20086
      IL1=IPLOC(ILOW,AMAT,IMAT)
      IF (.NOT.(IL1.GE.LMX-1)) GO TO 20089
      ILOW=ILOW+2
      IL1=IPLOC(ILOW,AMAT,IMAT)
20089 CONTINUE
      IPAGE=ABS(IMAT(LMX-1))
      GO TO 20087
20086 IL1=IHI+1
20087 CONTINUE
      IHI=IMAT(J+4)-(ILOW-IL1)
20092 IU1=MIN(LMX-2,IHI)
      IF (.NOT.(IU1.GE.IL1 .AND.MOD(IU1-IL1,2).EQ.0)) GO TO 20094
      RZ(J)=RZ(J)-AMAT(IL1)*DUALS(IMAT(IL1))
      IL1=IL1+1
20094 CONTINUE
      IF (.NOT.(IL1.GT.IU1)) GO TO 20097
      GO TO 20093
20097 CONTINUE
C
C     UNROLL THE DOT PRODUCT LOOP TO A DEPTH OF TWO.  (THIS IS DONE
C     FOR INCREASED EFFICIENCY).
      DO 40 I=IL1,IU1,2
      RZ(J)=RZ(J)-AMAT(I)*DUALS(IMAT(I))-AMAT(I+1)*DUALS(IMAT(I+1))
40    CONTINUE
      IF (.NOT.(IHI.LE.LMX-2)) GO TO 20100
      GO TO 20093
20100 CONTINUE
      IPAGE=IPAGE+1
      KEY=1
      CALL PRWPGE(KEY,IPAGE,LPG,AMAT,IMAT)
      IL1=NVARS+5
      IHI=IHI-LPG
      GO TO 20092
20093 PAGEPL=IHI.EQ.(LMX-2)
      RZ(J)=RZ(J)*CSC(J)
C
C     NONBASIC DEPENDENT VARIABLES (COLUMNS DEFINED IMPLICITLY)
      GO TO 20081
20080 PAGEPL=.TRUE.
      SCALR=-ONE
      IF(IND(J).EQ.2) SCALR=ONE
      I=J-NVARS
      RZ(J)=-SCALR*DUALS(I)
20081 CONTINUE
20078 CONTINUE
C
      RCOST=RZ(J)
      IF (MOD(IBB(J),2).EQ.0) RCOST=-RCOST
      IF (.NOT.(IND(J).EQ.3)) GO TO 20103
      IF(BU(J).EQ.BL(J)) RCOST=ZERO
20103 CONTINUE
      IF (IND(J).EQ.4) RCOST=-ABS(RCOST)
      CNORM=ONE
      IF (J.LE.NVARS) CNORM=COLNRM(J)
      IF (RCOST+ERDNRM*DULNRM*CNORM.LT.ZERO) NNEGRC=NNEGRC+1
      J=MOD(J,MRELAS+NVARS)+1
      IF (.NOT.(NNEGRC.GE.NPP .OR. J.EQ.JSTRT)) GO TO 20106
      GO TO 20076
20106 CONTINUE
      GO TO 20075
20076 JSTRT=J
20040 CONTINUE
      GO TO 20030
C
C     THIS IS NECESSARY ONLY FOR PRINTING OF INTERMEDIATE RESULTS.
20029 NPR003=20109
      GO TO 30003
20109 CONTINUE
20030 RETURN
C     PROCEDURE (TRANSLATE RIGHT HAND SIDE)
C
C     PERFORM THE TRANSLATION ON THE RIGHT-HAND SIDE.
30001 IF (.NOT.(IBAS.LE.NVARS)) GO TO 20110
      I=0
20113 CALL PNNZRS(I,AIJ,IPLACE,AMAT,IMAT,IBAS)
      IF (.NOT.(I.LE.0)) GO TO 20115
      GO TO 20114
20115 CONTINUE
      RHS(I)=RHS(I)-SCALR*AIJ*CSC(IBAS)
      GO TO 20113
20114 GO TO 20111
20110 I=IBAS-NVARS
      IF (.NOT.(IND(IBAS).EQ.2)) GO TO 20118
      RHS(I)=RHS(I)-SCALR
      GO TO 20119
20118 RHS(I)=RHS(I)+SCALR
20119 CONTINUE
20111 CONTINUE
      RHSNRM=MAX(RHSNRM,SASUM(MRELAS,RHS,1))
C      GO TO NPR001, (20009,20013,20017,20028)
      if(NPR001.eq.20009)go to 20009
      if(NPR001.eq.20013)go to 20013
      if(NPR001.eq.20017)go to 20017
      if(NPR001.eq.20028)go to 20028
C     PROCEDURE (COMPUTE NEW PRIMAL)
C
C     COPY RHS INTO WW(*), SOLVE SYSTEM.
30002 CALL SCOPY(MRELAS,RHS,1,WW,1)
      TRANS = .FALSE.
      CALL LA05BS(BASMAT,IBRC,LBM,MRELAS,IPR,IWR,WR,GG,WW,TRANS)
      CALL SCOPY(MRELAS,WW,1,RPRIM,1)
      RPRNRM=SASUM(MRELAS,RPRIM,1)
      GO TO 20038
C     PROCEDURE (COMPUTE NEW DUALS)
C
C     SOLVE FOR DUAL VARIABLES. FIRST COPY COSTS INTO DUALS(*).
30003 I=1
      N20121=MRELAS
      GO TO 20122
20121 I=I+1
20122 IF ((N20121-I).LT.0) GO TO 20123
      J=IBASIS(I)
      IF (.NOT.(J.LE.NVARS)) GO TO 20125
      DUALS(I)=COSTSC*COSTS(J)*CSC(J) + XLAMDA*PRIMAL(I+NVARS)
      GO TO 20126
20125 DUALS(I)=XLAMDA*PRIMAL(I+NVARS)
20126 CONTINUE
      GO TO 20121
C
20123 TRANS=.TRUE.
      CALL LA05BS(BASMAT,IBRC,LBM,MRELAS,IPR,IWR,WR,GG,DUALS,TRANS)
      DULNRM=SASUM(MRELAS,DUALS,1)
C     GO TO NPR003, (20042,20074,20109)
      if(NPR003.eq.20042)go to 20042
      if(NPR003.eq.20074)go to 20074
      if(NPR003.eq.20109)go to 20109
      END
