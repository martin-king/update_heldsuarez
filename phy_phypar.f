      SUBROUTINE PHYPAR (VOR1,DIV1,T1,Q1,PHI1,PSL1,
     &                   UTEND,VTEND,TTEND,QTEND)
C--
C--   SUBROUTINE PHYPAR (VOR1,DIV1,T1,Q1,PHI1,PSL1,
C--  &                   UTEND,VTEND,TTEND,QTEND)
C--
C--   Purpose: compute physical parametrization tendencies for U, V, T, Q 
C--   and add them to dynamical grid-point tendencies
C--   Input-only  arguments:   VOR1   : vorticity (sp)
C--                            DIV1   : divergence (sp)
C--                            T1     : temperature (sp)
C--                            Q1     : specific humidity (sp)
C--                            PHI1   : geopotential (sp)
C--                            PSL1   : log of sfc pressure (sp)
C--   Input-output arguments:  UTEND  : u-wind tendency (gp)
C--                            VTEND  : v-wind tendency (gp)
C--                            TTEND  : temp. tendency (gp)
C--                            QTEND  : spec. hum. tendency (gp)
C--   Modified common blocks:  PHYGR1, PHYGR2, PHYGR3, PHYTEN, FLUXES
C--

C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Constants + functions of sigma and latitude
      include "com_physcon.h"

C     Model variables, tendencies and fluxes on gaussian grid
      include "com_physvar.h"

C     Surface properties (time-inv.)
      include "com_surfcon.h"

C     Surface fields (daily averages)
      include "com_cli_sea.h"
      include "com_cli_land.h"
      include "com_var_sea.h"
      include "com_var_land.h"

C     Logical and coupling flags
      include "com_lflags.h"
      include "com_cpl_flags.h"

      COMPLEX VOR1(MX,NX,NLEV), DIV1(MX,NX,NLEV), T1(MX,NX,NLEV),
     &          Q1(MX,NX,NLEV), PHI1(MX,NX,NLEV), PSL1(MX,NX),
     &          UCOS(MX,NX), VCOS(MX,NX)

      REAL UTEND(NGP,NLEV), VTEND(NGP,NLEV), TTEND(NGP,NLEV),
     &     QTEND(NGP,NLEV)

      INTEGER IPTOP(NGP), ICLTOP(NGP,2), ICNV(NGP)
      REAL    RPS(NGP), GSE(NGP)

      iitest=0

C--   1. Compute grid-point fields

C     1.1 Convert model spectral variables to grid-point variables

      if (iitest.eq.1) print *, ' 1.1 in PHYPAR'

      DO K=1,NLEV

        CALL UVSPEC (VOR1(1,1,K),DIV1(1,1,K),UCOS,VCOS)
        CALL GRID   (UCOS,UG1(1,K),2)
        CALL GRID   (VCOS,VG1(1,K),2)

      ENDDO

      DO K=1,NLEV

        CALL GRID   (T1(1,1,K),  TG1(1,K),  1)
        CALL GRID   (Q1(1,1,K),  QG1(1,K),  1)
        CALL GRID   (PHI1(1,1,K),PHIG1(1,K),1)

      ENDDO

      CALL GRID (PSL1,PSLG1,1)

C     Remove negative humidity values
C     CALL QNEG (QG1)

C     1.2 Compute thermodynamic variables

      if (iitest.eq.1) print *, ' 1.2 in PHYPAR'

      DO J=1,NGP
       PSG(J)=EXP(PSLG1(J))
       RPS(J)=1./PSG(J)
      ENDDO

      DO K=1,NLEV
       DO J=1,NGP
c       remove when QNEG is implemented
	qg1(j,k)=max(qg1(j,k),0.)
        SE(J,K)=CP*TG1(J,K)+PHIG1(J,K)
       ENDDO
      ENDDO

      DO K=1,NLEV
       CALL SHTORH (1,NGP,TG1(1,K),PSG,SIG(K),QG1(1,K),
     &              RH(1,K),QSAT(1,K))
      ENDDO

C--   2. Precipitation 

C     2.1 Deep convection

cmpk      CALL CONVMF (PSG,SE,QG1,QSAT,
cmpk     &             IPTOP,CBMF,PRECNV,TT_CNV,QT_CNV)
cmpk
cmpk      DO K=2,NLEV
cmpk       DO J=1,NGP
cmpk        TT_CNV(J,K)=TT_CNV(J,K)*RPS(J)*GRDSCP(K)
cmpk        QT_CNV(J,K)=QT_CNV(J,K)*RPS(J)*GRDSIG(K)
cmpk       ENDDO
cmpk      ENDDO

cmpk      DO J=1,NGP
cmpk        ICNV(J)=NLEV-IPTOP(J)
cmpk      ENDDO

C     2.2 Large-scale condensation

cfk#if !defined(KNMI)
cmpk      CALL LSCOND (PSG,QG1,QSAT,
cmpk     &             IPTOP,PRECLS,TT_LSC,QT_LSC)
cfk#else
cmpk      CALL LSCOND (PSG,QG1,QSAT,TS,
cfk     &             IPTOP,PRECLS,SNOWLS,TT_LSC,QT_LSC)
cfk#endif

cmpk      DO K=2,NLEV
cmpk       DO J=1,NGP
cmpk        TTEND(J,K)=TTEND(J,K)+TT_CNV(J,K)+TT_LSC(J,K)
cmpk        QTEND(J,K)=QTEND(J,K)+QT_CNV(J,K)+QT_LSC(J,K)
cmpk       ENDDO
cmpk      ENDDO


C--   3. Radiation (shortwave and longwave) and surface fluxes

C     3.1 Compute shortwave tendencies and initialize lw transmissivity

cmpk      if (iitest.eq.1) print *, ' 3.1 in PHYPAR'

C     The sw radiation may be called at selected time steps

cmpk      IF (LRADSW) THEN

cmpk        DO J=1,NGP
cmpk          GSE(J) = (SE(J,NLEV-1)-SE(J,NLEV))/
cmpk     &             (PHIG1(J,NLEV-1)-PHIG1(J,NLEV))
cmpk        ENDDO

cmpk        CALL CLOUD (QG1,RH,PRECNV,PRECLS,IPTOP,GSE,FMASK1,
cmpk     &              ICLTOP,CLOUDC,CLSTR)

cmpk        DO J=1,NGP
cmpk          CLTOP(J)=SIGH(ICLTOP(J,1)-1)*PSG(J)
cmpk          PRTOP(J)=float(IPTOP(J))
cmpk        ENDDO

cmpk        CALL RADSW (PSG,QG1,ICLTOP,CLOUDC,CLSTR,
cmpk     &              SSRD,SSR,TSR,TT_RSW)

cmpk        DO K=1,NLEV
cmpk         DO J=1,NGP
cmpk          TT_RSW(J,K)=TT_RSW(J,K)*RPS(J)*GRDSCP(K)
cmpk         ENDDO
cmpk        ENDDO

cmpk      ENDIF

C     3.2 Compute downward longwave fluxes 

cmpk      CALL RADLW (-1,TG1,TS,
cmpk     &            SLRD,SLRU(1,3),
cmpk     &            SLR,OLR,TT_RLW)

C     3.3. Compute surface fluxes and land skin temperature

cmpk      if (iitest.eq.1) then 
cmpk         print *, ' 3.3 in PHYPAR'
cmpk         print *, 'mean(STL_AM) =', sum(STL_AM(:))/ngp
cmpk         print *, 'mean(SST_AM) =', sum(SST_AM(:))/ngp
cmpk      endif

cmpk      CALL SUFLUX (PSG,UG1,VG1,TG1,QG1,RH,PHIG1,
cmpk     &             PHIS0,FMASK1,STL_AM,SST_AM,SOILW_AM,SSRD,SLRD,
cmpk     &             USTR,VSTR,SHF,EVAP,SLRU,HFLUXN,
cmpk     &             TS,TSKIN,U0,V0,T0,Q0,.true.)
     
C--  3.3.1. Recompute sea fluxes in case of anomaly coupling

cmpk      IF (ICSEA .GT. 0) THEN 

cmpk         CALL SUFLUX (PSG,UG1,VG1,TG1,QG1,RH,PHIG1,
cmpk     &             PHIS0,FMASK1,STL_AM,SSTI_OM,SOILW_AM,SSRD,SLRD,
cmpk     &             USTR,VSTR,SHF,EVAP,SLRU,HFLUXN,
cmpk     &             TS,TSKIN,U0,V0,T0,Q0,.false.)

cmpk      ENDIF   

C     3.4 Compute upward longwave fluxes, convert them to tendencies 
C         and add shortwave tendencies

cmpk      if (iitest.eq.1) print *, ' 3.4 in PHYPAR'

cmpk      CALL RADLW (1,TG1,TS,
cmpk     &            SLRD,SLRU(1,3),
cmpk     &            SLR,OLR,TT_RLW)

cmpk      DO K=1,NLEV
cmpk       DO J=1,NGP
cmpk        TT_RLW(J,K)=TT_RLW(J,K)*RPS(J)*GRDSCP(K)
cmpk        TTEND (J,K)=TTEND(J,K)+TT_RSW(J,K)+TT_RLW(J,K)
cmpk       ENDDO
cmpk      ENDDO

C--   4. PBL interactions with lower troposphere

C     4.1 Vertical diffusion and shallow convection

cmpk      CALL VDIFSC (UG1,VG1,SE,RH,QG1,QSAT,PHIG1,ICNV,
cmpk     &             UT_PBL,VT_PBL,TT_PBL,QT_PBL)

C     4.2 Add tendencies due to surface fluxes 

cmpk      DO J=1,NGP
cmpk       UT_PBL(J,NLEV)=UT_PBL(J,NLEV)+USTR(J,3)*RPS(J)*GRDSIG(NLEV)
cmpk       VT_PBL(J,NLEV)=VT_PBL(J,NLEV)+VSTR(J,3)*RPS(J)*GRDSIG(NLEV)
cmpk       TT_PBL(J,NLEV)=TT_PBL(J,NLEV)+ SHF(J,3)*RPS(J)*GRDSCP(NLEV)
cmpk       QT_PBL(J,NLEV)=QT_PBL(J,NLEV)+EVAP(J,3)*RPS(J)*GRDSIG(NLEV)
cmpk      ENDDO

cmpk      DO K=1,NLEV
cmpk       DO J=1,NGP
cmpk        UTEND(J,K)=UTEND(J,K)+UT_PBL(J,K)
cmpk        VTEND(J,K)=VTEND(J,K)+VT_PBL(J,K)
cmpk        TTEND(J,K)=TTEND(J,K)+TT_PBL(J,K)
cmpk        QTEND(J,K)=QTEND(J,K)+QT_PBL(J,K)
cmpk       ENDDO
cmpk      ENDDO

C--   5. Store all fluxes for coupling and daily-mean output

      CALL DMFLUX (1)

C--   6. Random diabatic forcing 

      IF (LRANDF) THEN

C       6.1 Compute zonal-mean cross sections of diabatic forcing

        IF (LRADSW) THEN
          CALL XS_RDF (TT_LSC,TT_CNV,1)
          CALL XS_RDF (TT_RSW,TT_RLW,2)
        ENDIF

C--     6.2 Compute and store 3-D pattern of random diabatic forcing

        DO K=1,NLEV
         DO J=1,NGP
          TT_CNV(J,K)=TT_CNV(J,K)+TT_LSC(J,K)
         ENDDO
        ENDDO

        CALL SETRDF (TT_LSC)

        DO K=1,NLEV
         DO J=1,NGP
          TTEND(J,K)=TTEND(J,K)+TT_LSC(J,K)
         ENDDO
        ENDDO

      ENDIF

cmpk- copied cmpkhs modifications
cmpkhs Held-Suarez 1994 BAMS forcings; refer to table on p. 1826
cmpkhs begin
      DO K=1,NLEV
       DO J=1,NGP
cmpkhs index for latitude
        JLAT=(J+NLON-1)/NLON
cmpkhs
        temptemp=(SIG(K)-0.7d0)*3.333333333d0
        smallkt=2.893518519d-7+2.604166667d-6*max(0.d0,temptemp)
     &          *CLAT(JLAT)**4.d0
cmpkhs
        temptemp=(315.d0-60.d0*SLAT(JLAT)**2.d0-10.d0*log(SIG(K)
     &      *PSG(J))*CLAT(JLAT)**2.d0)*(SIG(K)*PSG(J))
     &      **2.857142857d-1
        Teq=max(200.,temptemp)
        TTEND(J,K)=TTEND(J,K)-smallkt*(TG1(J,K)-Teq)
cmpkhs no need to do moisture, so fixed to 0
        QTEND(J,K)=0.d0
       ENDDO
      ENDDO
cmpkhs
      DO K=1,NLEV
       DO J=1,NGP
        temptemp=(SIG(K)-0.7)*3.333333333
        smallkv=1.157407407d-5*max(0.d0,temptemp)
        VTEND(J,K)=VTEND(J,K)-smallkv*VG1(J,K)
        UTEND(J,K)=UTEND(J,K)-smallkv*UG1(J,K)
       ENDDO
      ENDDO
cmpkhs end

C--
      RETURN
      END

      include "phy_setrdf.f"

