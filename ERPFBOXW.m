function ACCEL = ERPFBOXW(ERM,ANT,GRD,REFF,YSAT,SUN,KAPPA,MONTH,SVN,MJD)
% NAME       :  ERPFBOXW
%
% PURPOSE    :  COMPUTATION OF EARTH RADIATION PRESSURE ACTING ON A
%               BOW-WING SATELLITE
%
% PARAMETERS :
%         IN :  ERM       : EARTH RADIATION MODEL
%                           0 = NONE
%                           1 = EARTH RADIATION PRESSURE (ANALYTICAL)
%                           2 = EARTH RADIATION PRESSURE (CERES DATA)
%               CERES DAT : CERES DATA (EXTERNAL ASCII FILES) -> SET PATH IN DATPATH
%               ANT       : 0 = NO ANTENNA THRUST
%                         : 1 = WITH ANTENNA THRUST
%               GRD       : 1 = 2.5 degree grid
%                         : 2 = 5.0 degree grid
%                         : 3 = 10  degree grid
%               REFF      : ACCELERATION IN REFERENCE FRAME:
%                           0 = INERTIAL
%                           1 = BODY FIXED (Z,Y,X)
%                           2 = SUN FIXED (D,Y,B)
%                           3 = ORBITAL (RADIAL, ALONG- AND CROSS-TRACK)
%               YSAT      : SATELLITE POSITION [m] AND VELOCITY [m/s] (INERTIAL), 
%                           (POSITION,VELOCITY) = (RX,RY,RZ,VX,VY,VZ)
%               SUN       : SUN POSITION VECTOR [m] (INERTIAL)
%               KAPPA     : ROTATION MATRIX FROM INERTIAL TO EARTH-FIXED
%                           REFERENCE FRAME (FOR ERM=2)
%               MONTH     : MONTH NUMBER (1 ... 12), FOR ERM=2
%               BLKID     : BLOCK Identifier
%                             1 = GPS-I
%                             2 = GPS-II
%                             3 = GPS-IIA
%                             4 = GPS-IIR
%                             5 = GPS-IIR-A
%                             6 = GPS-IIR-B
%                             7 = GPS-IIR-M
%                             8 = GPS-IIF
%                             9 = GPS-IIIA  (Updated from acc_albedo_propboxw.f)
%                           101 = GLONASS
%                           102 = GLONASS-M (Added TAH 190702)
%                           103 = GLONASS-K (Added TAH 190702)
%                           201 = Galileo (IOV) (Added from acc_albedo_propboxw.f)
%                           202 = Galileo (FOC) (Added from acc_albedo_propboxw.f)
%                           301 = BDS-2 GEO
%                           302 = BDS-2 IGSO
%                           303 = BDS-2 MEO
%                           304 = BDS-3 MEO-CAST
%                           305 = BDS-3 MEO-SECM
%                           306 = BDS-3 IGSO
%                           307 = BDS-3 GEO
%                           401 = QZSS-1
%                           402 = QZSS-2I
%                           403 = QZSS-2G
%                           NB: when adding a new block number, check last block in this function%
%               SVN       : SPACE VEHICLE NUMBER          
%               MJD       : MODIFIED JULIAN DAY
%
%        OUT : ACCEL      : ACCELERATION VECTOR [m/s^2]

      global Cslight
      
      % Initialization of force vector
      FORCE = [0.0, 0.0, 0.0];
      ACCEL = [0.0, 0.0, 0.0];

      % --------------------------
      % NOMINAL SATELLITE ATTITUDE
      % --------------------------

      ABSPOS = sqrt(YSAT(1)^2+YSAT(2)^2+YSAT(3)^2);
      for K=1:3
         RADVEC(K) = YSAT(K)/ABSPOS;
      end
      
      CRSVEC(1) = YSAT(2)*YSAT(6)-YSAT(3)*YSAT(5);
      CRSVEC(2) = YSAT(3)*YSAT(4)-YSAT(1)*YSAT(6);
      CRSVEC(3) = YSAT(1)*YSAT(5)-YSAT(2)*YSAT(4);
      ABSCRS = DSQRT(CRSVEC(1)^2+CRSVEC(2)^2+CRSVEC(3)^2);
      for K=1:3
         CRSVEC(K) = CRSVEC(K)/ABSCRS;
      end

      ALGVEC(1) = CRSVEC(2)*RADVEC(3)-CRSVEC(3)*RADVEC(2);
      ALGVEC(2) = CRSVEC(3)*RADVEC(1)-CRSVEC(1)*RADVEC(3);
      ALGVEC(3) = CRSVEC(1)*RADVEC(2)-CRSVEC(2)*RADVEC(1); 

      if(ERM > 0)
          % DISTANCE FROM SATELLITE TO SUN
          RSUN = sqrt((YSAT(1)-SUN(1))^2+(YSAT(2)-SUN(2))^2+(YSAT(3)-SUN(3))^2);

          % D VECTOR AND Z VECTOR
          for K=1:3
             D_SUN(K) = (SUN(K)-YSAT(K))/RSUN;
             Z_SAT(K) = -YSAT(K)/ABSPOS;
          end

          % Y VECTOR
          Y0SAT(1) = Z_SAT(2)*D_SUN(3)-Z_SAT(3)*D_SUN(2);
          Y0SAT(2) = Z_SAT(3)*D_SUN(1)-Z_SAT(1)*D_SUN(3);
          Y0SAT(3) = Z_SAT(1)*D_SUN(2)-Z_SAT(2)*D_SUN(1);
          ABSY0V = sqrt(Y0SAT(1)^2 + Y0SAT(2)^2 + Y0SAT(3)^2);
          for K=1:3
             Y_SAT(K) = Y0SAT(K)/ABSY0V;
          end

          % B VECTOR
          B_SUN(1) = Y_SAT(2)*D_SUN(3) - Y_SAT(3)*D_SUN(2);
          B_SUN(2) = Y_SAT(3)*D_SUN(1) - Y_SAT(1)*D_SUN(3);
          B_SUN(3) = Y_SAT(1)*D_SUN(2) - Y_SAT(2)*D_SUN(1);

          % X VECTOR
          X_SAT(1) = Y_SAT(2)*Z_SAT(3) - Y_SAT(3)*Z_SAT(2);
          X_SAT(2) = Y_SAT(3)*Z_SAT(1) - Y_SAT(1)*Z_SAT(3);
          X_SAT(3) = Y_SAT(1)*Z_SAT(2) - Y_SAT(2)*Z_SAT(1);

          for K=1:3
             ATTSURF(K,1) = Z_SAT(K);
             ATTSURF(K,2) = Y_SAT(K);
             ATTSURF(K,3) = X_SAT(K);
             ATTSURF(K,4) = D_SUN(K);
          end

          % =============================================================================
          % FORM A ORBIT-NORMAL ATTITUDE FOR BDS GEO SATELLITE (ONLY FOR TESTING PURPOSE).
          % BECAUSE THE IMPACT OF EARTH ALBEDO IS MOSTLY ON THE RADIAL DIRECTION,
          % THE CONTRIBUTIONS FROM X AND Y SIDES ARE ZERO NO MATTER WHAT THE
          % SATELLITE IS OPERATED WITH YS OR ON ATTITUDE. THIS IS BECAUSE THE
          % EARTH RADIATION DIRECTION (RADIAL DIRECTION) IS PERPENDICULAR TO THE
          % OTHER TWO AXES (Y AND X). THE ACCELERATION IS ONLY
          % RESULTED FROM THE Z SIDE AND THE SOLAR PANEL. IN A CASE OF SMALL BETA ANGLES, 
          % THE DIFFERENCE IN ACCELERATION AT THE SOLAR PANEL BETWEEN YS AND ON ATTITUDES 
          % IS A COSINE FUNCTION OF THE BETA ANGLE. AS A RESULT, THE IMPACT
          % OF EARTH ALBEDO ON YAW-STERRING IS SIMILAR TO THAT ON ORBIT-NORMAL AT
          % THE SMALL BETA ANGLES. (04-11-2019 Dr. TZUPANG TSENG)
          %
          % FOR BDS IGSO ECLIPSING PERIOD (BETA < 4 DEG), THE CONTRIBUTION OF
          % COSINE(4DEG) FROM THE SOLAR PANEL IN ON ATTITUDE IS SIMILAR TO THAT IN
          % YS ATTITUDE. (05-11-2019 Dr. TZUPANG TSENG)
          % =============================================================================
          if (BLKID == 301 || BLKID == 307)
              % Y VECTOR
              Y_SAT(1) = -CRSVEC(1); 
              Y_SAT(2) = -CRSVEC(2); 
              Y_SAT(3) = -CRSVEC(3);
 
              % X VECTOR
              X_SAT(1) = Y_SAT(2)*Z_SAT(3) - Y_SAT(3)*Z_SAT(2);
              X_SAT(2) = Y_SAT(3)*Z_SAT(1) - Y_SAT(1)*Z_SAT(3);
              X_SAT(3) = Y_SAT(1)*Z_SAT(2) - Y_SAT(2)*Z_SAT(1);
 
              % B VECTOR
              B_SUN(1) = Y_SAT(2)*D_SUN(3) - Y_SAT(3)*D_SUN(2);
              B_SUN(2) = Y_SAT(3)*D_SUN(1) - Y_SAT(1)*D_SUN(3);
              B_SUN(3) = Y_SAT(1)*D_SUN(2) - Y_SAT(2)*D_SUN(1);
 
              % D VECTOR RESULTED FROM ORBIT-NORMAL ATTITUDE
              D_SUN(1) = Y_SAT(2)*B_SUN(3) - Y_SAT(3)*B_SUN(2);
              D_SUN(2) = Y_SAT(3)*B_SUN(1) - Y_SAT(1)*B_SUN(3);
              D_SUN(3) = Y_SAT(1)*B_SUN(2) - Y_SAT(2)*B_SUN(1);
 
              for K=1:3
                 ATTSURF(K,1) = Z_SAT(K);
                 ATTSURF(K,2) = Y_SAT(K);
                 ATTSURF(K,3) = X_SAT(K);
                 ATTSURF(K,4) = D_SUN(K);
              end 
          end
      end

      % ---------------------------- 
      % OPTICAL PROPERTIES PER BLOCK
      % ----------------------------

      INDB=1;

      if (ERM > 0)
         if (BLKID <= 10)
             INDB = BLKID;
         elseif (BLKID > 100  &&  BLKID < 200 )
             INDB = BLKID - 90;
         elseif (BLKID > 200  &&  BLKID < 300 )   % MOD TAH 190722: Added Galileo
             INDB = BLKID - 180;
         elseif (BLKID > 300  &&  BLKID < 400 )   % MOD SCM 191219: Added BDS
             INDB = BLKID - 270;
         elseif (BLKID > 400)   % MOD SCM 191219: Added QZSS         
             INDB = BLKID - 360;
         end
         for II = 1:4
            for JJ = 1:2
               AREAS(II,JJ) = AREA(II,JJ,INDB);
               REFLS(II,JJ) = REFL(II,JJ,INDB);
               DIFUS(II,JJ) = DIFU(II,JJ,INDB);
               ABSPS(II,JJ) = ABSP(II,JJ,INDB);

               AREA2S(II,JJ) = AREA2(II,JJ,INDB);
               REFL2S(II,JJ) = REFL2(II,JJ,INDB);
               DIFU2S(II,JJ) = DIFU2(II,JJ,INDB);
               ABSP2S(II,JJ) = ABSP2(II,JJ,INDB);

               REFLIRS(II,JJ) = REFLIR(II,JJ,INDB);
               DIFUIRS(II,JJ) = DIFUIR(II,JJ,INDB);
               ABSPIRS(II,JJ) = ABSPIR(II,JJ,INDB);
            end
         end
      end

      % ----------------------
      % EARTH RADIATION MODELS
      % ----------------------
      if (ERM > 0)
         ABSSUN = sqrt(SUN(1)^2 + SUN(2)^2 + SUN(3)^2);
         for K=1:3
            ESUN(K) = SUN(K)/ABSSUN;
         end

         PSIDOT = ESUN(1)*RADVEC(1)+ESUN(2)*RADVEC(2)+ESUN(3)*RADVEC(3);
         if (abs(PSIDOT) > (1.0 - 1e-6))
            PSI = 0.0;
         else
            PSI = acos(PSIDOT);
         end
         S1 = S0*(AU/ABSSUN)^2;
      end

      %  ANALYTICAL MODEL
      if (ERM == 1)
         NCFVEC(1) = RADVEC(1);
         NCFVEC(2) = RADVEC(2);
         NCFVEC(3) = RADVEC(3);
         ALBFAC = (pi*TOA^2)*(S1/cslight)/(ABSPOS^2);
         PHASEVI = (2*ALB/(3*pi^2))*((pi-PSI)*DCOS(PSI)+DSIN(PSI));
         PHASEIR = (1-ALB)/(4*pi);
         ABSNCFVI = ALBFAC*PHASEVI;
         ABSNCFIR = ALBFAC*PHASEIR;

         FORCE = SURFBOXW(AREAS,REFLS,DIFUS,ABSPS,AREA2S,REFL2S,DIFU2S,ABSP2S, ...
                           REFLIRS,DIFUIRS,ABSPIRS,ABSNCFVI,ABSNCFIR,NCFVEC,ATTSURF);


      % NUMERICAL MODEL (CERES DATA)
      elseif (ERM == 2)

         ABSSUN = sqrt(SUN(1)^2 + SUN(2)^2 + SUN(3)^2);
         for K=1:3
            ESUN(K) = SUN(K)/ABSSUN;
         end

         for LATK = 1:LATKMX
            for LONK = 1:LONKMX
               D_AREA = D_AREA_ALL(LATK,LONK);
               V_NS(1) = V_NS_ALL(LATK,LONK,1);
               V_NS(2) = V_NS_ALL(LATK,LONK,2);
               V_NS(3) = V_NS_ALL(LATK,LONK,3);

               for II=1:3
                  V_INS(II)=0.0;
                  for JJ=1:3
                     V_INS(II) = V_INS(II) + KAPPA(JJ,II)*V_NS(JJ);
                  end
               end

               % Distance and direction from point in the Earth to satellite
               V_DIST(1) = YSAT(1)-TOA*V_INS(1); 
               V_DIST(2) = YSAT(2)-TOA*V_INS(2);
               V_DIST(3) = YSAT(3)-TOA*V_INS(3);
               DIST2 = V_DIST(1)^2 +V_DIST(2)^2 +V_DIST(3)^2;
               ABSDIST = sqrt(DIST2);
               V_SAT(1) = V_DIST(1)/ABSDIST;
               V_SAT(2) = V_DIST(2)/ABSDIST;
               V_SAT(3) = V_DIST(3)/ABSDIST;

               % Cosine of angles of incident and reflected radiation
               COS_IN = ESUN(1)*V_INS(1) + ESUN(2)*V_INS(2) + ESUN(3)*V_INS(3);
               COS_RE = V_SAT(1)*V_INS(1) +V_SAT(2)*V_INS(2) + V_SAT(3)*V_INS(3);

               if (COS_RE >= 0)
                   % Reflectivity and emissivity coefficients
                   REFL_CF = CERGRE(LATK,LONK);
                   EMIT_CF = CERGEM(LATK,LONK);

                   % Reflected Irradiance
                   if (COS_IN >= 0)
                       E_REFL=(REFL_CF/(pi*DIST2))*COS_RE*COS_IN*S1*D_AREA;
                   else
                       E_REFL=0.0;
                   end

                   % Emitted Irradiance
                   E_EMIT = (EMIT_CF/(4*pi*DIST2))*COS_RE*S1*D_AREA;

                   % Non-conservative force
                   ABSNCFVI = E_REFL/Cslight;
                   ABSNCFIR = E_EMIT/Cslight;
                   NCFVEC(1) = V_SAT(1);
                   NCFVEC(2) = V_SAT(2);
                   NCFVEC(3) = V_SAT(3);

                   FORCE_LL = SURFBOXW(AREAS,REFLS,DIFUS,ABSPS,AREA2S,REFL2S,DIFU2S,ABSP2S, ...
                                         REFLIRS,DIFUIRS,ABSPIRS,ABSNCFVI,ABSNCFIR,NCFVEC,ATTSURF);

                   FORCE(1) = FORCE(1) + FORCE_LL(1);
                   FORCE(2) = FORCE(2) + FORCE_LL(2);
                   FORCE(3) = FORCE(3) + FORCE_LL(3);
               end
            end
         end
      end


%     ANTENNA POWER OF GPS SATELLITES (IN WATTS)
%     IGS MODEL (JIM RAY, 2011)

%     NO ANTENNA POWER INFORMATION FOR GLONASS SATELLITES
%     initialisation of ANTPOW
      ANTPOW = 0.0;

%     GPS BLOCK IIA (ASSUMED THE SAME FOR BLOCK I AND II) 
      if (BLKID <= 3)
          ANTPOW = 76.D0;

%     GPS BLOCK IIR
      elseif ((BLKID >= 4) && (BLKID <= 6))
          ANTPOW = 85.D0;

%     GPS BLOCK IIR-M
      elseif (BLKID == 7)
          ANTPOW = 198.D0;

%     GPS BLOCK IIF
      elseif (BLKID == 8)
          ANTPOW = 249.D0;
          if ((SVN == 62) && (MJD >= 55656.D0))
              % NO M-CODE FOR SVN62/PRN25 STARTING 05APR2011; NANU 2011026
              ANTPOW = 154.D0;
          end
      end



%     NAVIGATION ANTENNA THRUST (SIMPLE MODEL)
      if ((ANT == 1) && (BLKID < 100))
          ANTFORCE = ANTPOW/Cslight;
          FORCE(1) = FORCE(1) + ANTFORCE*RADVEC(1);
          FORCE(2) = FORCE(2) + ANTFORCE*RADVEC(2);
          FORCE(3) = FORCE(3) + ANTFORCE*RADVEC(3);
      end

      if ((REFF > 0) && ((ERM > 0)||(ANT == 1)))
          for K=1:3
              FREF(K) = 0.D0;
          end

          % FORCE IN BODY-FIXED REFERENCE FRAME
          if (REFF == 1)
             for K=1:3
                FREF(1) = FREF(1) + FORCE(K)*Z_SAT(K);
                FREF(2) = FREF(2) + FORCE(K)*Y_SAT(K);
                FREF(3) = FREF(3) + FORCE(K)*X_SAT(K);
             end

          % FORCE IN SUN-FIXED REFERENCE FRAME
          elseif (REFF == 2)
              for K=1:3
                  FREF(1) = FREF(1) + FORCE(K)*D_SUN(K);
                  FREF(2) = FREF(2) + FORCE(K)*Y_SAT(K);
                  FREF(3) = FREF(3) + FORCE(K)*B_SUN(K);
              end

          %  FORCE IN ORBITAL REFERENCE FRAME
          elseif (REFF == 3)
              for K=1:3
                  FREF(1) = FREF(1) + FORCE(K)*RADVEC(K);
                  FREF(2) = FREF(2) + FORCE(K)*ALGVEC(K);
                  FREF(3) = FREF(3) + FORCE(K)*CRSVEC(K);
              end
          end

          for K=1:3
             FORCE(K) = FREF(K);
          end
      end

      % CONVERSION TO ACCELERATION
      for K=1:3
         ACCEL(K) = FORCE(K)/MASS;
      end
end