function [IECLIPS, yaw_angle, MU, MURATE] = ECLIPS(IDIR, IPRN, TTAG, SVBCOS, ANOON, ANIGHT, ...
                                                     NECLIPS, ECLSTM, ECLETM,xsv, santxyz, vsvc, beta, iblk)
%
%       NAME            ECLIPS (version De%  2013)
%
%     PURPOSE   DETECT ECLIPSING & YAW ROTATE ECLIP. SATELLITES
%                       (THE INPUT BODY-X UNIT VECTOR - SANTXYZ IS YAW-ROTATED
%                        BY PHI-YANGLE (THE ECL-NOMINAL) YAW ANGLE DifFERENCE)  
%
%       COPYRIGHT       GEODETI% SURVEY DIVISION, 2011.
%                       ALL RIGHTS RESERVED.
%                       ALL TERMS AND CONDITIONS APPLY AS DETAILED IN
%                       " TERMS AND CONDITIONS FOR SOFTWARE " 
%
%       CONTACT         kouba@geod.nrcan.gc.ca
%
%       UPDATE HISTORY: Aug. 23, 2011:1996-2008 W. AVERAGES of JPL REPROCESSING
%                                YRATE SOLUTIONS FOR ALL II/IIA CODED IN DATA
%                                STATEMENT, THIS ENABLES REPROCESSING FROM 1996
%                       Sep 26, 2011: 1. Corrected bug causing Block Iif shadow 
%                               CROSSING MANEVURE WITH 0.06 DEG/SE% RATE EVEN 
%                               FOR BETADG > 8 DEG
%                                     2. CORRECTED/improved IIA RECOVERY logic
% De% 12, 2013 - start
%                       De% 18, 2013: 1.corrected Iif night turns (USAF Doc.)
%                                     2.small neg beta Iif and small pos IIA noon turn
%                                    (wrong) directions  for |beta| < 0.9deg
%                                     3. PRN/SVN 23 IIA YBIAS= -0.5 deg up to Feb 2004 
%                                     4. All the above changes labeled "% De% 12, 2013"
% De% 12, 2013 - end
%
%                      Jan 24, 2014:  NOON RECOVERY CORRECTED if IIA/Iif's HAVE YBIAS=0
%                                    (NOT APPLICABLE CURRENTLY,POSSIBLE FOR FUTURE Iif?)
%                                    & SMALL IIA/Iif BETA if STATEMENT SIMPLifIED
%
%     PARAMETERS        DESCRIPTION
%
%        IDIR           DIRECTION OF PROCESSING (1=FORWARD, -1=BACKWARD)
%        IPRN           SV PRN NUMBER (<=32 FOR GPS,  > 32 FOR GLONASS)
%        TTAG           OBSERVATION EPOCH TIME TAG (EG SE% OF GPS WEEK)
%        SVBCOS         SV BETA ANGLE (BETWEEN SV AND SUN RADIUS) COSINE
%        ANOON          SV BETA ANGLE LIMIT (DEG) FOR A NOON TURN MANEVURE 
%        ANIGHT         SV BETA ANGLE LIMIT (DEG) FOR A NIGHT SHADOW CROSSING
%        NECLIPS        NUMBER OF ECLIPSING FOR THE PRN SATELLITE 
%        ECLSTM         ECLIPSE START TIME(EG SE% OF GPS WEEK)
%        ECLETM         ECLIPSE END TIME  ( "         "      )
%        IECLIPS        SV ECLIPSING (0=NO,1, 2=YES(1=night, 2=noon))
%        PI             = PI=3.1415926536D0
%        XSV(3)         SV X, Y, Z (m)(ITRF)
%        SANTXYZ        BODY-X UNIT VECTOR (ITRF)
%                       WARNING: THE IIA BODY-X ORIENTATION EXPECTED FOR THE IIR
%                        THE  BODY-X REVERSED FOR IIR (SEE BELOW) & RETURNED
%        VSV%           SV INERTIAL VELOCIRY VECTOR IN ITRF
%        BETA           90 DEG + THE SUN ANGLE WITH ORBITAL PLANE(IN RAD)
%        IBLK           SV BLOCK  1=I, 2=II, 3=IIA, IIR=(4, 5) Iif=6
%
%        INTERNAL PARAMETRS DESCRIPTION
% De% 12, 2013
%        YBIAS       IIA YAW BIAS= 0.5 deg SINCE NOV95,EXCEPT PRN/SVN23?
%        YANGLE         THE NOMINAL YAW ANGLE
%        PHI            THE ECLIPSING YAW ANGLE            
%
% REMARKS:
% SVBCOS, the COS of the angle between the sv radius vector and the sun 
% radius vector, (the dot product of the above respective unit vectors), 
% SVBCOS is used to test against CNOON=cos(ANOON), if sv is entering 
% a noon maneuverer (i.e., SVBCOS > CNOON). The ANOON limit 
% (typically < 5.7 deg), is determined within the subroutine and depends 
% on the  sv yaw rate (YRATE)and GNSS (GPS or GLONASS). 
%
% ANIGHT , the shadow limit is the "input " in the subroutine 
% call statement, it can be hard coded as it is constant for a GNSS type,
% e.g., 
% 180.D0+-13.25D0 for GPS (IPRN >=  32) and 180.D0+-14.20D0 for 
% GLONASS (IPRN > 32), resp.). CNIGHT=cos(ANIGHT) is used for testing.
% When SVBCOS < CNIGHT (CNIGHT is negative and close to -1), 
% the SV enters, or it is in the shadow. 
%
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
%

%
%     MAXSAT - MAXIMUM NUMBER OF SATELLITES, CHANGED if REQUIRED
%      
      MAXSAT = 64;
	  
% MAX YAW RATES OF CURRENT&PAST BLOCK II/IIA's,(AVER'D 1996-2008 JPL SOLUTIONS)  
% CHANGE if REQUIRED OR INPUT if ESTIMATED 
% PRN                
      YRATE = [ 0.1211d0, 0.1339d0, 0.123d0, 0.1233d0, 0.1180d0, 0.1266d0,0.1269d0, ...             % 01 ~ 07
                0.1033d0, 0.1278d0, 0.0978d0, 0.200d0, 0.199d0, 0.200d0, 0.0815d0, 0.1303d0, ...    % 08 ~ 15
                0.0838d0, 0.1401d0, 0.1069d0, 0.098d0, 0.103d0, 0.1366d0, 0.1025d0, 0.1140d0, ...   % 16 ~ 23
                0.1089d0, 0.1001d0, 0.1227d0, 0.1194d0, 0.1260d0, 0.1228d0, 0.1165d0, 0.0969d0, ... % 24 ~ 31
                0.1140d0, 32*0.250d0];  % 32  33-64: GLONASS RATES (DILSSNER 2010)                                  
          
%  CHECK FOR BLOCK IIR AND FIX TO NOMINAL YAW RATE
      if (IPRN <= 32 && IBLK(IPRN) >= 4)
          YRATE(IPRN) = 0.2D0;
      end
% THE NEW GPS BLK Iif YAW RATES ( DILSSNER (2010) INSIDE GNSS)
      if (IPRN <= 32 && IBLK(IPRN) > 5)
          YRATE(IPRN) = 0.11D0;
      end
% De% 12, 2013
%  YBIAS=-0.5 FOR IIA (IBLK<4) PRN23 (SVN23 UP TO FEB2004,AFTER  IIR SVN60
%   AND NOT USED), Iif (IBLK=6)=-0.5, USED FOR SMALL NEG BETA NOON TURNS ONLY!
      YBIAS=0.0D0;
      if(IBLK(IPRN) <= 3)
          YBIAS= 0.5D0;
      end
      if(IPRN == 23 || IBLK(IPRN) == 6)
          YBIAS=-0.5D0;
      end
%
      IECLIPS=0;
      TWOHR = 7200.D0;
      HALFHR = 1800.D0;
      DTR = pi/180.D0;
% compute the noon beta angle limit (beta zero) FOR A NOON TURN from YRATEs
% & THE ACTUAL SAT ORBIT ANGLE RATE (MURATE) (~0.00836 FOR GPS; ~ 0.00888 GLNS)
      MURATE= sqrt((VSVC(1)^2+VSVC(2)^2+VSVC(3)^2) / (xsv(1)^2+xsv(2)^2+xsv(3)^2))/DTR;
      ANOON=atan(MURATE/YRATE(IPRN))/DTR;  													
% 31/03/2017 CNOON limit extension
      CNOON = cos( (ANOON + HALFHR * MURATE) * DTR);
      CNIGHT = cos(ANIGHT*DTR);
	  
	  
% 22 August 2016																
% Attention: CNIGHT limit expansion for shadow exit recovery
% Expand the CNIGHT limit for shadow exit recovery   
      if (IBLK(IPRN)<=3)
          CNIGHT = cos((ANIGHT + HALFHR * MURATE) * DTR);
      end
	  
      NOON = false;
      NIGHT = false;
      BETADG = beta/DTR - 90.d0;
	  
	  
      if (IPRN > 32 && abs(BETADG) < ANOON)
          % GLONASS NOON TURN MODE ACORDING TO DILSSNER 2010 
          YAWEND=75.D0;
          %  ITERATION FOR YAWEND OF THE GLONASS  NOON TURN
          for J=1:3
              YAWEND=abs(atan2(-tan(BETADG*DTR),sin(PI-DTR*MURATE*YAWEND/YRATE(IPRN)))/DTR ...
                          - atan2(-tan(BETADG*DTR),sin(PI+DTR*MURATE*YAWEND/YRATE(IPRN)))/DTR)/2.D0;
          end
% UPDATE ANOON, CNOON FOR NEW GLONASS NOON TURN LIMITS
          ANOON = MURATE*YAWEND/YRATE(IPRN);
          CNOON = cos(ANOON*DTR);
      end 
% BLK IIR'S
      if(IBLK(IPRN) == 4  ||  IBLK(IPRN) == 5)

%         CNIGHT=cos((ANOON+180.d0)*DTR); 										% Night turn analogous to Noon has been deactivated
         for J=1:3
             % BODY-X U VECTOR REVERSAL FOR IIR ONLY
             SANTXYZ(J)=-SANTXYZ(J);
         end
      end
 
      if (SVBCOS < CNIGHT)
          NIGHT=true;
      end
      if (SVBCOS > CNOON)
          NOON=true;
      end
	  
%     if SV IN NIGHT SHADOW OR NOON TURN DURING FORWARD PASS
%     STORE START AND END TIME OF YAW MANEUVRE (FOR THE BACKWARD RUN)

% init PHI
      PHI = pi/2.d0;
% YAW ANGLE
      YANGLE = cos((santxyz(1)*vsvc(1) +   ...
                    santxyz(2)*vsvc(2) + santxyz(3)*vsvc(3)) / sqrt(vsvc(1)^2+vsvc(2)^2 + vsvc(3)^2))/DTR;
 
	 
	 
% IIR YANGLE has the same sign as beta, II/IIA has the opposite sign
      if(BETADG < 0.d0 && IBLK(IPRN) >= 4 && IBLK(IPRN) <= 5)
          YANGLE = -YANGLE;
      end
      if(BETADG > 0.d0&&IBLK(IPRN) ~= 4 && IBLK(IPRN) ~= 5)
          YANGLE = -YANGLE;
      end

      if ( (NIGHT  ||  NOON))
          DET=sqrt((180.d0-acos(svbcos)/DTR)^2-BETADG^2);
          PHI = pi/2.d0;
          % Check if already after a midnight or noon
          if(NIGHT)
              if(IBLK(IPRN) == 4 || IBLK(IPRN) == 5)
                  if(abs(YANGLE) > 90.d0)
                      DET = -DET;
                  end
                  if(DET ~= 0.d0)
                      PHI=atan2(tan(BETADG*DTR),-sin(-DET*DTR))/DTR;
                  end
              else
                  % BLK IIA & GLONASS TOO !
                  if(abs(YANGLE) < 90.d0)
                      DET = -DET;
                  end
                  if(DET ~= 0.d0)
                      PHI=atan2(-tan(BETADG*DTR), sin(-DET*DTR))/DTR;
                  end
              end
          end
          if(NOON)
              % 29 July 2016
              % 180. changed to 180.D0
              DET=sqrt((acos(svbcos) * 180.D0 /pi)^2-BETADG^2);
              if(IBLK(IPRN) == 4 || IBLK(IPRN) == 5)
                  if (abs(YANGLE) < 90.d0)
                      DET=-DET;
                  end
                  if (DET ~= 0.d0)
                      PHI=atan2(tan(BETADG*DTR),-sin(PI-DET*DTR))/DTR;
                  end
              else
                  % BLK IIA & GLONASS !
                  if(Dabs(YANGLE) > 90.d0)
                      DET=-DET;
                  end
                  if(DET ~= 0.d0)
                      PHI=atan2(-tan(BETADG*DTR),sin(PI-DET*DTR))/DTR;
                  end
              end
          end

          % 28 July 2016
          MU = DET;
	  
          % 18 August 2016
          % Replace YANGLE with the PHI value (as computed in the directly previous lines of code)
          % YANGLE:	Angle of the dot product of the Body-X and Velocity unit vectors
          % PHI: Nominal yaw angle based on the fundamental equation atan2(-tan(beta),sin(MU))
          YANGLE = PHI;
	
          % ONLY FORWARD
          if (IDIR > 0)
              %  INITIALIZE ECLIPSE START AND TIME TAG ARRAYS  
              if ( NECLIPS(IPRN) == 0)
                  NECLIPS(IPRN)=NECLIPS(IPRN)+1;
                  ECLSTM(IPRN,NECLIPS(IPRN))=TTAG+DET/MURATE;
                  % IIR MIDNIGHT/NOON TURN or II/IIA NOON TURN START
                  % for IIR/GLONAS NIGHT (turn) only makes sense when BETADG < ANOON!
                  % For IIA it gets here only when NOON is true and that happens  only when BETADG < ANOON!
                  YAWEND=Atan(MURATE/YRATE(IPRN))/DTR;
                  if((IBLK(IPRN) > 3 && IBLK(IPRN)<=5 || NOON) && abs(BETADG) < YAWEND)
                      % GLONASS
                      if( IPRN > 32)
                          % GLONASS NOON TURN MODE ACORDING TO DILSSNER ET AL 2010 
                          ECLSTM(IPRN,NECLIPS(IPRN))= ECLSTM(IPRN,NECLIPS(IPRN))-ANOON/MURATE;
                          ECLETM(IPRN,NECLIPS(IPRN))= ECLSTM(IPRN,NECLIPS(IPRN))+2.D0*ANOON/MURATE;
                      else
                          % GPS IIA/IIR/Iif NOON OR IIR MIDNIGHT TURNs
                          ECLSTM(IPRN,NECLIPS(IPRN))= ECLSTM(IPRN,NECLIPS(IPRN))-abs(BETADG)*sqrt(ANOON/abs(BETADG)-1.d0)/MURATE;
                          ECLETM(IPRN,NECLIPS(IPRN))= ECLSTM(IPRN,NECLIPS(IPRN))+2*abs(BETADG)*sqrt(ANOON/abs(BETADG)-1.d0)/MURATE;
                      end
                  end
                  %  II/IIA SHADOW START & END TIMES
                  if((IBLK(IPRN)<=3 || IBLK(IPRN) > 5) && NIGHT)
                      ECLSTM(IPRN,NECLIPS(IPRN))= ECLSTM(IPRN,NECLIPS(IPRN))-sqrt((ANIGHT-180.d0)^2-BETADG^2)/MURATE;
                      ECLETM(IPRN,NECLIPS(IPRN))= ECLSTM(IPRN,NECLIPS(IPRN))+2.d0*sqrt((ANIGHT-180.d0)^2-BETADG^2)/MURATE;
                  end
              end
	  
              % UPDATE SV COSINE AND TIME TAG ARRAYS
              % (TO BE USED DURING BACKWARDS RUN)
              if ((NIGHT && SVBCOS < CNIGHT) || (NOON && SVBCOS > CNOON))
                  DTTAG = abs(TTAG-ECLSTM(IPRN,NECLIPS(IPRN)));
                  % ECLIPSE TIME IS MORE THAN 2 HOURS, THIS IS A NEW ECLIPSE!
                  if(DTTAG > TWOHR)
                      NECLIPS(IPRN)=NECLIPS(IPRN)+1;
                      ECLSTM(IPRN,NECLIPS(IPRN))=TTAG+DET/MURATE;
                      % IIR MIDNIGHT/NOON TURN  or II/IIA NOON TURN START
                      %                                  AND GLONASS NOON
                      if(IBLK(IPRN) > 3 && IBLK(IPRN) <= 5 || NOON)
                          % GLONASS
                          if(IPRN > 32)
                              % GLONASS NOON TURN MODE ACORDING TO DILSSNER ET AL 2010 
                              ECLSTM(IPRN,NECLIPS(IPRN))= ECLSTM(IPRN,NECLIPS(IPRN))-ANOON/MURATE;
                              ECLETM(IPRN,NECLIPS(IPRN))= ECLSTM(IPRN,NECLIPS(IPRN))+2.D0*ANOON/MURATE;
                          else
                              % GPS TURNS ONLY
                              ECLSTM(IPRN,NECLIPS(IPRN))= ECLSTM(IPRN,NECLIPS(IPRN)) - ...
                                                            abs(BETADG)*sqrt(ANOON/abs(BETADG)-1.d0)/MURATE;
                              ECLETM(IPRN,NECLIPS(IPRN))= ECLSTM(IPRN,NECLIPS(IPRN)) + ...
                                                            2*abs(BETADG)*sqrt(ANOON/abs(BETADG)-1.d0)/MURATE;
                          end
                      end
                      %     II/IIA SHADOW START & END TIMES
                      %   & GLONASS & Iif AS WELL !
                      if((IBLK(IPRN)<=3 || IBLK(IPRN) > 5)&&NIGHT)
                          ECLSTM(IPRN,NECLIPS(IPRN))= ECLSTM(IPRN,NECLIPS(IPRN))- ...
                                                             sqrt((ANIGHT-180.d0)^2-BETADG^2)/MURATE;
                          ECLETM(IPRN,NECLIPS(IPRN))= ECLSTM(IPRN,NECLIPS(IPRN))+ ...
                                                             2.d0*sqrt((ANIGHT-180.d0)^2-BETADG^2)/MURATE;
                      end
                  end
              end
              %  END OF FORWARD LOOP (IDIR = 1)
          end
      end
  

	  
%     BOTH FWD (IDIR= 1) OR BWD (IDIR=-1)
%     SET ECLIPSE FLAG (1=NIGHT SHADOW, 2=NOON TURN) 
      if (NECLIPS(IPRN)  ~=  0)
          % CHECK if IPRN IS ECLIPSING AND WHICH SEQ NO (I)
          I=0;
          for J=1:NECLIPS(IPRN)
              if(TTAG >= ECLSTM(IPRN,J) && TTAG<=(ECLETM(IPRN,J)+HALFHR))
                  I= J;
              end
          end
          % CURRENTLY NOT ECLIPSING (i=0)
          if(I == 0) 
              return;
          end
          if (TTAG  >=  ECLSTM(IPRN,I) && TTAG <= (ECLETM(IPRN,I)+HALFHR))
              % velocity & radius unit vectors V & R
              for J=1:3
                  V(J)=VSVC(J)/sqrt(VSVC(1)^2+VSVC(2)^2+VSVC(3)^2);
                  R(J)=XSV(J)/sqrt(XSV(1)^2+XSV(2)^2+XSV(3)^2); 
              end
              % ORBIT ANGLE MU AT ECLIPSE/TURN START
              DET= MURATE*(ECLETM(IPRN,I)-ECLSTM(IPRN,I))/2.d0;
              % De% 12, 2013 - start
              % YAWEND HERE - Iif SHADOW YAW RATE
              YAWEND=(atan2(-tan(BETADG*DTR), sin( DET*DTR))- ...
                             atan2(-tan(BETADG*DTR), sin(-DET*DTR)))/DTR/(ECLETM(IPRN,I)-ECLSTM(IPRN,I));
              % De% 12, 2013 - end

              if (SVBCOS < 0)
                  % SHADOW CROSSING
                  % BLK IIA/Iif SHADOW CROSSING
                  if(IPRN<=32 && (IBLK(IPRN)<=3 || IBLK(IPRN) > 5))
                      if(TTAG<=ECLETM(IPRN,I))
                          % IIA NIGHT TURN
                          if(IBLK(IPRN) <= 3) 
                              PHI=atan2(-tan(BETADG*DTR), sin(-DET*DTR))/DTR ...
                                    +sign(YRATE(IPRN),BETADG)*(TTAG-ECLSTM(IPRN,I));     % Rev.12/09/2016: Case: IIA Beta<0
                          end
	                       
                          % Iif NIGHT TURN (DILSSNER  2010)
                          if(IBLK(IPRN) > 5)
                              PHI=atan2(-tan(BETADG*DTR), sin(-DET*DTR))/DTR ...
                                                          +YAWEND*(TTAG-ECLSTM(IPRN,I));
                                  % De% 12, 2013 - end
                          end
                      else
                          % ** WARNING
                          % IIA/Iif SHADOW EXIT RECOVERY: USING THE IIA DATA  DURING
                          % THE IIA RECOVERY (UP TO 30 MIN) IS NOT RECOMMENDED!
                          % ** WARNING
                          % GPS IIA  AT SHADOW EXIT
                          if (IBLK(IPRN) <= 3)
                              PHI=atan2(-tan(BETADG*DTR), sin(-DET*DTR))/DTR ...
                                                 +sign(YRATE(IPRN),BETADG)*(ECLETM(IPRN,I)-ECLSTM(IPRN,I));    % Rev.12/09/2016: Case: IIA Beta<0 
                          end
	                      
                          % GPS Iif AT SHADOW EXIT
                          % De% 12, 2013
                          % NO NEED FOR Iif RECOVERY ALREADY AT THE EXIT YAW!
                          if(IBLK(IPRN) > 5)
                              return;
                          end
                          % YAWEND- HERE THE ACTUAL YAW DifFERENCE  AT THE SHADOW EXIT
                          YAWEND= YANGLE- PHI;   
                          % 23 August 2016	% Rev.
                          % Deactivated to avoid error (opposite) SIGN
                          PHI=PHI+sign(YRATE(IPRN),YAWEND)*(TTAG-ECLETM(IPRN,I));

                          % SANTX- THE CURRENT ANGLE DifF, CONSISTENT WITH YAWEND
                          SANTX= YANGLE-PHI;

                          % 23 August 2016		% Rev.
                          % Deactivated to avoid error (opposite) SIGN				   
                          % STOP! THE NOMINAL YAW (YANGLE) REACHED!
                          if(abs(SANTX) > abs(YAWEND))
                              return;
                          end
                          if(YAWEND ~= 0.D0&&((SANTX)/YAWEND) < 0.D0)
                              return;
                          end
				              
                          % 19 September 2016	 % Rev.
                          % Deactivated during initial epochs of IIA exit recovery where may abs(Phi)>180				   
                          			  
                          % SET PHI <-180,+180>
                          %                   PHI= DMOD(PHI, 360.D0)
                          %                   if(abs(PHI) > 180.D0) PHI= PHI-360.D0*PHI/abs(PHI)
                      end
                  end

                  % GLONASS
                  if(IPRN > 32)
                      % GLONASS   NIGHT TURN (DILSSNER AT AL 2010 )
                      if(TTAG > ECLETM(IPRN,I))
                          return;
                      end
                      YAWEND=YRATE(IPRN);
                      PHI=atan2(-tan(BETADG*DTR), sin(-DET*DTR))/DTR ...
                                     +sign(YAWEND,BETADG)*(TTAG-ECLSTM(IPRN,I));
                      % YAWEND -YAW ANGLE AT THE (GLONASS) SHADOW EXIT
                      YAWEND=atan2(-tan(BETADG*DTR), sin( DET*DTR))/DTR;
                      if((YAWEND/PHI) >= 1.d0 || (PHI/YAWEND) < 0.d0)
                          PHI = YAWEND;
                      end
                  end

                  % GPS BLK Iif NIGHT YAW RATE(DILSSNER 2010):
                  if(IBLK(IPRN) > 3 && IBLK(IPRN) <= 5)
                      % BLK II R SHADOW (MIDNIGHT TURN) CROSSING
                      PHI=atan2(tan(BETADG*DTR),-sin(-DET*DTR))/DTR ...
                                 +sign(YRATE(IPRN),BETADG)*(TTAG-ECLSTM(IPRN,I));
                      if((PHI/YANGLE) >= 1.d0 || (PHI/YANGLE) < 0.d0)
                         return;				% Deactivate condition for return (GOTO 1)
                      end
                      % & BETADG, ECLETM(IPRN,I),I
                      IECLIPS=1;
                  else
                      % NOON TURNS 
                      PHI=atan2(-tan(BETADG*DTR),sin(PI-DET*DTR))/DTR ...
                                  -sign(YRATE(IPRN),BETADG)*(TTAG-ECLSTM(IPRN,I));
                      % De% 12, 2013 -start
                      % SMALL NEGATIVE BETA Iif OR SMALL POSIT. IIA NOON TURN PROBLEM
                      % Jan 24, 2014
                      if (IPRN <= 32 && (BETADG*sign(1.D0,YBIAS)) <= 0.1D0 && ...   % Rev.19/09/2016: Beta sign change 
                                                                  (BETADG*YBIAS) > 0.0D0)
                          PHI=atan2(-tan(BETADG*DTR),sin(PI-DET*DTR))/DTR ...
                                            +sign(YRATE(IPRN),YBIAS )*(TTAG-ECLSTM(IPRN,I));
                      end
                      
                      % De% 12, 2013 - end
                      if (IBLK(IPRN) > 3 && IBLK(IPRN)<=5)
                          % BLK IIR NOON TURNS ONLY
                          PHI=atan2(tan(BETADG*DTR),-sin(PI-DET*DTR))/DTR ...
                                       -sign(YRATE(IPRN),BETADG)*(TTAG-ECLSTM(IPRN,I));
                          % IIR END TURN CHECK
                          if ((YANGLE/PHI) >= 1.d0 || (PHI/YANGLE) < 0.d0)
                              return;
                          end
                      else
                          % GLONASS END TURN CHECK
                          if(IPRN  > 32 && TTAG > ECLETM(IPRN,I))
                              return;
                          end
                          % IIA OR Iif END TURN CHECK
                          % De% 12, 2013 -start
                          if (IPRN <= 32 && BETADG*sign(1.D0,YBIAS)<= 0.9D0 && ...
                                   BETADG*YBIAS > 0.0D0 && ...
                                   (((PHI-sign(1.d0,YBIAS)*360.D0)/YANGLE)<=1.d0 || ... 
                                   ((PHI-sign(1.D0,YBIAS)*360.D0)/YANGLE) < 0.d0))
                              return;
                          end
                          
                          if(IPRN <= 32 && (BETADG*sign(1.D0,YBIAS) > 0.9D0 || BETADG*YBIAS<=0.0D0) && ...
                                     ((PHI/YANGLE) >= 1.d0 || (PHI/YANGLE) < 0.d0))
                              return;
                          end
                          % Rev.03/04/2017: GPS IIA/Iif Noon turn end check  % Rev.03/04/2017 			 	 
                          if (IBLK(IPRN)<=3  ||  IBLK(IPRN)==6)
                              if (abs(YANGLE) > 170.D0 && abs(YANGLE/PHI) < 1.0D0)
                                  return;
                              end
                          end
                          
                          IECLIPS = 2;
                      end
                      
                      % ROTATE X-VECTOR TO ECLIPSING YAW ANGLE PHI 
                      % ECLIPSING (II/IIA) NOT TO BE USED  A HALF HR AFTER SHADOW !
                      SANTX=(cos((PHI-YANGLE)*DTR)*(V(2)-V(3)*R(2)/R(3))-cos(PHI*DTR)*  ...
                             (SANTXYZ(2)-SANTXYZ(3)*R(2)/R(3)))/(SANTXYZ(1)*V(2)-SANTXYZ(2)*v(1)+ ...
                             ((SANTXYZ(2)*V(3)-SANTXYZ(3)*V(2))*R(1)+(SANTXYZ(3)*V(1)-SANTXYZ(1)*V(3))*R(2))/R(3));
                      SANTY = (cos(PHI*DTR) - (V(1)-V(3)*R(1)/R(3))*SANTX)/(V(2)-V(3)*R(2)/R(3));
                      % THE BODY-X UNIT VECTOR ROTATED BY (PHI-YANGLE) RETURNED
                      SANTXYZ = [SANTX, SANTY, (-R(1)*SANTX-R(2)*SANTY)/R(3)];
                  end
              end
              
              % 27 June 2016
              % The output argument "yaw_angle" has been added for providing the computed 
              % values of the nominal yaw angle (YANGLE) and the eclipsing yaw angle (PHI)            
              yaw_angle(1) = YANGLE;
              yaw_angle(2) = PHI;

              % 28 July 2016
              % 'DET' and 'MURATE' have been included to the output arguments
          end
      end
end