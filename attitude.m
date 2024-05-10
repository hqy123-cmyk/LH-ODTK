function [eclipsf, beta, Mangle, Yangle, eBX_nom, eBX_ecl] = attitude(mjd, rsat, vsat, rSun, PRNsat, BLKsat)
% ----------------------------------------------------------------------
% Purpose:
%  GNSS yaw-attitude modelling during nominal and eclipse seasons.
%  This subroutine "attitude.f03" has been obtained by the partial conversion
% of the GNSS Yaw-attitude Tool (GYT) standalone program to subroutine  
% ----------------------------------------------------------------------
% Input arguments:
% ----------------------------------------------------------------------
% - mjd:			Modified Julian Day number (including the fraction of the day) in Terrestrial Time (TT)
% - rsat:			Satellite Position vector (m)   in inertial frame (ICRF)
% - vsat:			Satellite Velocity vector (m/s) in inertial frame (ICRF)
% - rSun:			Sun Position vector (m)   in inertial frame (ICRF)
% - PRNsat:		    GNSS PRN number (Constellation ID + Number) e.g. G03
% - BLKsat:			GNSS satellites block type
%
% Output arguments:
% - eclipsf:		Satellite attitude status flag:
%					0 = Nominal
%					1 = Non-nominal, midnight
%					2 = Non-nominal, noon
% - beta:			Sun angle with the orbital plane in degrees (see Note 1)
% - Yangle: 		Yaw angle array (in degrees) with 2 values: Yaw nominal and Yaw modelled
%					During nominal periods, the two values are equal
% - Mangle:			Orbit angle between satellite and orbit midnight (degrees)
% - eBX_nom:		Body-X axis unit vector based on noninal yaw-attitude
% - eBX_ecl:		Body-X axis unit vector based on modelled yaw-attitude during eclipse seasons (Set to eBX_nom when out of eclipse season) 
%
% ----------------------------------------------------------------------
% Note 1:
%  Beta angle is computed based on the input vectors in inertial system
%  Sign conventions for beta angle:
%   Postiive: above orbital plane
%   Negative: below orbital plane  	
% ----------------------------------------------------------------------
% Author :	Dr. Thomas Papanikolaou
%			Geoscience Australia, Frontier-SI
% Created:	July 2016
% ----------------------------------------------------------------------
% Last modified
% - Dr. Thomas Papanikolaou, 21 November 2018:
%	GYT program upgrade from Fortran 90 to Fortran 2003 for calling the modified 
%   subroutines (upgraded to F03) of the POD code package 
% - Dr. Thomas Papanikolaou, 25 March 2019:
%   The major part of the GYT program has been converted into subroutine 
%   for the integration to the Precise Orbit Determination (POD) package  
% ----------------------------------------------------------------------

  % ----------------------------------------------------------------------
  % GNSS Satellite Block Type
  % ----------------------------------------------------------------------
  % BLK_TYP :: Global variable in mdl_param
  % GPS case: Satellite Block ID:        1=I, 2=II, 3=IIA, IIR=(4, 5), IIF=6
  % BLKTYP = BLKsat
  
  satblk = 6;
  if(BLKsat == "GPS-I")
  	  satblk = 1;
  elseif (BLKsat == "GPS-II")
  	  satblk = 2;
  elseif (BLKsat=="GPS-IIA")
      satblk = 3;
  elseif (BLKsat=="GPS-IIR")
      satblk = 4;
  elseif (BLKsat=="GPS-IIR-A")
      satblk = 5;
  elseif (BLKsat=="GPS-IIR-B")
      satblk = 5;
  elseif (BLKsat=="GPS-IIR-M")
      satblk = 5;
  elseif (BLKsat=="GPS-IIF")
      satblk = 6;
  end
  % ----------------------------------------------------------------------
  % Beidou case: 'IGSO', 'MEO'
  % 1. BDSorbtype = 'IGSO'
  % 2. BDSorbtype = 'MEO'
  % 3. BDSorbtype = 'IGSO'
  BDSorbtype = "MEO";
  if (BLKsat=="BDS-2G" || BLKsat == "BDS-3G")
      BDSorbtype = "GEO";
  end
  if (BLKsat == "BDS-2I" || BLKsat == "BDS-3I" || ...
     BLKsat == "BDS-3SI-SECM" || BLKsat == "BDS-3SI-CAST")
      BDSorbtype = "IGSO";
  end
  if (BLKsat == "BDS-2M" || BLKsat == "BDS-3M" || ...
     BLKsat == "BDS-3M-SECM" || BLKsat == "BDS-3M-CAST")
      BDSorbtype = "MEO";
  end
  % ----------------------------------------------------------------------
  
  
  % ----------------------------------------------------------------------
  % Satellite Position & Velocity vector in ICRF
  r_CRS = rsat;
  v_CRS = vsat;
  % Sun position vector in ICRF
  r_sun_crs = rSun;
  % ----------------------------------------------------------------------
  
  
  % ----------------------------------------------------------------------
  % GNSS Yaw-attitude Tool (GYT) program configuration
  % ----------------------------------------------------------------------
  % Start of INPUT Configuration parameters:
  % ----------------------------------------------------------------------
  
  % ----------------------------------------------------------------------
  % Yaw attitude model/method
  % ----------------------------------------------------------------------	    
  % 1. Yaw attitude modelling applied based on the models adopted for each GNSS constellation
  %    Set yaw_mod variable to the GNSSid (GNSS constellation id letter)
  %    GNSS constellation id letters
  %    G: GPS
  %    R: GLONASS
  %    E: Galileo
  %	   C: BDS (BeiDou)
  % ----------------------------------------------------------------------
  % 2. Dynamic Yaw-steering method 
  %    The yaw_mod variable should set to 'D' for special cases of  
  %    Galileo and Beidou satellites that apply the dynamic yaw-steering method
  % ----------------------------------------------------------------------
  %yaw_mod = 'G'
  %yaw_mod = 'R'
  %yaw_mod = 'E'
  %yaw_mod = 'C'
  %yaw_mod = 'D'
  
  read (PRNsat, FMT='(A1)' , IOSTAT=ios) GNSSid;
  yaw_mod = GNSSid;
  %print *,"yaw_mod", yaw_mod
  % ----------------------------------------------------------------------
  
  
  % ----------------------------------------------------------------------
  % Beidou case only:
  % ----------------------------------------------------------------------
  % Approach for the beta angle computation at the latest epoch that the orbital angle M is equal to 90 degrees
  % ----------------------------------------------------------------------
  % BetaP = 1 : Orbit backward prediction based on Keplerian elements  
  % BetaP = 2 : Orbit backward computation based on numerical interpolation of the precise orbit data sp3 
  BetaP = 1;
  % ----------------------------------------------------------------------
  
  
  % ----------------------------------------------------------------------	  
  % Satellite Body-fixed frame definition type
  % ----------------------------------------------------------------------	  
  % satbf = 1 : Body-fixed frame according to the IGS Conventions; Cases: GPS Block II,IIA,IIF, GLONASS, BeiDou  
  % satbf = 2 : Body-fixed frame X,Y axes reversal; Cases: Galileo, GPS Block IIR 
  satbf = 1;
  % ----------------------------------------------------------------------	
  
  
  % ----------------------------------------------------------------------
  % Dynamic yaw-steering method
  % ----------------------------------------------------------------------
  % yawdyn.f90 subroutine input arguments to be configured
  % ----------------------------------------------------------------------	  
  % Beta angle condition limit for the "dynamic yaw steering" method
  beta0 = 2.0;
  % Design parameter
  dparam = 258;
  % ----------------------------------------------------------------------	  
  
  
  % ----------------------------------------------------------------------
  % End of INPUT Configuration
  % ----------------------------------------------------------------------
  
  
  
  % ----------------------------------------------------------------------
  % Satellite PRN number 
  % ----------------------------------------------------------------------
  % Consisted by the GNSS constellation id letter (GNSSid) and the number (PRN_num)  
  % ----------------------------------------------------------------------
  % GNSSid:		GNSS constellation id letter
  % 				G: GPS
  % 				R: GLONASS
  % 				E: Galileo
  % 				C: BDS (BeiDou)
  % 				J: QZSS
  % 				S: SBAS
  % ----------------------------------------------------------------------
  % PRN_num:		PRN numbering as adopted in the sp3 format (Numerical value after the GNSS constellation id letter)
  % ----------------------------------------------------------------------
  %fmt_line = '(A1,I2.2)'
  READ (PRNsat, FMT='(A1,I2.2)' , IOSTAT=ios) GNSSid, PRN_num;
  
  % ----------------------------------------------------------------------
  % PRN_eclips:	PRN number as adopted in the eclips.f subroutine numbering
  %			 	GPS: 		1-32 
  % 			GLONASS:	33-64
  %				Galileo:	65-100
  %				BeiDou:		101-136
  % ----------------------------------------------------------------------
  if (GNSSid == 'G')
      PRN_eclips = PRN_num;
  elseif  (GNSSid == 'R')
      PRN_eclips = PRN_num + 32;
  elseif  (GNSSid == 'E')
      PRN_eclips = PRN_num + 64;
  elseif  (GNSSid == 'C')
      PRN_eclips = PRN_num + 100;
  end
  % ----------------------------------------------------------------------
  
  
  % ----------------------------------------------------------------------
  % Beta angle
  beta = beta_angle(r_CRS, v_CRS, r_sun_crs);
  % ----------------------------------------------------------------------
  
  
  % ----------------------------------------------------------------------
  % Time System of input epoch:
  % ----------------------------------------------------------------------
  % 1. TT
  % 2. GPS time
  % 3. UTC
  % 4. TAI
  time_in = TT_time;
  % ----------------------------------------------------------------------
  
  % ----------------------------------------------------------------------
  % "Time Systems" transformation											 
  % ----------------------------------------------------------------------
  if (time_in == TT_time)
      [mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC] = time_TT(mjd);
  elseif (time_in == GPS_time)
      [mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC] = time_GPS(mjd);
  elseif (time_in == UTC_time)
      [mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC] = time_UTC(mjd);
  elseif (time_in == TAI_time)
      [mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC] = time_TAI(mjd);
  end 
  % ----------------------------------------------------------------------
  
  
  
  % ----------------------------------------------------------------------
  % Yaw angle computation based on the configured attitude model
  % ----------------------------------------------------------------------	    
  
  % -----------------------------------	    
  % GPS and GLONASS yaw-attitude model
  % -----------------------------------	    
  if (yaw_mod == 'G'  ||  yaw_mod == 'R')
      % Eclipsing attitude modelling is computed based on the modified version of the eclips.f
      % orbdir: Direction of processing (1=FORWARD, -1=BACKWARD); eclips.f input argument
      orbdir = 1; 
      
      % Nominal and Eclipsing Yaw-attitude
      [eclipsf, eBX_nom, eBX_ecl, Yangle, Mangle, Mangle_e, Mrate_e, Ynom_e] = yaw_attitude(mjd_GPS, r_CRS, ...
                                                            v_CRS, r_sun_crs, beta, PRN_eclips, satblk, orbdir, satbf);
      				   
      if (eclipsf == 0)
          Yangle(2) = Yangle(1);
      end
  % ----------------------	    
  % Galileo attitude law
  % ----------------------	    
  elseif (yaw_mod == 'E')
      % Yaw-steering as provided by the European GNSS Service Centre 
      [eclipsf, eBX_nom, eBX_ecl, Yangle, Mangle] = yaw_gal(mjd_GPS, r_CRS, v_CRS, r_sun_crs, ....
                                                                         beta, PRNsat, BLKsat, satbf);
  % ---------------------- 
  % BeiDou attitude law
  % ----------------------    
  elseif (yaw_mod == 'C')
      % Beta angle based on interpolation is cancelled
      NPint = 0;
      [beta, beta_t0, eclipsf, eBX_nom, eBX_ecl, Yangle, Mangle] = yaw_bds(mjd_GPS, r_CRS, v_CRS, r_sun_crs, ...
                                                                                    BDSorbtype, satbf, BetaP, NPint);
  % -----------------------------	    
  % Dynamic Yaw-steering method
  % -----------------------------	    
  elseif (yaw_mod == 'D')
      [beta, eclipsf, eBX_nom, eBX_ecl, Yangle, Mangle] = yawdyn(mjd_GPS, r_CRS, v_CRS, r_sun_crs, satbf, beta0, dparam);    
  end
end