function GMST = gmst_iers(mjd, ut1_utc)
% ----------------------------------------------------------------------
% Purpose:
%  Greenwich Mean Sidereal Time (GMST) computation based on the IERS 
%  Convetions 2010 and IAU 2006 resolutions
% ----------------------------------------------------------------------
% Input arguments
% - mjd:			Modified Julain Day (MJD) number in Terrestrial Time (TT) 
% 					scale including the fraction of the day 
% - ut1_utc:   		Time difference between UT1 and UTC (seconds)
%
% Output arguments:
% - GMST :			Greenwich mean sidereal time (radians)
% ----------------------------------------------------------------------
% Author :	Dr. Thomas Papanikolaou, Cooperative Research Centre for Spatial Information, Australia
% Created:	27 October 2017
% ----------------------------------------------------------------------

  % Time Systems transformation											 
  [mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC] = time_TT(mjd);
  % ----------------------------------------------------------------------
  
  % ----------------------------------------------------------------------
  % TT
        TT1 = 2400000.5D0;
        TT2 = mjd_TT;
  % ----------------------------------------------------------------------
  % UTC
        mjd_UTC_day = round(mjd_UTC);
  % ----------------------------------------------------------------------
  % UT1
        mjd_UT1 = mjd_UTC + UT1_UTC / 86400.D0;
        TT1_UT1 = 2400000.5D0;
        TT2_UT1 = mjd_UT1;  
  % ----------------------------------------------------------------------
  
  
  % ----------------------------------------------------------------------
  % IAU SOFA subroutines
  % Model consistent with IAU 2000 resolutions	  
  GMST00 = iau_GMST00 (TT1_UT1, TT2_UT1, TT1, TT2);
  
  % IERS Conventions 2010 Equation 5.32
  % Model consistent with IAU 2006 resolutions	  
  GMST06 = iau_GMST06 (TT1_UT1, TT2_UT1, TT1, TT2);
  % ----------------------------------------------------------------------
  
  
  % ----------------------------------------------------------------------
  % Greenwich Mean Sidereal Time (GMST) in radians
  %GMST = GMST00
  GMST = GMST06;
end