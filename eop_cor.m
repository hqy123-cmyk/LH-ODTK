function [EOP_cr] = eop_cor(mjd, EOP_days)
% ----------------------------------------------------------------------
% Subroutine:  eop_cor.f90
% ----------------------------------------------------------------------
% Purpose:
%  IERS Earth Orientation Parameters (EOP) data reading and processing.
%  Corrections computing due to tidal variations
% ----------------------------------------------------------------------
% Input arguments:
% - mjd:			Modified Julian Day number of the epoch (in TT)
% - EOP_days:		EOP data array of the days (data points aplied for interpolation)
%
% Output arguments:
% - EOP_cr:		Array of the EOP data after applying the tidal corrections
%   EOP_cr = [MJD xp_cor yp_cor UT1-UTC_cor LOD_cor dX dY]  1x6 matrix
%   MJD:     		MJD of the input epoch (in TT)
%   xp_cor,yp_cor:	Polar motion coordinates (arcsec) at input epoch considering 
%					corrections due to tidal variations (ocean tidal and libration effects)  	
%   UT1UTC_cor:		UT1-UTC difference (sec) at input epoch considering corrections
%					due to tidal variations	(ocean tidal and libration effects)  
%   LOD:			Length of Day (LOD)
%   dX, dY:			Corrections to the Precession-Nutation model (arcsec)
% ----------------------------------------------------------------------

   % ------------------------------
   % Time Systems transformation
   % ------------------------------
   [mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC] = time_TT(mjd);
   
   % UTC
   mjd_UTC_day = round(mjd_UTC);
   
   % EOP data array
   sz1_EOP = size(EOP_days,1);
   sz2_EOP = size(EOP_days,2);
   
   % Ensure LOD, dX_eop, dY_eop are set to something sensible
   LOD = 0.0;
   dX_eop = 0.0;
   dY_eop = 0.0;
   
   for i = 1:sz1_EOP
       MJDint_ar(i) = EOP_days(i,1);
       % xp,yp (arcsec)
       xint_ar(i) = EOP_days(i,2);
       yint_ar(i) = EOP_days(i,3);

       % UT1-UTC (sec)
       UT1int_ar(i) = EOP_days(i,4);
   
       if (mjd_UTC_day == EOP_days(i,1))
           LOD = EOP_days(i,5);
           dX_eop = EOP_days(i,6);
           dY_eop = EOP_days(i,7);
       end
   end
   
   % Diurnal and semi-diurnal tidal corrections to EOP data by IERS
   % "Polar motion" and "UT1" corrections due to ocean tidal and libration effects at interpolation epoch
   [x_int, y_int, ut1_int] = INTERP_iers(MJDint_ar, xint_ar, yint_ar, UT1int_ar, mjd_UTC);
   
   % xp,yp corrected (arcsec) 											
   xp_cor = x_int;
   yp_cor = y_int;
   
   % UT1-UTC corrected (sec) 													
   UT1UTC_cor = ut1_int;
   
   EOP_cr(1) = mjd;
   EOP_cr(2) = xp_cor;
   EOP_cr(3) = yp_cor;
   EOP_cr(4) = UT1UTC_cor;
   EOP_cr(5) = LOD;
   EOP_cr(6) = dX_eop;
   EOP_cr(7) = dY_eop;
end