function [EOP_days] = eop_data(mjd, n_interp)
% ----------------------------------------------------------------------
% Subroutine:  eop_data.m
% ----------------------------------------------------------------------
% Purpose:
%  IERS Earth Orientation Parameters (EOP) data reading and processing.
%  Corrections computing due to tidal variations
% ----------------------------------------------------------------------
% Input arguments:
% - mjd:			Modified Julian Day number of the epoch (in TT)
% - EOP_fname:		EOP data file name  e.g. 'eopc04_IAU2000.62-now'
% - EOP_sol:		EOP solution type (see Note 1)
% - n_interp:		number of points to be used for EOP data interpolation 
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
   
   global eopdata;

   % -------------------------------
   % Time Systems transformation	
   % -------------------------------
   [mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC] = time_TT(mjd);
   
   % UTC
   mjd_UTC_day = round(mjd_UTC);
   
   % Epochs center definition
   if (n_interp/2 - round(n_interp/2) == 0)
       n_cent = n_interp/2;
   else
       n_cent = round(n_interp/2) + 1;
   end
   
   mjd_day_int = mjd_UTC_day;
   
   % --------------------------
   % EOP data reading          
   % --------------------------
   kk = 1;
   for i = 1:n_interp
       MJD_i = mjd_day_int - n_cent + i;
       for j = kk:size(eopdata,2)
           if (MJD_i == eopdata(4,j))
               EOP_days(i,1) = MJD_i;
               EOP_days(i,2) = eopdata(5,j);
               EOP_days(i,3) = eopdata(6,j);
               EOP_days(i,4) = eopdata(7,j);
               EOP_days(i,5) = eopdata(8,j);
               EOP_days(i,6) = eopdata(11,j);
               EOP_days(i,7) = eopdata(12,j);
               kk = j;
           end
       end
   end

end