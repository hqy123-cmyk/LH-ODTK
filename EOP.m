function [EOP_cr, CRS2TRS, TRS2CRS, d_CRS2TRS, d_TRS2CRS] = EOP(mjd)
% ----------------------------------------------------------------------
% Purpose:
%  IERS Earth Orientation Parameters (EOP) data reading and processing.
%  Corrections computing due to tidal variations
%  Computation of the ICRF-ITRF transformation matrix (direct/inverse) 
%  and its derivatives.
%  Earth Orientation Matrix computation is based on the EOP data corrected 
%  due to tidal variations.
% ----------------------------------------------------------------------
% Input arguments:
% - mjd:			Modified Julian Day number of the epoch (in TT scale)
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
%   dX, dY:			Corrections to the Precession-Nutation model (arcsec)  岁差-章动模型的改正
% - CRS2TRS:		GCRS to ITRS transformation matrix  天球参考系(GCRS)至国际地球参考系(ITRS)的旋转矩阵
% - TRS2CRS:		ITRS to GCRS transformation matrix  国际地球参考系(ITRS)至天球参考系(GCRS)的旋转矩阵
% - d_CRS2TRS:		Derivative of GCRS to ITRS transformation matrix   天球参考系(GCRS)至国际地球参考系(ITRS)的旋转矩阵的偏导数
% - d_TRS2CRS:		Derivative of ITRS to GCRS transformation matrix   国际地球参考系(ITRS)至天球参考系(GCRS)的旋转矩阵的偏导数
% ----------------------------------------------------------------------

   global n_interp iau_pn_model;
   
   mjd_TT = mjd;
   
   % EOP data
   EOP_day_glb = eop_data(mjd_TT, n_interp);
   
   % EOP by IERS data: EOP reading, interpolation and corrections
   EOP_cr = eop_cor(mjd_TT, EOP_day_glb);
   
   %SCM 20190606 Removed the nutation model correction terms from EOP array,
   % These terms were adding noise to the orbit fits to IGS products.
   % Add flag to turn on/off in the future
   EOP_cr(6) = 0.0;
   EOP_cr(7) = 0.0;
   
    
   % -----------------------------
   % Tranformation: ITRF to ICRF
   % -----------------------------
   % ICRF-ITRF transformation matrix (including derivatives)
   [CRS2TRS, TRS2CRS, d_CRS2TRS, d_TRS2CRS] = crs_trs(mjd_TT, EOP_cr, iau_pn_model);
   
end