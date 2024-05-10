function [CRS2TRS, TRS2CRS, d_CRS2TRS, d_TRS2CRS] = crs_trs(mjd, EOP_ar, iau_model)
% ---------------------------------------------------------------------------------------
% Purpose:
%  ICRF-ITRF transformation matrix (direct/inverse) for position/velocity vectors
%  Earth Orientation Matrix based on the EOP data corrected due to tidal variations
%
% Remark:
%  The matrices required for the ICRF-ITRF direct and inverse transformation
%  are computed for the position and velocity vectors individually.
% ---------------------------------------------------------------------------------------
% Input arguments:
% - mjd:			Modified Julian Day number of the epoch (in TT scale)
% - iau_model:		Precession-Nutation model by International Astronomical Union (IAU)
%					IAU = 2 refers to the IAU 2000A model
%					IAU = 3 refers to the IAU 2006/2000A model
% - EOP_ar:			Array of the EOP values for the input epoch
%   				eop = [MJD xp yp UT1_UTC LOD dX dY]  1x6 matrix
%   				MJD:     MJD of the input epoch (in TT) 
%   				x,y:     Polar motion coordinates (arcsec) 
%   				UT1_UTC: Difference between UT1 and UTC (sec)
%   				LOD_cor: Length of Day (LOD)
%					dX,dY:   Corrections to Precession-Nutation model (arcsec)
%
% Output arguments:
% - CRS2TRS:		GCRS to ITRS transformation matrix
% - TRS2CRS:		ITRS to GCRS transformation matrix 
% - d_CRS2TRS:		Derivative of GCRS to ITRS transformation matrix
% - d_TRS2CRS:		Derivative of ITRS to GCRS transformation matrix
% -----------------------------------------------------------------------------------------

  arcsec2rad = pi / (3600.0D0 * 180.0D0);

  % EOP values
  xp = EOP_ar(2);
  yp = EOP_ar(3);
  UT1_UTC = EOP_ar(4);
  dX_eop = EOP_ar(6);
  dY_eop = EOP_ar(7);

  % ----------------------------------------------------------------------
  % Time Systems transformation											 
  % ----------------------------------------------------------------------
  [mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC] = time_TT(mjd);
  % ----------------------------------------------------------------------
  % TT
  TT1 = 2400000.5;
  TT2 = mjd_TT;
  % ----------------------------------------------------------------------
  % TAI
  TAI1 = 2400000.5;
  TAI2 = mjd_TAI;
  % ----------------------------------------------------------------------
  % UTC
  mjd_UTC_day = round(mjd_UTC);
  % ----------------------------------------------------------------------
  % UT1
  mjd_UT1 = mjd_UTC + UT1_UTC / 86400;
  TT1_UT1 = 2400000.5;
  TT2_UT1 = mjd_UT1;
  % ----------------------------------------------------------------------

  % ----------------------------------------------------------------------
  % CIO based transformation
  % ----------------------------------------------------------------------
  
  
  % ----------------------------------------------------------------------
  % Q(t)
  % ----------------------------------------------------------------------
  %-----------------------------------------------------------------------
  % variable initialisation (X_iau, Y_iau, s_iau)
  X_iau = 0.d0;
  Y_iau = 0.d0;
  s_iau = 0.d0;

  % ----------------------------------------------------------------------
  % Precession-Nutation model:  X, Y (in radians)
  % ----------------------------------------------------------------------
  % 1. IAU 2000A
  if (iau_model == "IAU2000")
      [X_iau00A, Y_iau00A, s_iau00A] = iau_XYS00A(TT1, TT2);
      X_iau = X_iau00A;
      Y_iau = Y_iau00A;
      s_iau = s_iau00A;
  % 2. IAU 2006/2000A
  elseif (iau_model == "IAU2006")
      [X_iau06, Y_iau06] = iau_XY06(TT1, TT2);

      % CIO locator s (in radians) | Function iau_S06
      s_iau06 = iau_S06(TT1, TT2, X_iau06, Y_iau06);

      %or   CALL XYS06A
      [X_iau06A, Y_iau06A, s_iau06A] = iau_XYS06A(TT1, TT2);
      X_iau = X_iau06A;
      Y_iau = Y_iau06A;
      s_iau = s_iau06A;
  end
  
  % EOP: dX,dY (arcsec to radians)
  X_pn = X_iau + dX_eop * arcsec2rad;
  Y_pn = Y_iau + dY_eop * arcsec2rad;
  
  % Q(t)^-1
  RC2I = iau_C2IXYS(X_pn, Y_pn, s_iau);
  RC2I_T = RC2I';
  
  % ----------------------------------------------------------------------
  % R(t)
  % ----------------------------------------------------------------------
  %- ERA00: era (in radians)
  era = iau_ERA00(TT1_UT1, TT2_UT1);
  % ----------------------------------------------------------------------
  % Rotation Matrices initialization: set them to singular I3x3
  % R3(ERA)
  R_era(1,1) = 1.0;
  R_era(1,2) = 0.0;
  R_era(1,3) = 0.0;
  R_era(2,1) = 0.0;
  R_era(2,2) = 1.0;
  R_era(2,3) = 0.0;
  R_era(3,1) = 0.0;
  R_era(3,2) = 0.0;
  R_era(3,3) = 1.0;

  % R3(-ERA)
  R_era_n = R_era;

  % R(t)^T = R3(ERA)	  
  R_era = iau_RZ (era, R_era);

  % R(t) = R3(-ERA)
  era_n = -1.0D0 * era;
  R_era_n = iau_RZ(era_n, R_era_n);
  % ----------------------------------------------------------------------
  
  % ----------------------------------------------------------------------
  % W(t)
  % ----------------------------------------------------------------------
  % xp, yp arcsec to radians
  xp_rad = xp * arcsec2rad;
  yp_rad = yp * arcsec2rad;
  
  % s': TIO locator (radians)
  sp_iau00 = iau_SP00(TT1, TT2);
  
  % W(t)^-1
  RPOM = iau_POM00(xp_rad, yp_rad, sp_iau00);

  % W(t)
  RPOM_T = RPOM';
  
  % ----------------------------------------------------------------------
  % Earth Orientation Matrix (EOM)
  % ----------------------------------------------------------------------
  % [GCRS] = Q(t) * R(t) * W(t) * [ITRS] = EOM * {ITRS}
  Qt = RC2I_T;
  Qt_inv = RC2I;
  
  Rt = R_era_n;
  Rt_inv = R_era;
  
  Wt = RPOM_T;
  Wt_inv = RPOM;
  
  
  
  % ----------------------------------------------------------------------
  % GCRS to ITRS transformation matrix
  % ----------------------------------------------------------------------
  %SUBROUTINE iau_RXR ( A, B, ATB )
  GCRS2TIRS = R_era * RC2I;
  CRS2TRS = RPOM * GCRS2TIRS;
  % ----------------------------------------------------------------------
  
  % ----------------------------------------------------------------------
  % ITRS to GCRS transformation matrix (Inverse transformation)
  TRS2CRS = CRS2TRS';
  % ----------------------------------------------------------------------
  
  
  % ----------------------------------------------------------------------
  % Alternative approach
  % ----------------------------------------------------------------------
  
  % ----------------------------------------------------------------------
  % ITRS to GCRS transformation matrix
  % ----------------------------------------------------------------------
  % [GCRS] = Q(t) * R(t) * W(t) * [ITRS] = EOM * {ITRS}
  % ----------------------------------------------------------------------
  % TRS2CRS = Qt * Rt * Wt
  QR = Qt * Rt;
  TRS2CRS = QR * Wt;
  % ----------------------------------------------------------------------
  
  % ----------------------------------------------------------------------
  % GCRS to ITRS transformation matrix
  % ----------------------------------------------------------------------
  % [ITRS] = W(t)^T * R(t)^T * Q(t)^T * [GCRS] = W(t)^T * R3(ERA) * Q(t)^T * [GCRS]
  %*        [TRS]  =  RPOM * R_3(ERA) * RC2I * [CRS]
  % ----------------------------------------------------------------------
  % Inverse transformation matrix as transpose matrix (orthogonality)
   %CALL iau_TR ( TRS2CRS, CRS2TRS )
  
  % or (analytical form)
  
  % CRS2TRS = Wt_inv * Rt_inv * Qt_inv
  % CRS2TRS = RPOM * R_era * RC2I
  Wi_Ri = Wt_inv * Rt_inv;
  CRS2TRS = Wi_Ri * Qt_inv;
  % ----------------------------------------------------------------------
  
  
  
  
  % ----------------------------------------------------------------------
  % Derivatives of CRS-TRS transformation matrix (direct/inverse)
  % ----------------------------------------------------------------------
  
  % ----------------------------------------------------------------------
  % Earth angular velocity
  % ----------------------------------------------------------------------
  % IERS 2010;
  %      dtheta = 0.7292115D0 * 10.0D-04
  % ----------------------------------------------------------------------
  % Omega as derivative of Earth Rotation Angle (ERA), dERA/dt in rad/s
  dtheta = 2.0D0 * pi * 1.00273781191135448D0 * (1.0D0 / 86400D0);
  % ----------------------------------------------------------------------
  
  % ----------------------------------------------------------------------
  % d_TRS2CRS
  % ----------------------------------------------------------------------
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Computation of derivative of TRS2CRS matrix
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % d_TRS2CRS = d(TRS2CRS)/dt = Qt * (dtheta * P * Rt) * Wt
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ----------------------------------------------------------------------
  %P = [ 0  -1   0
  %      1   0   0
  %      0   0   0 ] ;
  % ----------------------------------------------------------------------
  P_dR (1,1) =  0.0;
  P_dR (1,2) = -1.0;
  P_dR (1,3) =  0.0;
  P_dR (2,1) =  1.0;
  P_dR (2,2) =  0.0;
  P_dR (2,3) =  0.0;
  P_dR (3,1) =  0.0;
  P_dR (3,2) =  0.0;
  P_dR (3,3) =  0.0;
  % ----------------------------------------------------------------------
  % d_TRS2CRS = dtheta * Qt * P * Rt * Wt
  QP = Qt * P_dR;  
  QPR = QP * Rt;  
  QPRW = QPR * Wt;  
  d_TRS2CRS = dtheta * QPRW;
  % ----------------------------------------------------------------------
  
  
  % ----------------------------------------------------------------------
  % d_CRS2TRS 
  % ----------------------------------------------------------------------
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Computation of derivative of CRS2TRS matrix
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % d_CRS2TRS = Wt_inv * (omega * P^T * Rt_inv) * Qt_inv
  % d_CRS2TRS = omega * W(-t) * P^T * R(-t) * Q(-t) ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  P_dR_T = P_dR;
  Ri_Qi = Rt_inv * Qt_inv; 
  PT_Ri_Qi = P_dR_T * Ri_Qi; 
  Wi_PT_RiQi = Wt_inv * PT_Ri_Qi;
  d_CRS2TRS = dtheta * Wi_PT_RiQi;
  % ----------------------------------------------------------------------

end