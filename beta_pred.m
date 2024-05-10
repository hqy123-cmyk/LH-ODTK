function beta = beta_pred(mjd, r_sat, v_sat, fti, ft0)
% ----------------------------------------------------------------------
% Purpose:
%  Computation of the beta angle (Sun elevation angle w.r.t. orbtial plane) at a previous of future epoch
%  based on the Keplerian elements for the satellite orbit propagation and the DE Ephemeris for the Sun position
% ----------------------------------------------------------------------
% Input arguments:
% - mjd:		Modified Julian Day number in GPS Time (including the fraction of the day)
% - r_sat: 		Satellite position vector (m) in ITRF
% - v_sat: 		Satellite velocity vector (m/sec) in ITRF
% - fti:		Orbit angle (degrees) at the input epoch ti 
% - ft0:		Orbit angle (degrees) at the required epoch t0 (forward or backward)
%
% Output arguments:
% - beta:		Sun elevation angle with respect to the orbital plane at the epoch t0 (degrees) 
% ----------------------------------------------------------------------
% Note 1:
% - kepler:		Kepler elements array (6)
%   a:       	semi-major axis  (m)
%   e:      	eccentricity
%   i:      	inclination (degrees)
%   Omega:  	right ascension of the ascending node (degrees)
%   omega:  	argument of perigee (degrees)
%   E:      	eccentric anomaly (degrees)
% ----------------------------------------------------------------------
  
  global const;
  % ----------------------------------------------------------------------
  % Current epoch in seconds since start of the day (0h)
  ti = (mjd - round(mjd)) * 86400.D0;
  % ----------------------------------------------------------------------


  % ----------------------------------------------------------------------
  % Computation of the forward/backward epoch t0
  % ----------------------------------------------------------------------
  % Orbital angle rate (in degrees/sec)
  frate = sqrt((v_sat(1)^2+v_sat(2)^2+v_sat(3)^2) / (r_sat(1)^2+r_sat(2)^2+r_sat(3)^2)) * (180.0D0 / pi);
  
  t0 = ti + (ft0 - fti)/ frate;
  % ----------------------------------------------------------------------

  % ----------------------------------------------------------------------
  % Computation of the satellite position and velocity vectors at epoch t0
  % through Keplerian orbit consideration
  % ----------------------------------------------------------------------
  % Gravitational constant
  GM = const.GM_Earth;
  
  % Kepler elements at ti
  kepler_i = kepler_z2k(r_sat, v_sat, GM);
  
  a_semiaxis = kepler_i(1);
  ec = kepler_i(2);
  
  % Eccentric anomaly at ti (radians)
  Ei_rad = kepler_i(7) * (pi / 180.D0);
  
  % Mean anomaly  at ti (radians)
  Mi_rad = Ei_rad - ec * sin(Ei_rad);
  
  % mean motion (rad/sec)
  n_motion = sqrt(GM / a_semiaxis^3);
  
  % ----------------------------------------------------------------------
  % Time difference since current epoch
  dt = t0 - ti;
  % ----------------------------------------------------------------------

  % ----------------------------------------------------------------------
  % Mean anomaly (radians) at epoch t0 
  Mo_rad = Mi_rad + n_motion * dt;
  
  % Mean anomaly (in degrees) at epoch t0
  Mo_deg = Mo_rad * (180.D0/pi);
  
  % Angle reduction within the range {0 - 360} degrees
  if (abs(Mo_deg) >= 360.D0)
      Mo_deg = Mo_deg - round(Mo_deg / 360.D0) * 360.D0;
  end

  if (Mo_deg < 0.0D0)
      Mo_deg = Mo_deg + 360.D0;
  end
  % ----------------------------------------------------------------------


  % Kepler Equation 
  % Eccentric anomaly (degrees) at epoch t0 obtained by solving Kepler's Equation
  Eo_deg = kepler_eq(Mo_deg, ec);
  
  % Kepler elements at epoch t0
  kepler_0 = [kepler_i(1), kepler_i(2), kepler_i(3), kepler_i(4), kepler_i(5), Eo_deg];
  
  % State vector (Cartesian coordinates) at epoch t0
  [rsat_t0, vsat_t0] = kepler_k2z(kepler_0, GM);
  % ----------------------------------------------------------------------

  % ----------------------------------------------------------------------
  % Sun position vector at reference epoch t0
  % ----------------------------------------------------------------------
  
  % Time System transformation: GPS to TT
  mjd_t0 = mjd + dt/86400.0D0;
  [mjd_TT, ~, ~, ~] = time_GPS (mjd_t0);
  
  [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, rsun_t0, ~] = JPL_Eph_DE440(mjd_TT);
  % ----------------------------------------------------------------------
  % Sun elevation angle at epoch t0
  % ----------------------------------------------------------------------
  % Beta angle (in degrees)
  beta = beta_angle(rsat_t0, vsat_t0, rsun_t0);
  % ----------------------------------------------------------------------
end