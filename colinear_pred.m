function delta_t = colinear_pred(mjd, r_sat, v_sat, r_sun, epsilon_t0)
% ----------------------------------------------------------------------
% Purpose:
%  Computation of the time of the colinearity angle (Sun elevation angle w.r.t. orbtial plane) at a previous epoch
%  based on the Keplerian elements for the satellite orbit propagation and the DE Ephemeris for the Sun position
% ----------------------------------------------------------------------
% Input arguments:
% - mjd:		Modified Julian Day number in GPS Time (including the fraction of the day)
% - r_sat: 		Satellite position vector (m) in ITRF
% - v_sat: 		Satellite velocity vector (m/sec) in ITRF
% - r_sun:			Sun position vector (m)
% - epsilon_t0:	colinearity angle (epsilon) at the current epoch
%
% Output arguments:
% - delta_t:		Sun elevation angle with respect to the orbital plane at the epoch t0 (degrees) 
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

  % ----------------------------------------------------------------------
  % Current epoch in seconds since start of the day (0h)
  ti = (mjd - round(mjd)) * 86400.D0;
  % ----------------------------------------------------------------------
  
  % ----------------------------------------------------------------------
  % Colinearity angle at input epoch ti
  [epsilon_ti, epsilon_dot] = colinearity_angle(r_sat, v_sat, r_sun);
  % ----------------------------------------------------------------------
  
  % ----------------------------------------------------------------------
  % Time span
  delta_t = (epsilon_ti - epsilon_t0)  / epsilon_dot;
  % ----------------------------------------------------------------------

end