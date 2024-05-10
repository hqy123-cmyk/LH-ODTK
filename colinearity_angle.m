function [epsilon_angle, epsilon_dot] = colinearity_angle(r_sat, v_sat, r_sun)
% ----------------------------------------------------------------------
% Purpose:
%  Computation of the colinearity angle between the orbital plane and the Sun vector 
% ----------------------------------------------------------------------
% Input arguments:
% - r_sat: 			Satellite position vector (m) 
% - v_sat: 			Satellite velocity vector (m/sec)
% - r_sun:			Sun position vector (m)
%
% Output arguments:
% - epsilon_angle:		Colinearity angle in degrees
% - epsilon_angle:		Colinearity angle partial derivative w.r.t. time 
% ----------------------------------------------------------------------

  % ----------------------------------------------------------------------
  % Unit Vectors
  % ----------------------------------------------------------------------
  % Xsat:		Satellite position unit vector
  Xsat = (1D0/sqrt(r_sat(1)^2 + r_sat(2)^2 + r_sat(3)^2) ) * r_sat;
  	  
  % dXsat:	Satellite Velocity unit vector
  dXsat = (1D0/sqrt(v_sat(1)^2 + v_sat(2)^2 + v_sat(3)^2) ) * v_sat;
  	  
  % Xsun:		Sun position unit vector
  Xsun = (1D0/sqrt(r_sun(1)^2 + r_sun(2)^2 + r_sun(3)^2) ) * r_sun;
  % ----------------------------------------------------------------------

  % ----------------------------------------------------------------------
  % Orbit normal vector (cross-track)
  % ----------------------------------------------------------------------
  % cross product of r,v : angular momentum per unit mass
  h_angmom = cross(r_sat, v_sat);
  h_dot = dot(h_angmom, h_angmom);
  h_norm = sqrt(h_dot);
  
  % Cross-track or normal component
  en = (1D0 / h_norm) * h_angmom;
  % ----------------------------------------------------------------------

  % ----------------------------------------------------------------------
  % Colinearity angle: epsilon
  % ----------------------------------------------------------------------
  x_cross = cross(en, Xsun);
  y_cross = cross(en, x_cross);
  ry_dot = dot(Xsat, y_cross);
  vy_dot = dot(dXsat, y_cross);

  acos_ry = acos(ry_dot) * (180.0D0/pi);
  
  % Epsilon angle and partial derivative of epsilon w.r.t. time
  if (acos_ry <= 90.0D0)
      epsilon_angle = acos_ry;
  	  epsilon_dot = (-1.0D0/sin(epsilon_angle)) * vy_dot;
  else
      epsilon_angle = 180.0D0 - acos_ry;
  	  epsilon_dot = (1.0D0/sin(epsilon_angle)) * vy_dot;
  end
end