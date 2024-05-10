function [fx,fy,fz] = force_ant(r)
% ----------------------------------------------------------------------
% Purpose:
% Acceleration due to the antenna thrust effect  由于卫星天线辐射推力影响的加速度
% ----------------------------------------------------------------------
% Input arguments:
% - r            : satellite position vector (m)   卫星位置矢量
% 
% Output arguments:
%
% - fx,fy,fz:	 : Accelerations in the inertial frame (m/s^2)  在惯性系中的加速度
% 
% ----------------------------------------------------------------------
  global POWER MASS cslight 

  % Radial vector of satellite
  rsat = sqrt(r(1)^2+r(2)^2+r(3)^2);
  er(1) = r(1)/sqrt(r(1)^2+r(2)^2+r(3)^2);
  er(2) = r(2)/sqrt(r(1)^2+r(2)^2+r(3)^2);
  er(3) = r(3)/sqrt(r(1)^2+r(2)^2+r(3)^2);

  % The acceleration caused by the satellite antenna thrust
  f_ant = POWER/(MASS*cslight);

  % forces in the inertial frame
  fx = (f_ant)*er(1);
  fy = (f_ant)*er(2);
  fz = (f_ant)*er(3);
  
end