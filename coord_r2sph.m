function [phi,lamda,radius] = coord_r2sph(r)
% ! ----------------------------------------------------------------------
% ! Purpose:
% ! Geocentric spherical coordinates 大地坐标系
% !  Computation of the (geocentric) spherical coordinates i.e. longitude
% !  and latitude, from position vector components (Cartesian coordinates)
% ! ----------------------------------------------------------------------
% ! Input arguments:
% ! - r:				Position vector 位置矢量  r = [x y z]
% !
% ! Output arguments:
% ! - phi:			Latitude  (radians) 纬度(rad)
% ! - lamda:			Longitude (radians) 经度(rad)
% ! - radius:			卫星到地心的距离(m)
% ! ---------------------------------------------------------------------------
      x = r(1);
      y = r(2);
      z = r(3);

   % 详见《卫星轨道——模型、方法和应用》中公式2.46，p22
   % radius: 卫星到地心的距离
      radius = sqrt(x^2 + y^2 + z^2);       % r=sqrt(x^2+y^2+z^2)

   % Phi computation: analytical  
      phi = atan( z / sqrt(x^2 + y^2) );    % phi=z/sqrt(x^2+y^2)

   % Lamda computation
      lamda = arctan (y, x);                   % lamda=arctan(y/x)

end