function [fx,fy,fz] = force_gm(GM, r)
% ----------------------------------------------------------------------
% Purpose:
% Acceleration due to the central Earth Gravity Field 由地球重力场中心引起的(地球的中心引力)加速度
%  Computation of satellite's acceleration based on Newton's law of gravity
%  considering Earth as a point mass  基于牛顿万有引力定律计算卫星的加速度
% ----------------------------------------------------------------------
% Input arguments:
% - r:				position vector (m)  位置矢量
% 
% Output arguments:
% - fx,fy,fz:		Acceleration's cartesian components (m/s^2)  加速度的笛卡尔坐标分量
% ----------------------------------------------------------------------
   % 计算大地坐标系的坐标
   [phi,lamda,radius] = coord_r2sph(r);

   % Gradient of geopotential V  地球重力势的梯度 V
   fr = - GM / (radius ^ 2);   % 详见《卫星轨道——模型、方法和应用》中公式3.2，p13
   
   %     GM: 地球引力系数 单位：m^3/sec^2     
   %     radius: 卫星到地心的距离，m
   %     phi：大地纬度，rad
   %     lamda：大地经度，rad
   % Cartesian counterparts (fx,fy,fz) of acceleration fr   fr加速度的笛卡尔坐标分量
   fx = fr * cos(phi) * cos(lamda);    % 详见《卫星轨道——模型、方法和应用》中公式2.45，p22
   fy = fr * cos(phi) * sin(lamda);
   fz = fr * sin(phi);

end