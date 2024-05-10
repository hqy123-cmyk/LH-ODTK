function [fx,fy,fz] = force_gfm(GM, ae, r, n_max, m_max, Cnm, Snm)
% ! Purpose:
% ! Earth Gravity field effects
% !  Computation of the acceleration based on a global Earth Gravity Field  
% !  Model expressed by a set of spherical harmonic coefficients
% ! ----------------------------------------------------------------------
% ! Input arguments:
% ! - GM:				Earth gravity constant  (m^3/sec^2)
% ! - ae:				Earth radius (meters)
% ! - r:				Position vector (m) in Terrestrial Reference System (ITRS)
% !   				r = [x y z]
% ! - n_max:          maximum degree expansion
% ! - m_max:          maximum order expansion (m_max <= n_max)
% ! - Cnm, Snm:		Spherical Harmonic Coefficients (degree n, order m); dynamic allocatable arrays
% ! 
% ! Output arguments:
% ! - fx,fy,fz:		Acceleration's cartesian components in ITRS (m)

    % computation of spherical coordinates in radians  计算大地坐标系坐标
    %     radius: 卫星到地心的距离, m
    %     phi：大地纬度, rad
    %     lamda：大地经度, rad 
    [phi,lamda,radius] = coord_r2sph(r);
    
    % computation of normalized associated Legendre functions  计算归一化后的缔合勒让德函数
    Pnm_norm = legendre (phi, n_max);

    % First-order derivatives of normalized associated Legendre functions  归一化后的缔合勒让德函数的一阶导数
    dPnm_norm = legendre_drv1 (phi, n_max);

    % Computation of acceleration vector components in a local tangent frame:
    % 在一个局部切线坐标框架中，计算加速度矢量的分量:
    comp_option = 2;

    % 1st approach:
    if (comp_option == 1)
        % ----------------------------------------------------------------------
        % - fr     : radius component  卫星与地心连线方向的分量
        % - flamda : geocentric longitude component  大地经度方向的分量
        % - ftheta : geocentric latitude component (theta(余纬度) = 90-latitude)  大地纬度方向的分量
        % ----------------------------------------------------------------------
        fr = 0.0D0;
        ftheta = 0.0D0;
        flamda = 0.0D0;
        % ----------------------------------------------------------------------
        for n = 2:n_max
            if (n > m_max)
                m_limit = m_max;
            else
                m_limit = n;
            end
            for m = 0:m_limit
             fr = fr + ( -1.0D0*(n+1) * ((ae/radius)^n) * Pnm_norm(n+1,m+1) ...         % ae：地球参考椭球体的赤道半径
			         * (Cnm(n+1,m+1)*cos(m*lamda) + Snm(n+1,m+1)*sin(m*lamda)) );
	 
             ftheta = ftheta + ((ae/radius)^n) * dPnm_norm(n+1,m+1) ...
			                 * (Cnm(n+1,m+1)*cos(m*lamda)+Snm(n+1,m+1)*sin(m*lamda));

             flamda = flamda + ((ae/radius)^n) * (1.0D0/cos(phi)) ...
     		                 * Pnm_norm(n+1,m+1) * m * (Snm(n+1,m+1)*cos(m*lamda)-Cnm(n+1,m+1)*sin(m*lamda));
            end
        end

        fr = fr * (GM/radius^2) - GM / (radius^2);
        ftheta = ftheta * (GM/radius^2); 
        flamda = flamda * (GM/radius^2);

        % Cartesian counterparts in the Earth fixed system (ITRF) 在地心地固系(ITRF)中，相应的笛卡尔坐标分量
        fx = fr * cos(phi)*cos(lamda) + ftheta * sin(phi)*cos(lamda) - flamda * sin(lamda);
        fy = fr * cos(phi)*sin(lamda) + ftheta * sin(phi)*sin(lamda) + flamda * cos(lamda);
	    fz = fr * sin(phi) - ftheta * cos(phi);
        
    % 2nd approach:
    else  
        % Partial derivatives of potential with respect to spherical coordinates :
        % - dV_r     : partial derivative of geopotential to radius
        % - dV_phi   : partial derivative of geopotential to latitude
        % - dV_lamda : partial derivative of geopotential to longtitude
        dV_r = 0.0D0;
        dV_phi = 0.0D0;
        dV_lamda = 0.0D0;

        for n = 2:n_max
            if (n > m_max)
                m_limit = m_max;
            else
                m_limit = n;
            end
            for m = 0:m_limit
                dV_r = dV_r + ( -1.0D0*(n+1) * ((ae/radius)^n)*Pnm_norm(n+1,m+1) * ...
                              ( Cnm(n+1,m+1) * cos(m*lamda)+ Snm(n+1,m+1) * sin(m*lamda)) );

                dV_phi = dV_phi + ((ae/radius)^n) * dPnm_norm(n+1,m+1) * ...
                                  ( Cnm(n+1,m+1)*cos(m*lamda) + Snm(n+1,m+1)*sin(m*lamda) );
	 
                dV_lamda = dV_lamda + m * ((ae/radius)^n) * Pnm_norm(n+1,m+1) * ...
                           ( Snm(n+1,m+1)*cos(m*lamda) - Cnm(n+1,m+1)*sin(m*lamda) );
            end
        end
        dV_r = - GM / (radius^2) + (GM/radius^2) * dV_r;
        dV_phi = (GM / radius) * dV_phi;
        dV_lamda = (GM / radius) * dV_lamda;

        % Computation of Cartesian counterparts of the acceleration  计算加速度相应的笛卡尔坐标分量
        fx = dV_r * cos(phi)*cos(lamda) + dV_phi * (1.0D0/radius)*sin(phi)*cos(lamda) + ... 
             dV_lamda * (-1.0D0/(radius*cos(phi)))*sin(lamda);
        fy = dV_r * cos(phi)*sin(lamda) + dV_phi * (1.0D0/radius)*sin(phi)*sin(lamda) + ... 
             dV_lamda * (1.0D0/(radius*cos(phi)))*cos(lamda);
        fz = dV_r * sin(phi) + dV_phi * (-1.0D0/radius)*cos(phi);
    end

end