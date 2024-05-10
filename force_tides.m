function [fx, fy, fz] = force_tides(r, GM, ae, n_max, m_max, dCnm, dSnm)
% ! ----------------------------------------------------------------------
% ! Purpose:
% ! Tides acceleration
% !   Tides acceleration components are computed as partial derivatives of
% !   spherical harmonics expansion of the tides harmonics coefficients.
% ! ----------------------------------------------------------------------
% ! Input arguments:
% ! - r:			position vector (m) in Terrestrial Reference System (ITRS)
% !   				r = [x y z]
% !
% ! Output arguments:
% ! - fx,fy,fz:		Tides acceleration's Cartesian components (m)
% ! ----------------------------------------------------------------------

  % Tides corrections Nmax
  n1 = size(dCnm, 1);
  Nmax_tide = n1 - 1;

  % computation of spherical coordinates in radians
  [phi, lamda, l] = coord_r2sph (r);

  % computation of normalized associated Legendre functions
  Pnm_norm = legendre (phi, n_max);

  % First-order derivatives of normalized associated Legendre functions
  dPnm_norm = legendre_drv1 (phi, n_max);

  % Partial derivatives of potential w.r.t. the spherical coordinates (2nd approach):
  % - dV_r     : partial derivative w.r.t. radius
  % - dV_phi   : partial derivative w.r.t. latitude
  % - dV_lamda : partial derivative w.r.t. longtitude
  dV_r = 0.0;
  dV_phi = 0.0;
  dV_lamda = 0.0;
  for n = 2:n_max
      if (n > m_max)
          m_limit = m_max;
      else
          m_limit = n;
      end

      for m = 0:m_limit    
          dV_r = dV_r + (-1.0D0 * (GM/l^2)) * (n+1)*1.0D0*((ae/l)^n) * Pnm_norm(n+1,m+1) ...
  		              * (dCnm(n+1,m+1) * cos(m*lamda) + dSnm(n+1,m+1) * sin(m*lamda)); 
          dV_phi = dV_phi + (GM / l) * ((ae/l)^n) * dPnm_norm(n+1,m+1) ...
  		                  * (dCnm(n+1,m+1)*cos(m*lamda) + dSnm(n+1,m+1) * sin(m*lamda)); 
          dV_lamda = dV_lamda + (GM / l) * m * ((ae/l)^n) * Pnm_norm(n+1,m+1) ... 
  		                      * (dSnm(n+1,m+1) * cos(m*lamda) - dCnm(n+1,m+1) * sin(m*lamda)); 
      end
  end

  % Partial derivatives of (r,phi,lamda) with respect to (x,y,z)
  PDVrx(1,1:3) = [cos(phi)*cos(lamda) , cos(phi)*sin(lamda) , sin(phi)]; 
  PDVrx(2,1:3) = [( 1.0D0/l)*sin(phi)*cos(lamda) , ( 1.0D0/l)*sin(phi)*sin(lamda) , (-1.0D0/l)*cos(phi)];
  PDVrx(3,1:3) = [( -1.0D0/(l*cos(phi)) )*sin(lamda) , ( 1.0D0/(l*cos(phi)) )*cos(lamda) , 0.0];
  
  % Computation of Cartesian counterparts of the acceleration
  dV_rpl = [dV_r, dV_phi, dV_lamda];
  PDVrx_transp = PDVrx';
  fx = PDVrx_transp(1,1) * dV_r + PDVrx_transp(1,2) * dV_phi + PDVrx_transp(1,3) * dV_lamda;
  fy = PDVrx_transp(2,1) * dV_r + PDVrx_transp(2,2) * dV_phi + PDVrx_transp(2,3) * dV_lamda;
  fz = PDVrx_transp(3,1) * dV_r + PDVrx_transp(3,2) * dV_phi + PDVrx_transp(3,3) * dV_lamda;


end