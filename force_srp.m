function [fx, fy, fz] = force_srp(lambda, eBX_ecl, eclipsf, GM, GNSSid, srpid, r, v, r_sun)
% ! ----------------------------------------------------------------------
% ! Purpose:
% ! Acceleration due to the solar radiation pressure 
% ! Computation of SRP acceleration using various SRP models, 
% ! such as a simply cannonball model (srpid=1),a box-wing model (srpid=2) 
% ! and a ECOM model (srpid=3) 
% ! ----------------------------------------------------------------------
% ! Input arguments:
% ! - GM           : the earth gravitational constant  地球引力常数
% ! - srpid        : =1: a simply cannonball model;  简易模型
% !                  =2: box-wing model;
% !                  =3: ECOM model
% ! - GNSSid       : id of satellite constellation
% ! - r            : satellite position vector (m)  卫星位置矢量
% ! - v            : satellite velocity vector   卫星速度矢量
% ! - r_sun        : Sun position vector  太阳位置矢量
% ! - lambda       : shadow coefficient  阴影系数
% ! - eBX_ecl      : dynamic ex of satellite body frame
% ! 
% ! Output arguments:
% ! - fx,fy,fz:		Acceleration's cartesian components (m)
% ! ----------------------------------------------------------------------  

  % Numerical Constants
    Cr = 1.5; % SRP coefficient ranges from 1.3 to 1.5
  % ---------------------------------------------------------------------
    ex_i = 0;   % change the definition of the unit vector ex 
                % ex_i = 0 (default)
                %      = 1 (using dynamic ex vector from attitude routine)
    att_ON = 1; % att_ON = 1 : use the orbit-normal attitude for BDS satellite
                %              when the beta < 4 deg
                %        = 0 : use the yaw-steering attitude for BDS satellite 
                %              for all beta angles 
  % ---------------------------------------------------------------------
  
  global const;
  global BLKTYP
  global srpcoef ECOMNUM yml_SRP_mode SRP;
       
  % initialize the SRP force
  fsrp = zeros(3);
      
  % Initialise
  fx = 0.d0;
  fy = 0.d0;
  fz = 0.d0;
  [X_SIDE, Z_SIDE, A_SOLAR, F0] = apr_srp(GNSSid, BLKTYP);
  
  % The unit vector ez SAT->EARTH
  er(1)=r(1)/sqrt(r(1)^2+r(2)^2+r(3)^2);
  er(2)=r(2)/sqrt(r(1)^2+r(2)^2+r(3)^2);
  er(3)=r(3)/sqrt(r(1)^2+r(2)^2+r(3)^2);
  ez(1)=-er(1);
  ez(2)=-er(2);
  ez(3)=-er(3);
  
  % The unit vector ed SAT->SUN (where is just opposite to the solar radiation vector)
  Ds=sqrt((r_sun(1)-r(1))^2+(r_sun(2)-r(2))^2+(r_sun(3)-r(3))^2);
  ed(1)=((r_sun(1)-r(1))/Ds);
  ed(2)=((r_sun(2)-r(2))/Ds);
  ed(3)=((r_sun(3)-r(3))/Ds);

  % The unit vector ey = ez x ed/|ez x ed|, parallel to the rotation axis of solar panel
  % If the Y-axis is reversed for IIR satellites, the result is identical to the
  % IGS-defined. We have done the test (14/10/2020, Tzupang Tseng)
  yy = cross(ez,ed);
  ey(1)=yy(1)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
  ey(2)=yy(2)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
  ey(3)=yy(3)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);

  % The unit vector of x side, which is always illuminated by the sun.
  ex = cross(ey,ez);

  % Using the BODY-X univector from Kouba routine to redefine the ey for the satellite eclipsed
  if (ex_i > 0)
      ex(1)=eBX_ecl(1)/sqrt(eBX_ecl(1)^2+eBX_ecl(2)^2+eBX_ecl(3)^2);
      ex(2)=eBX_ecl(2)/sqrt(eBX_ecl(1)^2+eBX_ecl(2)^2+eBX_ecl(3)^2);
      ex(3)=eBX_ecl(3)/sqrt(eBX_ecl(1)^2+eBX_ecl(2)^2+eBX_ecl(3)^2);
      ey = productcross (ez, ex);
  end
  % The unit vector eb = ed x ey
  eb = cross(ed,ey);

  % the orbit normal vector  轨道的法向量
  ev(1)=v(1)/sqrt(v(1)^2+v(2)^2+v(3)^2);
  ev(2)=v(2)/sqrt(v(1)^2+v(2)^2+v(3)^2);
  ev(3)=v(3)/sqrt(v(1)^2+v(2)^2+v(3)^2);
  en = cross(er, ev);

  % computation of the satellite argument of latitude and orbit inclination
  kepler = kepler_z2k(r,v,GM);
  u_sat = kepler(9)*pi/180.d0;
  i_sat = kepler(3)*pi/180.d0;
  omega_sat = kepler(4)*pi/180.d0;
  
  % compute the sun position in the satellite orbit plane by rotating omega_sat and i_sat,
  % allowing us for the computation of u_sun and sun elevation angles (beta)    
  for i=1:3
      for j=1:3
          R33(i,j)=1.0d0;
      end
  end
  R33(1,3)=0.0d0;
  R33(2,3)=0.0d0;
  R33(3,1)=0.0d0;
  R33(3,2)=0.0d0;
  R33(3,3)=1.0d0;

  for i=1:3
     for j=1:3
         R11(i,j)=1.0d0;
     end
  end
  R11(1,1)=1.0d0;
  R11(1,2)=0.0d0;
  R11(1,3)=0.0d0;
  R11(2,1)=0.0d0;
  R11(3,1)=0.0d0;

  % rotate big omega to make the X-axis of the celestial frame consistent with the
  % direction of right ascension of ascending node
  R33 = R3(omega_sat);
  
  % rotate inclination to make the XY plane of the celestial frame consistent with
  % the orbit plane
  R11 = R1(i_sat);
  
  % convert the sun position in the celestial frame to the orbit plane
  r_sun1 = R33 * r_sun;
  r_sun2 = R11 * r_sun1;

  beta = atan2(r_sun2(3),sqrt(r_sun2(1)^2+r_sun2(2)^2));   % in rad 
  % write (*,*) beta*180.0d0/Pi
  u_sun = atan2(r_sun2(2),r_sun2(1));                       % in rad

  %  del_u = u_sat - u_sun
  if (yml_SRP_mode == const.ECOM1 || yml_SRP_mode == const.BERNE || ...
                        yml_SRP_mode ~= const.SBOXW)
      del_u = u_sat;
  else
      del_u = u_sat - u_sun;
  end

  if (del_u*180/pi > 360.0d0)
      del_u = del_u-2*Pi;
  elseif(del_u*180/pi < 0.0d0)
      del_u=del_u+2*pi;
  end

  % Implement the orbit-normal attitude for BDS satellites when the beat < 4 deg
  if (att_ON == 1)  
     if(BLKID == 301 || BLKID == 307)
         yy = cross (ez,ev);
         ey(1)=yy(1)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
         ey(2)=yy(2)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
         ey(3)=yy(3)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);

         yy = cross (ed,ey);
         eb(1)=yy(1)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
         eb(2)=yy(2)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
         eb(3)=yy(3)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);

         yy = cross (ey,eb);
         ed(1)=yy(1)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
         ed(2)=yy(2)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
         ed(3)=yy(3)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);

     elseif(BLKID == 302 || BLKID == 303 || BLKID == 304 || ...
                            BLKID == 305 || BLKID == 306)
         if (eclipsf == 3)
             yy = cross (ez,ev);
             ey(1)=yy(1)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
             ey(2)=yy(2)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
             ey(3)=yy(3)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);

             yy = cross (ed,ey);
             eb(1)=yy(1)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
             eb(2)=yy(2)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
             eb(3)=yy(3)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);

             % Switch on the (ed,on) does not improve the orbit accuracy!!
             % Only for testing purpose!      
             yy = cross (ey,eb);
             ed(1)=yy(1)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
             ed(2)=yy(2)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
             ed(3)=yy(3)/sqrt(yy(1)^2+yy(2)^2+yy(3)^2);
         end
     end
  end
  
  % ========================================
  %  angles between ed and each surface(k)
  %  k=1: +X
  %    2: +Yi
  %    3: +Z
  %    4: solar panels
  % =======================================

  cosang(1)=ed(1)*ex(1)+ed(2)*ex(2)+ed(3)*ex(3);
  cosang(2)=ed(1)*ey(1)+ed(2)*ey(2)+ed(3)*ey(3);
  cosang(3)=ed(1)*ez(1)+ed(2)*ez(2)+ed(3)*ez(3);
  cosang(4)=ed(1)*ed(1)+ed(2)*ed(2)+ed(3)*ed(3);


  % A scaling factor is applied to ECOM model
  sclfa=(AU/Ds)^2;
  alpha = 1.0d0;   % srpid == SRP_NONE

  if (srpid == const.SRP_CANNONBALL)
      fx=-F0/MASS*ed(1)*lambda;
      fy=-F0/MASS*ed(2)*lambda;
      fz=-F0/MASS*ed(3)*lambda;
      alpha = F0/MASS;
  elseif (srpid == const.SRP_SIMPLE_BW)
      xmul = 0.02;
      zmul = 0.02;
      solarmul = 1.7;
      if(cosang(3)<0.d0)
          ez = ez*(-1.d0);
      end

      fxo=Ps/MASS*(xmul*X_SIDE*cosang(1)*ex(1)+zmul*Z_SIDE*cosang(3)*ez(1)+solarmul*A_SOLAR*cosang(4)*ed(1));
      fyo=Ps/MASS*(xmul*X_SIDE*cosang(1)*ex(2)+zmul*Z_SIDE*cosang(3)*ez(2)+solarmul*A_SOLAR*cosang(4)*ed(2));
      fzo=Ps/MASS*(xmul*X_SIDE*cosang(1)*ex(3)+zmul*Z_SIDE*cosang(3)*ez(3)+solarmul*A_SOLAR*cosang(4)*ed(3));
      alpha = sqrt(fxo^2+fyo^2+fzo^2);
      fx=-fxo*lambda;
      fy=-fyo*lambda;
      fz=-fzo*lambda;
  elseif (srpid == const.SRP_FULL_BW)
      REFF = 0;    % 0: inertial frame,  1: satellite body-fixed frame, 
                   % 2: sun-fixed frame, 3: orbital frame 
      YSAT(1:3) = r;
      YSAT(4:6) = v;
      ACCEL = SRPFBOXW(REFF,YSAT,R_SUN,SVNID);
      alpha = sqrt(ACCEL(1)^2+ACCEL(2)^2+ACCEL(3)^2);
      fx=ACCEL(1)*lambda;
      fy=ACCEL(2)*lambda;
      fz=ACCEL(3)*lambda;
  end

  D0 = alpha;
  srpcoef= zeros(ECOMNUM);
  % ECOM model
  if (yml_SRP_mode ~= const.ECOM_NONE && yml_SRP_mode ~= const.SBOXW)
        % ie ECOM1, ECOM2 or ECOM_HYBRID
        PD_Param_ID = 0;
        if (SRP.bias_D)
            PD_Param_ID = PD_Param_ID + 1;
            srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
            if (lambda < 1)
                srpcoef(PD_Param_ID) = lambda*srpcoef(PD_Param_ID);
            end
            for i=1:3
               fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ed(i)*alpha;
            end
            %print*,'ECOM1-caused accelerations'
        end

        if (SRP.bias_Y)
            PD_Param_ID = PD_Param_ID + 1;
            srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
            for i=1:3
                fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ey(i)*alpha;
            end
        end

        if (SRP.bias_B)
            PD_Param_ID = PD_Param_ID + 1;
            srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
            for i=1:3
                fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*eb(i)*alpha;
            end
        end

        % DC term
        if(SRP.cpr_DC)
            PD_Param_ID = PD_Param_ID + 1;
            srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
            for i=1:3
                fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*cos(del_u)*ed(i)*alpha;
            end
        end

        % DS term
        if(SRP.cpr_DS)
           PD_Param_ID = PD_Param_ID + 1;
           srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
           for i=1:3
               fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*sin(del_u)*ed(i)*alpha;
           end
        end
        
        % YC term
        if(SRP.cpr_YC)
           PD_Param_ID = PD_Param_ID + 1;
           srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
           for i=1:3
               fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*cos(del_u)*ey(i)*alpha;
           end
        end

        % YS term
        if(SRP.cpr_YS)
           PD_Param_ID = PD_Param_ID + 1;
           srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
           for i=1:3
               fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*sin(del_u)*ey(i)*alpha;
           end
        end
        
        % BC term
        if(SRP.cpr_BC)
           PD_Param_ID = PD_Param_ID + 1;
           srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
           for i=1:3
               fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*cos(del_u)*eb(i)*alpha;
           end
        end

        % BS term
        if(SRP.cpr_BS)
           PD_Param_ID = PD_Param_ID + 1;
           srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
           for i=1:3
               fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*sin(del_u)*eb(i)*alpha;
           end
        end

        % D2C term
        if(SRP.cpr_D2C)
           PD_Param_ID = PD_Param_ID + 1;
           srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
           for i=1:3
               fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*cos(2*del_u)*ed(i)*alpha;
           end
        end

        % D2S term
        if(SRP.cpr_D2S)
           PD_Param_ID = PD_Param_ID + 1;
           srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
           for i=1:3
               fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*sin(2*del_u)*ed(i)*alpha;
           end
        end

        % D4C term
        if(SRP.cpr_D4C)
            PD_Param_ID = PD_Param_ID + 1;
            srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
            for i=1:3
                fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*cos(4*del_u)*ed(i)*alpha;
            end
        end

        % D4S term
        if(SRP.cpr_D4S)
            PD_Param_ID = PD_Param_ID + 1;
            srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
            for i=1:3
                fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*sin(4*del_u)*ed(i)*alpha;
            end
        end

        if (ECOMNUM ~= PD_Param_ID)

            fprintf('%s\n',"THE NUMBER OF FORCE PARAMETERS IS NOT CONSISTENT");
            fprintf('%s %d\n',"ECOMNUM     =",ECOMNUM);
            fprintf('%s %d\n',"PD_Param_ID =",PD_Param_ID);
            fprintf('%s\n',"PROGRAM STOP AT force_srp.m");
            return;
        end
   % SIMPLE BOX WING
   elseif (yml_SRP_mode == const.SBOXW)
         
       for PD_Param_ID = 1:7
           if (PD_Param_ID == 1)
               srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
               if (lambda < 1) 
                   srpcoef(PD_Param_ID) = lambda*srpcoef(PD_Param_ID);
               end
               for i=1:3
                   fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ex(i);
               end
           elseif (PD_Param_ID == 2)
               srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
               if (lambda < 1) 
                   srpcoef(PD_Param_ID) = lambda*srpcoef(PD_Param_ID);
               end
               for i=1:3
                   fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ez(i);
               end
           elseif (PD_Param_ID == 3)
               srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
               if (lambda < 1) 
                   srpcoef(PD_Param_ID) = lambda*srpcoef(PD_Param_ID);
               end
               for i=1:3
                   fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ed(i);
               end
           elseif (PD_Param_ID == 4)
               srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
               for i=1:3
                   fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ey(i);
               end
           elseif (PD_Param_ID == 5)
               srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
               for i=1:3
                   fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*eb(i);
               end
           elseif (PD_Param_ID == 6)
               srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
               for i=1:3
                   fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*cos(del_u)*eb(i);
               end
           elseif (PD_Param_ID == 7)
               srpcoef(PD_Param_ID) = ECOM_accel_glb(PD_Param_ID);
               for i=1:3
                   fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*sin(del_u)*eb(i);
               end
           end
       end
   end
   
   fx=-fsrp(1);
   fy=-fsrp(2);
   fz=-fsrp(3);

end