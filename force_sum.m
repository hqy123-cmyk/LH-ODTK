function [SFx, SFy, SFz] = force_sum(mjd, t_sec, rsat, vsat, integr_stage)
%  Purpose:
% !        基于被考虑的动力学模型的卫星加速度分量
% ! ----------------------------------------------------------------------
% ! Input arguments:
% ! - mjd:			简化儒略日数(包括天的小数部分) 地球时(TT)
% !                       
% ! - t_sec:		从这一天00h开始的秒
% ! - rsat:			在惯性系(ICRF)下，卫星的位置矢量(m)
% ! - vsat:			在惯性系(ICRF)下，卫星的速度矢量(m/s)
% ! 
% ! Output arguments:
% ! - fx,fy,fz:		加速度的笛卡尔坐标分量 (m/s^2)
% ! ----------------------------------------------------------------------
    % Global variables
    global const;
    global GFM_GM GFM_ae ppn_be ppn_ga cslight;
    global GFM_Nmax GFM_Mmax GFM_Cnm GFM_Snm;
    global OCEAN_Nmax OCEAN_Mmax;
    global AuxParam;
    global grav_effects non_grav_effects;
    global yml_apriori_srp yml_EMP_mode yml_pulses;

    GMearth = GFM_GM;
    aEarth = GFM_ae;

    % Relativistic parameters
    beta_ppn = ppn_be;
    gama_ppn = ppn_ga;
    c_light = cslight;

    if (FMOD_GRAVFIELD == 0)
        GMearth = GM_global;
        aEarth = Earth_radius;
    end

    % 在ICRF坐标系中的状态矢量
    rsat_icrf = rsat;
    vsat_icrf = vsat;
    rsat_itrf = 0.d0;

   % Init var
   vSun = 0.d0;

   % EOP数据和ITRF与ICRF的旋转矩阵
   if(FMOD_GRAVFIELD > 0 || yml_planetary_perturbations_enabled ...
                         || yml_tidal_effects_enabled || yml_rel_effects_enabled)

       [CRS2TRS, TRS2CRS, d_CRS2TRS, d_TRS2CRS] = EOP(mjd, EOP_cr);
       % 在计算历元时刻EOP数据的改正
       % 极移坐标系(由xy平面坐标定义，其x轴指南且与格林尼治平子午线一致，y轴指西)的坐标
       % xp和yp是真天极的角度值
       xp = EOP_cr(2); 
       yp = EOP_cr(3);
       % UT1-UTC
       ut1_utc = EOP_cr(4);

       % State Vector in ITRF  在ITRF坐标系中的状态矢量
       % r in ITRF
       rsat_itrf = CRS2TRS * rsat_icrf;

       % v in ITRF 
       v_TRS_1 = CRS2TRS*vsat_icrf;
       v_TRS_2 = d_CRS2TRS * rsat_icrf;
       vsat_itrf = v_TRS_1 + v_TRS_2;
   end

   % Gravitational Effects(地球引力效应)

   % Gravity Field  重力场
   if(grav_effects.yml_gravity == 0)
       % Setting to Zero
       fx = 0.0;
       fy = 0.0;
       fz = 0.0;	  
       Fgrav_icrf = [fx, fy, fz];
   else
       % 地球重力场的中心项
       if(FMOD_GRAVFIELD == 0)
           [fx,fy,fz] = force_gm(GM, r);   % 地球中心引力加速度
           Fgrav_icrf = [fx, fy, fz];
       else
           % 地球重力场模型摄动加速度分量(球谐函数表达)
           n_max = GFM_Nmax;
           m_max = GFM_Mmax;
           % 地球中心引力+地球非球形引力摄动 => 地球重力场模型摄动加速度
           [fx, fy, fz] = force_gfm(GMearth, aEarth, rsat_itrf, n_max, m_max, GFM_Cnm, GFM_Snm);
           Fgrav_itrf = [fx, fy, fz];
      
           % 将加速度矢量从地球坐标系转换至惯性坐标系
           Fgrav_icrf = TRS2CRS * Fgrav_itrf;
       end
   end

   % 行星轨道摄动
   if(grav_effects.yml_planetary_perturbations)

       % 计算mjd历元的行星位置, 中心天体: 地球
      [r_Mercury, r_Venus, ~, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, ...
                                                       r_Pluto, r_Moon, r_Sun, ~] = JPL_Eph_DE440(mjd);
      
      % Planets perturbation loop
      aPlanets_icrf = [0.0, 0.00, 0.00];

      % Luni-solar perturbations
      aPlanets_icrf = aPlanets_icrf + force_gm3rd(r_Sun, rsat_icrf, const.GM_Sun);    % Sun
      aPlanets_icrf = aPlanets_icrf + force_gm3rd(r_Moon, rsat_icrf, const.GM_Moon);  % Moon

      % Venus
      if (AuxParam.Venus)
          aPlanets_icrf = aPlanets_icrf + force_gm3rd(r_Venus, rsat_icrf, const.GM_Venus);
      end

      % Jupiter
      if (AuxParam.Jupiter)
          aPlanets_icrf = aPlanets_icrf + force_gm3rd(r_Jupiter, rsat_icrf, const.GM_Jupiter);
      end

      % Mars
      if (AuxParam.Mars)
          aPlanets_icrf = aPlanets_icrf + force_gm3rd(r_Mars, rsat_icrf, const.GM_Mars);
      end
      
      % Other Planetary perturbations
      if (AuxParam.planets)
          aPlanets_icrf = aPlanets_icrf + force_gm3rd(r_Mercury, rsat_icrf, const.GM_Mercury);
          aPlanets_icrf = aPlanets_icrf + force_gm3rd(r_Saturn, rsat_icrf, const.GM_Saturn);
          aPlanets_icrf = aPlanets_icrf + force_gm3rd(r_Uranus, rsat_icrf, const.GM_Uranus);
          aPlanets_icrf = aPlanets_icrf + force_gm3rd(r_Neptune, rsat_icrf, const.GM_Neptune);
          aPlanets_icrf = aPlanets_icrf + force_gm3rd(r_Pluto, rsat_icrf, const.GM_Pluto);
      end

       % Indirect J2 effect J2项摄动的间接影响
       % Sun and Moon position vectors in terrestrial reference frame  在惯性参考系中太阳和月球的位置矢量
       % Transformation GCRS to ITRS  将GCRS转换至ITRS
       rMoon_ITRS = CRS2TRS * r_Moon';
       rSun_ITRS = CRS2TRS * r_Sun';
       
       if (FMOD_GRAVFIELD > 0)
           % C20 spherical harmonic coefficient of geopotential C20 地球重力势的球谐系数
           C20 = GFM_Cnm(2+1,0+1);
       elseif (FMOD_GRAVFIELD == 0)
           C20 = -4.841694552725e-04;
       end

       % Earth Radius 地球半径
       Re = aEarth;
  
       % Indirect J2 effect of Sun and Moon  太阳和月球的间接J2项影响
       a_iJ2 = indirectJ2(C20, Re, GM_moon, rMoon_ITRS, GM_sun, rSun_ITRS);

       % Acceleration vector transformation from terrestrial to inertial frame  将加速度矢量从地球坐标系转换至惯性坐标系
       a_iJ2_icrf = TRS2CRS * a_iJ2;

       % Planetary/Lunar perturbation accleration   行星摄动加速度(第三体引力摄动)
       Fplanets_icrf = aPlanets_icrf + a_iJ2_icrf;

   elseif (~ grav_effects.yml_planetary_perturbations)
       Fplanets_icrf = [0.0, 0.0, 0.0];
   end
   
   % Tidal effects    潮汐影响
   if (grav_effects.yml_tidal_effects)

       % Solid Earth Tides   地球固体潮
       % Frequency-independent 计算与频率无关的球谐系数改正 : Step 1 (IERS Conventions 2010)
       dCnm_solid1 = 0.0;
       dSnm_solid1 = 0.0;
       if (grav_effects.yml_solid_nonfreq)
           [dCnm_solid1, dSnm_solid1] = tides_solid1(rMoon_ITRS, rSun_ITRS, GMearth, aEarth, GM_moon, GM_sun);
       end

       % Frequency-dependent 计算与频率有关的二阶项改正 : Step 2 (IERS Conventions 2010)
       dCnm_solid2 = 0.0;
       dSnm_solid2 = 0.0;
       if (grav_effects.yml_solid_freq)
           [dCnm_solid2, dSnm_solid2] = tides_solid2(mjd, ut1_utc);
       end
 
       % Zero Tide term  0潮汐项——计算由永久潮汐引起的二阶项改正
       dC20_perm = tide_perm();
 
       % Tide system of the input gravity field model is set in the GFM_tide global variable in the module mdl_param.f03
       % 输入重力场模型的潮汐系统在模块mdl_param.f03的GFM_tide全局变量中被设置
       if (GFM_tide == "zero_tide")
           dCnm_solid1(2+1,0+1) = dCnm_solid1(2+1,0+1) - dC20_perm;
       end
 
       % Pole Tide  极潮
       % Solid Earth pole tide  地球固体极潮
       dC21_pse = 0.0;
       dS21_pse = 0.0;
       if (grav_effects.yml_solid_pole)
           [dC21_pse, dS21_pse] = tide_pole_se(mjd,xp,yp);  
       end
 
       % Ocean pole tide  海洋极潮
       dC21_poc = 0.0;
       dS21_poc = 0.0;

       if (grav_effects.yml_ocean_pole)
           [dC21_poc, dS21_poc] = tide_pole_oc(mjd,xp,yp);
       end

       % dCnm, dSnm arrays sum
       sz_tides = size(dCnm_solid1,1);
	   
       dCnm_tides = dCnm_solid1;
       dSnm_tides = dSnm_solid1;

       for i = 1:3
          for j = 1:3
              dCnm_tides(i,j) = dCnm_tides(i,j) + dCnm_solid2(i,j); 
              dSnm_tides(i,j) = dSnm_tides(i,j) + dSnm_solid2(i,j); 
          end
       end

      dCnm_tides(2+1,1+1) = dCnm_tides(2+1,1+1) + dC21_pse + dC21_poc; 
      dSnm_tides(2+1,1+1) = dSnm_tides(2+1,1+1) + dS21_pse + dS21_poc;

      % Acceleration cartesian components 加速度的笛卡尔坐标分量
      a_solidtides = [0.0, 0.0, 0.0];
      [ax, ay, az] = force_tides(rsat_itrf, GMearth, aEarth, sz_tides-1, sz_tides-1, dCnm_tides, dSnm_tides);
      a_solidtides(1) = ax;
      a_solidtides(2) = ay;
      a_solidtides(3) = az;
      
      % Ocean Tides  海洋潮汐影响
      if (grav_effects.yml_ocean_tides)
          [dCnm_ocean, dSnm_ocean] = tides_ocean(OCEAN_Nmax, OCEAN_Mmax, mjd, ut1_utc);
          % Acceleration cartesian components  加速度的笛卡尔坐标分量
      	  [ax,ay,az] = force_tides(rsat_itrf, GMearth, aEarth, OCEAN_Nmax, OCEAN_Mmax, dCnm_ocean, dSnm_ocean);
          a_ocean(1) = ax;
          a_ocean(2) = ay;
          a_ocean(3) = az;
      elseif (~ grav_effects.yml_ocean_tides)
          a_ocean = [0.0, 0.0, 0.0];
      end

      % Tides acceleration cartesian components(潮汐加速度的笛卡尔坐标分量): Overall effects(潮汐摄动整体的影响)
      a_tides = a_solidtides + a_ocean;

      % Acceleration vector transformation from terrestrial to inertial frame
      % 将加速度矢量从地球坐标系转换至惯性坐标系
      Ftides_icrf = TRS2CRS * a_tides;
   elseif (~ grav_effects.yml_tidal_effects)
      Ftides_icrf = [0.0, 0.0, 0.0];
   end
   % End of Tidal effects  潮汐影响的结束

   % Relativistic effects   相对论效应影响
   if (grav_effects.yml_rel_effects)
   
      % Satellite State Vector in GCRS  在GCRS坐标系中，卫星的状态矢量
      Zsat_GCRS = [rsat_icrf(1), rsat_icrf(2), rsat_icrf(3), vsat_icrf(1), vsat_icrf(2), vsat_icrf(3)]; 
      
      % Earth state vector w.r.t. Sun  地球的状态矢量
      Zearth_GCRS = -1.0 * [rSun(1), rSun(2), rSun(3), vSun(1), vSun(2), vSun(3)];
      
      % Schwarzschild terms
      a_Schwarzschild = rel_schwarzschild(Zsat_GCRS, Zearth_GCRS, GMearth, beta_ppn, gama_ppn, c_light);
      
      % Lense-Thirring effects	  
      %a_LenseThirring = rel_LenseThirring(Zsat_GCRS, Zearth_GCRS,  GMearth, gama_ppn, c_light);
      % FIXME: why call function if you are going to ignore the result?
      a_LenseThirring = [0.0, 0.0, 0.0];
      
      % de Sitter effect or geodesic precesssion 
      a_deSitter = rel_deSitter(Zsat_GCRS, Zearth_GCRS, GMearth, GM_sun, gama_ppn, c_light);

      Frelativity_icrf = a_Schwarzschild + a_LenseThirring + a_deSitter;

   elseif (~ grav_effects.yml_rel_effects)
       Frelativity_icrf = [0.0, 0.0, 0.0];
   end

   % Summary of gravitational effects  总的引力摄动影响
   SFgrav = Fgrav_icrf + Fplanets_icrf + Ftides_icrf + Frelativity_icrf; 
   % End of Gravitational Effects   引力摄动的结束

   
   % Attitude model :: GNSS Yaw-attitude models   姿态模型:: GNSS卫星的偏航姿态模型
   % Variables "satblk" and "BDSorbtype" are temporary manually configured 
   % in the main program file (main_pod.f03) through setting the global variables: 
   
   % GPS case: Satellite Block ID:        1=I, 2=II, 3=IIA, IIR=(4, 5), IIF=6
   satblk = SATblock_glb;
   
   % Beidou case: 'IGSO', 'MEO'
   BDSorbtype = BDSorbtype_glb;
   
   fmt_line = '(A1,I2.2)';
   READ (PRN, fmt_line , IOSTAT=ios) GNSSid, PRN_no;
   
   % Yaw-attitude model 偏航姿态模型
   PRN_GNSS = PRN;
   BLKsat = BLKTYP;					 
   [eclipsf, beta, Mangle, Yangle, eBX_nom, eBX_ecl] = attitude(mjd, rsat_icrf, vsat_icrf, rSun, PRN_GNSS, BLKsat);
   % ----------------------------------------------------------------------
   
   
   % ----------------------------------------------------------------------
   % shadow model :: A conical model is applied.   阴影函数模型
   %                 lambda = 1     : In SUN LIGHT AREA   在太阳光照中
   %                        = 0     : In UMBRA AREA (full eclipse)  本影
   %                 0 < lambda < 1 : In PENUMBRA AREA  半影
   % ----------------------------------------------------------------------
   [lambda, ~] = shadow(rsat_icrf, rSun, rMoon);
   
   % ----------------------------------------------------------------------
   % Non-Gravitational Effects  非引力摄动的影响
   % ----------------------------------------------------------------------

   % Solar Radiation  太阳光压摄动力
   if (non_grav_effects.yml_solar_radiation)
       % SRP model
       srpid = yml_apriori_srp;
       [fx, fy, fz] = force_srp(lambda, eBX_ecl, eclipsf, GMearth, GNSSid, srpid, rsat_icrf, vsat_icrf, rSun);
       Fsrp_icrf = [fx, fy, fz];
   else
       Fsrp_icrf = [0.0, 0.0, 0.0];
   end

   % Earth radiation pressure  地球反照辐射压摄动力
   if (non_grav_effects.yml_earth_radiation)
       [fx, fy, fz] = force_erp(mjd, rsat_icrf, vsat_icrf, rSun);
       Ferp_icrf = [fx, fy, fz];
   else
       Ferp_icrf = [0.0, 0.0, 0.0];
   end
   
   % Antenna thrust effect  天线辐射推力摄动
   if (non_grav_effects.yml_antenna_thrust)
       [fx, fy, fz] = force_ant(rsat_icrf);
       Fant_icrf = [fx, fy, fz];
   else
       Fant_icrf = [0.0, 0.0, 0.0];
   end

   % Summary of non-gravitational effects   总的非引力摄动的影响
   SFnongrav = Fsrp_icrf + Ferp_icrf + Fant_icrf;

   % End of non-Gravitational Effects  非引力摄动影响的结束

   % Empirical forces  经验力
   if (yml_EMP_mode)
       Frame_EmpiricalForces = Frame_EmpiricalForces_glb;      
       Yawangle = Yangle(2);
       [SFemp, ~, ~, ~] = pd_empirical(rsat_icrf, vsat_icrf, GMearth, Yawangle,Frame_EmpiricalForces);
   else
       SFemp = [0.0, 0.0, 0.0];
   end

   % Pseudo-stochastic pulses  伪随机脉冲
   if (yml_pulses)
       [SFpulses, ~, ~, ~] = pulses_force(rsat_icrf, vsat_icrf, mjd, t_sec, integr_stage); 
   else
       SFpulses = [0.0, 0.0, 0.0];
   end

   % Acceleration satellites of the force model   卫星总的加速度
   SF = SFgrav + SFnongrav + SFemp + SFpulses;
   SFx = SF(1); 
   SFy = SF(2);
   SFz = SF(3);

end