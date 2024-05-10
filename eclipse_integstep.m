function [integstep_flag, integstep_initial, integstep_reduced] = eclipse_integstep(EQMfname, VEQfname, mjd, r_sat, v_sat)
% ----------------------------------------------------------------------
% Purpose:
%  Reduce the step of the numerical integration method applied for orbit 
%  determination during eclipse seasons 
% ----------------------------------------------------------------------
% Input arguments:
% - EQMfname: 	Input configuration file name for the orbit parameterization 
% - VEQfname: 	Input configuration file name for the orbit parameterization 
% - mjd:		Modified Julian Day (MJD) in Terrestrial Time (including the fraction of the day)
% - r_sat: 		Satellite position vector (m) in ICRF
% - v_sat: 		Satellite velocity vector (m/sec) in ICRF
%
% Output arguments:
% - integstep_flag: 	Flag regarding the change of the orbit integration step 
% - integstep_initial:	Initial value of the orbit integration method 
% - integstep_reduced:	Reduced value of the orbit integration method 
% ----------------------------------------------------------------------

  integstep_flag = false;
  found = false;
  
  reductions(1) = 100;
  reductions(2) = 75;
  reductions(3) = 50;
  reductions(4) = 25;
  reductions(5) = 10;
  reductions(6) = 5;
  
  % Sun (NTARG) Cartesian coordinates w.r.t. Earth (NCTR)
  [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, r_sun, ~] = JPL_Eph_DE440(mjd);
  
  % Beta angle (degrees)
  beta = beta_angle(r_sat, v_sat, r_sun);
  fprintf('%s %d %d %d %s %s %d\n',"day of year", YR, DOY, PRN, "  ", "beta", beta);
  
  % criteria for changing orbit integration stepsize
  if (abs(beta) <= 14)
      integstep_initial = integstep;
      integstep_int = round(integstep);
      if (integstep_initial >= 100.0)
          for i = 1:6
              integstep_reduced = 1.0d0 * reductions(i);
              integstep_int2 = integstep_int/reductions(i);
              integstep_int2 = integstep_int2*reductions(i);
              if (integstep_int == integstep_int2)
                  integstep_flag = true;
                  fprintf('%s %d %s %d %s\n',"eclipsing satellite ", PRN, ", reduced stepsize of ", reductions(i), " chosen");
                  break;
              end
          end
      end
  
      if (integstep_flag && ~ yaml_found)
          write (fname_id, *) "_INT";
          param_id = "integrator_step";
          write (param_value, *) integstep_reduced;
          Call write_prmfile (EQMfname, fname_id, param_id, param_value)
          Call write_prmfile (VEQfname, fname_id, param_id, param_value)
      end
  end

  for i = 1:prn_override_count
      if (yml_prn_overrides(i).name == TRIM(PRN))
          if (integstep_flag)
              yml_prn_overrides(i).integ.veq_enabled = true;
              yml_prn_overrides(i).integ.eqm_enabled = true;
              yml_prn_overrides(i).integ.veq_stepsize = integstep_reduced;
              yml_prn_overrides(i).integ.eqm_stepsize = integstep_reduced;
          else
              yml_prn_overrides(i).integ.veq_enabled = false;
              yml_prn_overrides(i).integ.eqm_enabled = false;
          end
          found = true;
      end
  end

  if (~ found && integstep_flag)
      call new_prn_override(PRN);
      i = prn_override_count;
      yml_prn_overrides(i).integ.veq_enabled = true;
      yml_prn_overrides(i).integ.eqm_enabled = true;
      yml_prn_overrides(i).integ.veq_stepsize = integstep_reduced;
      yml_prn_overrides(i).integ.eqm_stepsize = integstep_reduced;
  end
end