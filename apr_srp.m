function [X_SIDE, Z_SIDE, A_SOLAR, F0] = apr_srp(GNSSid, BLKTYP)
% ! Purpose:
% ! TBD
% ! ----------------------------------------------------------------------
% ! Input arguments:
% ! - GNSSid       : id of satellite constellation  
% ! - blktyp       : satellite block type
% ! 
% ! Output arguments:
% ! - X_SIDE       : area of X-sdie in body-fixed frame
% ! - Z_SIDE       : area of Z-side in body-fixed frame
% ! - A_SOLAR      : area of solar panel
% ! - F0           : net force (computed from cannoball) in D direction
% ! ----------------------------------------------------------------------

  % default init
  Z_SIDE = 0.d0;
  X_SIDE = 0.d0;
  A_SOLAR = 0.d0;
  F0 = 0.d0;
  
  if (GNSSid == 'G')
      % GPS constellation
      if (TRIM(BLKTYP)=="GPS-I")  % I
          Z_SIDE = 3.020D0;
          X_SIDE = 1.728D0;
          A_SOLAR= 6.053D0;
          F0 = 4.54e-5;
      elseif (TRIM(BLKTYP)=="GPS-II" || TRIM(BLKTYP)=="GPS-IIA")  % II and IIA
          Z_SIDE = 2.881D0;
          X_SIDE = 2.893D0;
          A_SOLAR= 11.871D0;
          F0 = 8.695e-5;
      elseif (TRIM(BLKTYP)=="GPS-IIF")  % IIF
           Z_SIDE = 5.05D0;
           X_SIDE = 4.55D0;
           A_SOLAR= 22.25D0;
           F0 = 16.7e-5;
  
      elseif (TRIM(BLKTYP)=="GPS-IIR" || TRIM(BLKTYP)=="GPS-IIR-A" || ...
                    TRIM(BLKTYP)=="GPS-IIR-B" || TRIM(BLKTYP)=="GPS-IIR-M")  % IIR
           Z_SIDE = 4.25D0;
           X_SIDE = 4.11D0;
           A_SOLAR= 13.92D0;
           F0 = 11.15e-5;
  
      elseif (TRIM(BLKTYP)=="GPS-IIIA")  % III
           Z_SIDE = 4.38D0;
           X_SIDE = 6.05D0;
           A_SOLAR= 22.25D0;
           F0 = 11.0e-5;
  
      else
           fprintf('%s %d %s\n',"apr_srp - Unknown block type: ", BLKTYP, "GNSSid = 'G'");
           return;     
      end
  
  elseif (GNSSid == 'R')   % GLONASS constellation
     if(TRIM(BLKTYP)=="GLO" || TRIM(BLKTYP)=="GLO-M" || TRIM(BLKTYP)=="GLO-M+" || ...
              TRIM(BLKTYP)=="GLO-K1A" || TRIM(BLKTYP)=="GLO-K1B")
        Z_SIDE = 1.6620D0;
        X_SIDE = 4.200D0;
        A_SOLAR= 23.616D0;
        if(TRIM(BLKTYP)=="GLO-K1A" || TRIM(BLKTYP)=="GLO-K1B")  % GLONASS-K
            F0 = 10.0e-5;
        end
  
        if(TRIM(BLKTYP)=="GLO" || TRIM(BLKTYP)=="GLO-M" || TRIM(BLKTYP)=="GLO-M+")  % GLONASS-M
            F0 = 20.9e-5;
        end

     else
        fprintf('%s %d %s\n',"apr_srp - Unknown block type: ", BLKTYP, "GNSSid = 'R'");
        return;    
     end
  
  elseif (GNSSid == 'E')   % GALILEO constellation

      if (TRIM(BLKTYP)=="GAL-1"||TRIM(BLKTYP)=="GAL-2")
          Z_SIDE = 3.002D0;
          X_SIDE = 1.323D0;
          A_SOLAR= 11.0D0;
          F0 = 8.35e-5;
         
      else
          fprintf('%s %d %s\n',"apr_srp - Unknown block type: ", BLKTYP, "GNSSid = 'E'");
          return;    
      end
  
  elseif (GNSSid == 'C')   % BDS constellation
      Z_SIDE = 3.96D0;
      X_SIDE = 4.5D0;
      A_SOLAR= 22.44D0;

      % BDS GEO
      if(TRIM(BLKTYP)=="BDS-2G")
          F0 = 50.1e-5;
      end

      if(TRIM(BLKTYP)=="BDS-3G")
          F0 = 50.1e-5;
      end

      % BDS MEO
      if(TRIM(BLKTYP)=="BDS-2M") 
          F0 = 8.35e-5;
      end

      if(TRIM(BLKTYP)=="BDS-3M-CAST")
          F0 = 8.35e-5;
      end

      if(TRIM(BLKTYP)=="BDS-3M-SECM-A")
          F0 = 8.35e-5;
      end

      if(TRIM(BLKTYP)=="BDS-3M-SECM-B")
          F0 = 8.35e-5;
      end
      % BDS IGSO
      if(TRIM(BLKTYP)=="BDS-2I")
          F0 = 17.0e-5;
      end

      if(TRIM(BLKTYP)=="BDS-3I")
          F0 = 50.1e-5;
      end

      if (BLKTYP(1:3)~="BDS")
          fprintf('%s %d %s\n',"apr_srp - Unknown block type: ", BLKTYP, "GNSSid = 'C'");
          return;    
      end
  elseif (GNSSid == 'J')  % QZSS constellation

      if (TRIM(BLKTYP)=="QZS-1")  % QZSS-1
          Z_SIDE = 5.6D0;
          X_SIDE = 9.0D0;
          A_SOLAR= 45.0D0;
          F0 = 35.0e-5;
      elseif (TRIM(BLKTYP)=="QZS-2I") % QZSS-2I         
          Z_SIDE = 5.6D0;
          X_SIDE = 10.1D0;
          A_SOLAR= 29.8D0;
          F0 = 25.0e-5;
  
      elseif (TRIM(BLKTYP)=="QZS-2G") % ! QZSS-2G
          Z_SIDE = 5.6D0;
          X_SIDE = 10.1D0;
          A_SOLAR= 29.8D0;
          F0 = 25.0e-5;
  
      else
         fprintf('%s %d %s\n',"apr_srp - Unknown block type: ", BLKTYP, "GNSSid = 'J'");
         return;
      end
  end

end