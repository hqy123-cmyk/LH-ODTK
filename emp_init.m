function [EMPNUM, empcoef, yml_EMP_mode, EMP] = emp_init(EMP_model, EMP_cpr_NO)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明
   yml_EMP_mode = EMP_model;
   if (yml_EMP_mode == 1)
       EMPNUM = 3;
       disp('EMPIRICAL MODEL IS ACTIVATED');
       % Bias accelerations per radial, along-track and cross-track directions
       EMP.bias_R = 1;
       EMP.bias_T = 1;
       EMP.bias_N = 1;
       % Cycle-per-revolution accelerations per radial, along-track and cross-track directions
       if (EMP_cpr_NO)
           EMPNUM = EMPNUM + 6;
           EMP.cpr_R  = 1;
           EMP.cpr_T  = 1;
           EMP.cpr_N  = 1;
       else
           EMP.cpr_R  = 0;
           EMP.cpr_T  = 0;
           EMP.cpr_N  = 0;
       end
       empcoef(EMPNUM) = 0.0;
   end
end