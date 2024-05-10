function [ECOMNUM, srpcoef, yml_SRP_mode, SRP] = ecom_init(SRP_model)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明
    global const;

    if (SRP_model == const.ECOM1)
       ECOMNUM = 5;
       srpcoef(ECOMNUM) = 0.0;
       yml_SRP_mode = SRP_model;
       SRP = set_srpparam(SRP_model);
       disp('ECOM1 SRP MODEL IS ACTIVATED');
   elseif (SRP_model == const.ECOM2)
       ECOMNUM = 9;
       srpcoef(ECOMNUM) = 0.0;
       yml_SRP_mode = SRP_model;
       SRP = set_srpparam(SRP_model);
       disp('ECOM2 SRP MODEL IS ACTIVATED');
   elseif (SRP_model == const.SBOXW)
       ECOMNUM = 7;
       srpcoef(ECOMNUM) = 0.0;
       yml_SRP_mode = SRP_model;
       disp('SIMPLE BOX WING IS ACTIVATED');
   
   elseif (SRP_model == const.BERNE)
       ECOMNUM = 9;
       srpcoef(ECOMNUM) = 0.0;
       yml_SRP_mode = SRP_model;
       SRP = set_srpparam(SRP_model);
       disp('BERNE SRP MODEL IS ACTIVATED');
   elseif(SRP_model == const.ECOM_HYBRID)
       ECOMNUM = 13;
       srpcoef(ECOMNUM) = 0.0;
       yml_SRP_mode = SRP_model;
       SRP = set_srpparam(SRP_model);
       disp('ECOM1+ECOM2 HYBRID MODEL IS ACTIVATED');
   end
end