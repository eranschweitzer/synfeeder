function b = b_from_c(c,omega)
%%% takes capacitance, c, in micro-Farad and returns susceptance, b, with
%%% unit Siemens
%%% 
%%% INPUTS:
%%%         c: capacitance [uF]
%%%         omega: angular frequence [rad/sec] i.e. 2*pi*f, with f in Hz
%%% OUTPUTS:
%%%         b: susceptance [S]

b = omega*c*1e-6;