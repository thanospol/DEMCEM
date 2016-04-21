function [coef,coefm] = coefficients_const(ko)
%% coefficients_const

%  Licensing: This code is distributed under the GNU LGPL license. 

%  Modified:  20 September 2011

%  Author:    Athanasios Polimeridis

% References

% A. G. Polimeridis and T. V. Yioultsis, “On the direct evaluation of weakly singular
% integrals in Galerkin mixed potential integral equation formulations,” IEEE Trans.
% Antennas Propag., vol. 56, no. 9, pp. 3011-3019, Sep. 2008.

% A. G. Polimeridis and J. R. Mosig, “Complete semi-analytical treatment of weakly
% singular integrals on planar triangles via the direct evaluation method,” Int. J.
% Numerical Methods Eng., vol. 83, pp. 1625-1650, 2010.

% A. G. Polimeridis, J. M. Tamayo, J. M. Rius and J. R. Mosig, “Fast and accurate
% computation of hyper-singular integrals in Galerkin surface integral equation
% formulations via the direct evaluation method,” IEEE Trans.
% Antennas Propag., vol. 59, no. 6, pp. 2329-2340, Jun. 2011.

% A. G. Polimeridis and J. R. Mosig, “On the direct evaluation of surface integral
% equation impedance matrix elements involving point singularities,” IEEE Antennas
% Wireless Propag. Lett., vol. 10, pp. 599-602, 2011.

% INPUT DATA
% ko = wavenumber

% OUTPUT DATA
% coef   
% coefm
%

coef  = zeros(3);
coefm = zeros(3);
%
j = sqrt(-1);
% 
c(1) = 1;
%
t5 = j ^ 2;
t8 = ko ^ 2;
t10 = c(1) / t5 / t8;
%
coef(1) = -c(1) / j / ko;
coef(2) = -t10;
coef(3) = t10;
%
cm(1) = 1;
%
t1 = j ^ 2;
t4 = ko ^ 2;
t6 = cm(1) / t1 / t4;
%
coefm(1) = -t6;
coefm(2) = t6;
coefm(3) = -cm(1) / j / ko;