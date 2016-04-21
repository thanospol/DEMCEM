function [coef,coefm] = coefficients_f1_f1(r1,r2,r3,r4,ko)
%% coefficients_f1_f1

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
% r1,r2,r3, r4 = point vectors of the triangular element's vertices
% Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
% Inner triangle Q:(rq1,rq2,rq3)=(r2,r1,r4)
% ko = wavenumber

% OUTPUT DATA
% coef   
% coefm
%

c  = zeros(1);
cm = zeros(1);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
j = sqrt(-1);
% 
c(1) = cc(1) * dd(2) * bb(3) / 0.3e1 - cc(1) * dd(3) * bb(2) / 0.3e1 + cc(2) * dd(3) * bb(1) / 0.3e1 - cc(2) * dd(1) * bb(3) / 0.3e1 + cc(3) * dd(1) * bb(2) / 0.3e1 - cc(3) * dd(2) * bb(1) / 0.3e1;
%
t1 = (j ^ 2);
t4 = (ko ^ 2);
t7 = 3 * c(1) / t1 / t4;
%
coef(1) = -t7;
coef(2) = 3 * c(1) / j / ko;
coef(3) = c(1);
coef(4) = t7;
%
cm(1) = cc(1) * dd(2) * bb(3) / 0.3e1 - cc(1) * dd(3) * bb(2) / 0.3e1 + cc(2) * dd(3) * bb(1) / 0.3e1 - cc(2) * dd(1) * bb(3) / 0.3e1 + cc(3) * dd(1) * bb(2) / 0.3e1 - cc(3) * dd(2) * bb(1) / 0.3e1;
%
t1 = (j ^ 2);
t4 = (ko ^ 2);
t7 = 3 * cm(1) / t1 / t4;
%
coefm(1) = -t7;
coefm(2) = 3 * cm(1) / j / ko;
coefm(3) = cm(1);
coefm(4) = t7;