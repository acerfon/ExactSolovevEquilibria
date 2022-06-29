clear all

% The construction of the flux functions that we use here is presented in
% detail in A.J. Cerfon and J.P. Freidberg, ``One size fits all" analytic
% solutions to the Grad-Shafranov equation, Physics of Plasmas 17, 032502 (2010)


% Copyright (C) 2010-2012: Antoine Cerfon
% Contact: cerfon@cims.nyu.edu
% 
% This program is free software; you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation; either version 2 of the License, or 
% (at your option) any later version.  This program is distributed in 
% the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
% PARTICULAR PURPOSE.  See the GNU General Public License for more 
% details. You should have received a copy of the GNU General Public 
% License along with this program; 
% if not, see <http://www.gnu.org/licenses/>.
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Specify equilibrium parameters of interest (in this example, we took
%   the NSTX spherical torus as an example)
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% NSTX PARAMETERS %%%%%%%%%%%%%%
epsilon = 0.68/0.85;%Inverse aspect ratio
kappa = 2;%Elongation
delta  = 0.35;%Triangularity
xsep = 0.5/0.85;%x-location of the separatrix
ysep = -1.5/0.85;%y-location of the separatrix
betaparam = -0.1435;%betaparam is the parameter determining the beta regime of interest (Note: it's called A in the article)
D=[0,-0.015,-0.04,-0.075,-0.115,-0.16,-0.2,-0.23,-0.25];%Values of the contours for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = asin(delta);%alpha as defined in the article
slope1 = 0;%outer equatorial point slope
slope2 = 0;%inner equatorial point slope
curv1 = -(1+alpha)^2/(epsilon*kappa^2);%curvature at the outboard midplane
curv2 = -kappa/(epsilon*(cos(alpha))^2);%curvature at the top
curv3 = (1-alpha)^2/(epsilon*kappa^2);%curvature at the inboard midplane

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Construct the matrix A of the boundary conditions for the funtions
%   which are solutions to the homogeneous equation
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = [psi1(1+epsilon,0) psi2(1+epsilon,0) psi3(1+epsilon,0) psi4(1+epsilon,0) psi5(1+epsilon,0) psi6(1+epsilon,0) psi7(1+epsilon,0) psi8(1+epsilon,0) psi9(1+epsilon,0) psi10(1+epsilon,0) psi11(1+epsilon,0) psi12(1+epsilon,0) %outer equatorial point
    psi1(1-epsilon,0) psi2(1-epsilon,0) psi3(1-epsilon,0) psi4(1-epsilon,0) psi5(1-epsilon,0) psi6(1-epsilon,0) psi7(1-epsilon,0) psi8(1-epsilon,0) psi9(1-epsilon,0) psi10(1-epsilon,0) psi11(1-epsilon,0) psi12(1-epsilon,0) %inner equatorial point
    psi1(1-epsilon*delta,kappa*epsilon) psi2(1-epsilon*delta,kappa*epsilon) psi3(1-epsilon*delta,kappa*epsilon) psi4(1-epsilon*delta,kappa*epsilon) psi5(1-epsilon*delta,kappa*epsilon) psi6(1-epsilon*delta,kappa*epsilon) psi7(1-epsilon*delta,kappa*epsilon) psi8(1-epsilon*delta,kappa*epsilon) psi9(1-epsilon*delta,kappa*epsilon) psi10(1-epsilon*delta,kappa*epsilon) psi11(1-epsilon*delta,kappa*epsilon) psi12(1-epsilon*delta,kappa*epsilon) %upper high point
    psi1(xsep,ysep) psi2(xsep,ysep) psi3(xsep,ysep) psi4(xsep,ysep) psi5(xsep,ysep) psi6(xsep,ysep) psi7(xsep,ysep) psi8(xsep,ysep) psi9(xsep,ysep) psi10(xsep,ysep) psi11(xsep,ysep) psi12(xsep,ysep) %lower X point    
    slope1*psi1x(1+epsilon,0)+psi1y(1+epsilon,0) slope1*psi2x(1+epsilon,0)+psi2y(1+epsilon,0) slope1*psi3x(1+epsilon,0)+psi3y(1+epsilon,0) slope1*psi4x(1+epsilon,0)+psi4y(1+epsilon,0) slope1*psi5x(1+epsilon,0)+psi5y(1+epsilon,0) slope1*psi6x(1+epsilon,0)+psi6y(1+epsilon,0) slope1*psi7x(1+epsilon,0)+psi7y(1+epsilon,0) slope1*psi8x(1+epsilon,0)+psi8y(1+epsilon,0) slope1*psi9x(1+epsilon,0)+psi9y(1+epsilon,0) slope1*psi10x(1+epsilon,0)+psi10y(1+epsilon,0) slope1*psi11x(1+epsilon,0)+psi11y(1+epsilon,0) slope1*psi12x(1+epsilon,0)+psi12y(1+epsilon,0) %outer equatorial point slope
    slope2*psi1x(1-epsilon,0)+psi1y(1-epsilon,0) slope2*psi2x(1-epsilon,0)+psi2y(1-epsilon,0) slope2*psi3x(1-epsilon,0)+psi3y(1-epsilon,0) slope2*psi4x(1-epsilon,0)+psi4y(1-epsilon,0) slope2*psi5x(1-epsilon,0)+psi5y(1-epsilon,0) slope2*psi6x(1-epsilon,0)+psi6y(1-epsilon,0) slope2*psi7x(1-epsilon,0)+psi7y(1-epsilon,0) slope2*psi8x(1-epsilon,0)+psi8y(1-epsilon,0) slope2*psi9x(1-epsilon,0)+psi9y(1-epsilon,0) slope2*psi10x(1-epsilon,0)+psi10y(1-epsilon,0) slope2*psi11x(1-epsilon,0)+psi11y(1-epsilon,0) slope2*psi12x(1-epsilon,0)+psi12y(1-epsilon,0) %inner equatorial point slope
    psi1x(1-epsilon*delta,kappa*epsilon) psi2x(1-epsilon*delta,kappa*epsilon) psi3x(1-epsilon*delta,kappa*epsilon) psi4x(1-epsilon*delta,kappa*epsilon) psi5x(1-epsilon*delta,kappa*epsilon) psi6x(1-epsilon*delta,kappa*epsilon) psi7x(1-epsilon*delta,kappa*epsilon) psi8x(1-epsilon*delta,kappa*epsilon) psi9x(1-epsilon*delta,kappa*epsilon) psi10x(1-epsilon*delta,kappa*epsilon) psi11x(1-epsilon*delta,kappa*epsilon) psi12x(1-epsilon*delta,kappa*epsilon) %upper high point maximum
    psi1x(xsep,ysep) psi2x(xsep,ysep) psi3x(xsep,ysep) psi4x(xsep,ysep) psi5x(xsep,ysep) psi6x(xsep,ysep) psi7x(xsep,ysep) psi8x(xsep,ysep) psi9x(xsep,ysep) psi10x(xsep,ysep) psi11x(xsep,ysep) psi12x(xsep,ysep) %By = 0 at lower X-point
    psi1y(xsep,ysep) psi2y(xsep,ysep) psi3y(xsep,ysep) psi4y(xsep,ysep) psi5y(xsep,ysep) psi6y(xsep,ysep) psi7y(xsep,ysep) psi8y(xsep,ysep) psi9y(xsep,ysep) psi10y(xsep,ysep) psi11y(xsep,ysep) psi12y(xsep,ysep) %Bx = 0 at lower X-point
    curv1*psi1x(1+epsilon,0)+psi1yy(1+epsilon,0) curv1*psi2x(1+epsilon,0)+psi2yy(1+epsilon,0) curv1*psi3x(1+epsilon,0)+psi3yy(1+epsilon,0) curv1*psi4x(1+epsilon,0)+psi4yy(1+epsilon,0) curv1*psi5x(1+epsilon,0)+psi5yy(1+epsilon,0) curv1*psi6x(1+epsilon,0)+psi6yy(1+epsilon,0) curv1*psi7x(1+epsilon,0)+psi7yy(1+epsilon,0) curv1*psi8x(1+epsilon,0)+psi8yy(1+epsilon,0) curv1*psi9x(1+epsilon,0)+psi9yy(1+epsilon,0) curv1*psi10x(1+epsilon,0)+psi10yy(1+epsilon,0) curv1*psi11x(1+epsilon,0)+psi11yy(1+epsilon,0) curv1*psi12x(1+epsilon,0)+psi12yy(1+epsilon,0)%curvature condition at outer equatorial point
    curv3*psi1x(1-epsilon,0)+psi1yy(1-epsilon,0) curv3*psi2x(1-epsilon,0)+psi2yy(1-epsilon,0) curv3*psi3x(1-epsilon,0)+psi3yy(1-epsilon,0) curv3*psi4x(1-epsilon,0)+psi4yy(1-epsilon,0) curv3*psi5x(1-epsilon,0)+psi5yy(1-epsilon,0) curv3*psi6x(1-epsilon,0)+psi6yy(1-epsilon,0) curv3*psi7x(1-epsilon,0)+psi7yy(1-epsilon,0) curv3*psi8x(1-epsilon,0)+psi8yy(1-epsilon,0) curv3*psi9x(1-epsilon,0)+psi9yy(1-epsilon,0) curv3*psi10x(1-epsilon,0)+psi10yy(1-epsilon,0) curv3*psi11x(1-epsilon,0)+psi11yy(1-epsilon,0) curv3*psi12x(1-epsilon,0)+psi12yy(1-epsilon,0)%curvature condition at inner equatorial point
    curv2*psi1y(1-epsilon*delta,kappa*epsilon)+psi1xx(1-epsilon*delta,kappa*epsilon) curv2*psi2y(1-epsilon*delta,kappa*epsilon)+psi2xx(1-epsilon*delta,kappa*epsilon) curv2*psi3y(1-epsilon*delta,kappa*epsilon)+psi3xx(1-epsilon*delta,kappa*epsilon) curv2*psi4y(1-epsilon*delta,kappa*epsilon)+psi4xx(1-epsilon*delta,kappa*epsilon) curv2*psi5y(1-epsilon*delta,kappa*epsilon)+psi5xx(1-epsilon*delta,kappa*epsilon) curv2*psi6y(1-epsilon*delta,kappa*epsilon)+psi6xx(1-epsilon*delta,kappa*epsilon) curv2*psi7y(1-epsilon*delta,kappa*epsilon)+psi7xx(1-epsilon*delta,kappa*epsilon) curv2*psi8y(1-epsilon*delta,kappa*epsilon)+psi8xx(1-epsilon*delta,kappa*epsilon) curv2*psi9y(1-epsilon*delta,kappa*epsilon)+psi9xx(1-epsilon*delta,kappa*epsilon) curv2*psi10y(1-epsilon*delta,kappa*epsilon)+psi10xx(1-epsilon*delta,kappa*epsilon) curv2*psi11y(1-epsilon*delta,kappa*epsilon)+psi11xx(1-epsilon*delta,kappa*epsilon) curv2*psi12y(1-epsilon*delta,kappa*epsilon)+psi12xx(1-epsilon*delta,kappa*epsilon)];%curvature condition at top

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Construct the matrix B of the boundary conditions for the particular
%   solutions to the equation
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = -[betaparam*psipart1(1+epsilon,0)+(1-betaparam)*psipart2(1+epsilon,0) % outer equatorial point
    betaparam*psipart1(1-epsilon,0)+(1-betaparam)*psipart2(1-epsilon,0) % inner equatorial point
    betaparam*psipart1(1-epsilon*delta,kappa*epsilon)+(1-betaparam)*psipart2(1-epsilon*delta,kappa*epsilon) %upper high point
    betaparam*psipart1(xsep,ysep)+(1-betaparam)*psipart2(xsep,ysep) %lower X-point
    betaparam*(slope1*psipart1x(1+epsilon,0)+psipart1y(1+epsilon,0))+(1-betaparam)*(slope1*psipart2x(1+epsilon,0)+psipart2y(1+epsilon,0)) %outer equatorial point slope
    betaparam*(slope2*psipart1x(1-epsilon,0)+psipart1y(1-epsilon,0))+(1-betaparam)*(slope2*psipart2x(1-epsilon,0)+psipart2y(1-epsilon,0)) %inner equatorial point slope
    betaparam*psipart1x(1-epsilon*delta,kappa*epsilon)+(1-betaparam)*psipart2x(1-epsilon*delta,kappa*epsilon) %upper high point maximum
    betaparam*psipart1x(xsep,ysep)+(1-betaparam)*psipart2x(xsep,ysep) %By = 0 at lower X-point
    betaparam*psipart1y(xsep,ysep)+(1-betaparam)*psipart2y(xsep,ysep) %Bx = 0 at lower X-point
    betaparam*(curv1*psipart1x(1+epsilon,0)+psipart1yy(1+epsilon,0))+(1-betaparam)*(curv1*psipart2x(1+epsilon,0)+psipart2yy(1+epsilon,0))%curvature condition at outer equatorial point
    betaparam*(curv3*psipart1x(1-epsilon,0)+psipart1yy(1-epsilon,0))+(1-betaparam)*(curv3*psipart2x(1-epsilon,0)+psipart2yy(1-epsilon,0))%curvature condition at inner equatorial point
    betaparam*(curv2*psipart1y(1-epsilon*delta,kappa*epsilon)+psipart1xx(1-epsilon*delta,kappa*epsilon))+(1-betaparam)*(curv2*psipart2y(1-epsilon*delta,kappa*epsilon)+psipart2xx(1-epsilon*delta,kappa*epsilon))];%curvature condition at top

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Solve the linear system for the coefficients C of the general
%   solution to the equation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = A\B;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot the solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y] = meshgrid(0:.05:1+epsilon+0.5,ysep-0.075:.05:kappa*epsilon+0.075);
Z = psitot(X,Y,C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9),C(10),C(11),C(12),betaparam);
contour(X, Y, Z, D,'LineWidth',1.5);
line([0 0],[-1.5*kappa*epsilon 1.5*kappa*epsilon], 'Color','black','LineWidth',0.5,'LineStyle','--')
axis equal
xlim([0 1+epsilon+0.5])
ylim([ysep-0.075 kappa*epsilon+0.075])
colormap jet