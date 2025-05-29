function [fc,P,flag] = fuel_model_Compact_simplified(v,a,g,project)
%FUEL_MODEL_COMPACT_SIMPLIFIED
%   Version 3.1
%   (C) 2023/03/10 by CIRCLES project energy team

%   Simplified fuel consumption model for vehicle:
%   Compact Sedan, a small sedan of mass 1543 kg.
%   The function outputs the instaneous fuel consumption rate (and
%   equivalent power and feasibility state), given instaneous velocity,
%   equivalent power and feasibility state), given instaneous velocity,
%   acceleration, and road grade.
%   Inputs may be given as vectors or matrices, resulting in outputs of
%   the same size.

%   Inputs:
%   v = vehicle velocity (in m/s)
%   a = vehicle acceleration (in m/s^2)
%   g = road grade (in radians)
%   project = true for projecting infeasible acceleration to maximum feasible acceleration (boolean)
%   [Note: a 3% downhill slope is to be inputted as g = atan(-0.03);
%   however, for any realistic slopes, the atan can be omitted]

%   Outputs:
%   fc = fuel consumption rate (in grams/sec)
%   P = equivalent power (in KW)
%   flag = 0 if no warning
%   flag = 1 if input infeasible (impossible to realize by the vehicle),
%   flag = 2 if velocity negative
%   [Note: The model extends gracefully outside of its feasibility
%   region: negative velocities are treated like zero; and non-realizable
%   (v,a,g)-requests are assigned fc and P values that simply extrapolate
%   the fitted polynomials if project = false.
%========================================================================
% Vehicle-specific model parameters
%========================================================================
fc_idle =0.0972;
b1 = 3.3605;
b2 = 41.6037;
b3 =  0.00021189;
b4 =  8.9362;
b5 = 3.9757;
b6 = 0.24476;
C0 = 0.15918;
C1 = 0.013463;
C2 = 0;
C3 = 3.1889e-05;
p0 = 0.047828;
p1 = 0.086975;
p2 = 6.825e-08;
q0 = 0.0025557;
q1 = 0.019099;
z0 = 0.13285;
z1 = 0.77984;
z2 = 0.0019733;
vc = 5.04;
beta0 = 0.0972;
a0 = -0.26981;
a1 = -0.0023996;
a2 = -9.0623;
a3 = -0.00029215;
a4 = -0.011899;
%========================================================================
% Other constants
%========================================================================
gs2kW = 42.36;                  % grams/sec to kW conversion factor:
%========================================================================
% Checks of inputs
%========================================================================
 if nargin<2, a = v*0; end       % if no accel given, assume cruising
 if nargin<3, g = v*0; end       % if no grade given, assume flat road
 if ~all(size(v)==size(a)&size(a)==size(g))
     error('Inputs must be of identical size.')
 end
%========================================================================
% Evaluate model
%========================================================================
% boundary of fitted feasibility region
 ma = min(b1,b2./max(v,1e-12)-b3*v.^2) - min(b4,b5+b6.*v).*g; 
 flag = zeros(size(v));         % initialize flag with zeros
 flag = flag + 2*(v<0) +...     % assign 2 where v<0,
        (v>=0).*(a>ma);         % and 1 where v>=0 but input infeasible

 if project, a = min(a,ma); end % projection of infeasible accel to max feasible accel
 v = max(v,0);                  % treat negative velocities as zero

% acceleration at which q(v) is min
aplus = max(a,-(p0 + p1*v +p2*v.^2)./(2*(q0 + q1.*max(v,1e-12))));

% Evaluate piecewise polyomial fuel consumption rate formula 
 fc = C0 + C1*v + C2*v.^2 + C3*v.^3 +... % cruising terms
     p0*a + p1*a.*v + p2*a.*v.^2 +...    % linear acceleration terms
     q0*aplus.^2 + q1*aplus.^2.*v +...   % quadrative acceleration terms
     z0*g + z1*g.*v + z2*g.*v.^2;        % road grade terms

 fc_min = (v<=vc)*beta0;        % minimum possible fc (fuel cut vel. vc)
 fc = max(fc,fc_min);           % cap fitted fc from below
 a_brake = a0 + a1.*v + a2.*g +...
           a3.*v.^2 + a4*v.*g;  % acceleration at fuel cut
 fc (v>vc & a<=a_brake)=0; 

 fc(v<.1&abs(a)<.01) = fc_idle; % set fc to idle for v<.1m/s&|a|<.01m/s^2

 P = fc*gs2kW;                  % calculate equivalent power

