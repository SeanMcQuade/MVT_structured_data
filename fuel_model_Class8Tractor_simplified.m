function [fc,P,flag] = fuel_model_Class8Tractor_simplified(v,a,g,project)
%FUEL_MODEL_CLASS8TRACTOR_SIMPLIFIED
%   Version 3.1
%   (C) 2023/03/10 by CIRCLES project energy team

%   Simplified fuel consumption model for vehicle:
%   Class8Tractor, a truck with manual transmission and a trailer and total mass 25104 kg.
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
fc_idle =0.2384;
b1 = 2.4229;
b2 = 8.4463;
b3 =  0.00026269;
b4 =  9.7395;
b5 = 8.6171;
b6 = 0.15762;
C0 = 0.59448;
C1 = 0.082609;
C2 = 0;
C3 = 0.00027278;
p0 = 0.20476;
p1 = 1.1962;
p2 = 0.019119;
q0 = 0;
q1 = 0.14424;
z0 = 0.88147;
z1 = 11.1899;
z2 = 0.18836;
h0 = 0.49109;
h1 = 0.025152;
%========================================================================
% Other constants
%========================================================================
gs2kW =  42.47;                 % grams/sec to kW conversion factor:
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

 H = h0 + h1*v;                      % minimum possible fitted fc
 fc = max(fc,H);                     % cap fitted fc from below

 fc(v<.1&abs(a)<.01) = fc_idle; % set fc to idle for v<.1m/s&|a|<.01m/s^2

 P = fc*gs2kW;                  % calculate equivalent power

