function [T, Gd, Gs, Igd, Igs, Stat, Qn] = ved4sphere(P, tmax, sw)
% =====================================================================
%  Diffusion on 3-sphere, with drift variability, pools within 16 orthants.
%  Separate across-trial drift variability parameters.
%  NEW CONVENTION: returns double-single not correct-error.
%
%  [T, Gd, Gs, Igd, Igs, Stat, Qn] = vd4sphere(P, tmax, sw)
%        P = [v1, v2, v3, v4, eta1, eta2, eta3, eta4, sigma, a, Ter, st]
%  Dependencies:
%       ved4sphere --> ved4orth300 --> ved4circle300
%  Building C-code:
%            mex ved4orth300.c -lgsl -lgslcblas  -lm 
%  Stat = [P1, P2, M1, M2]
%  sw: 0 or omitted, no Stat or Qn
%      : 1 Calculate Stat and Qn
%      : 2 Calculate Stat Qn and execution time
%
%  Orthonant: target +, nontarget -, response rule: S - single; D - double
%  Alternating differential form property of mu4 - NOT simple factorial.
%
%  (Corrected)
%  1 + + + +  D
%  2 + + - +  D
%  3 + + - -  D
%  4 + + + -  D
%  5 + - + +  D
%  6 + - - +  D
%  7 + - - -  S
%  8 + - + -  D
%  9 - + + +  D
% 10 - + - +  D
% 11 - + - -  S
% 12 - + + -  D
% 13 - - + +  D
% 14 - - - +  S
% 15 - - - -  S
% 16 - - + -  S
% =======================================================================

name = 'VED4SPHERE: ';
% Set me ! must agree with C-mex spacing.
% -------------------------------------- 
NS = 300;
north = 16; % Number of orthants
NP = 12;
badix = 5;

if nargin < 3
   sw = 0;
   trace = 0;
end;
if sw == 2
   trace = 1;
   t0 = clock;
else
   trace = 0;
end;
if length(P) ~= NP
   disp([name, ' Incorrect length parameter vector, exiting...']);
   length(P)
   return;
end;

h = tmax / (NS - 1);
      
epsilon = .0001;
Stat = zeros(1,4);
Gd = zeros(1, NS);
Gs = zeros(1, NS);
Igd = zeros(1, NS);
Igs = zeros(1, NS);

v1 = P(1);
v2 = P(2);
v3 = P(3);
v4 = P(4);
eta1 = P(5);
eta2 = P(6);
eta3 = P(7);
eta4 = P(8);
sigma = P(9);
a = P(10);
Ter = P(11);
st = P(12);

% Corrected
S = [7,11,14,15,16];             % Single target orthants
D = [1,2,3,4,5,6,8,9,10,12,13];  % Double target orthants


P4 = P(1:10);
[T, Gto, Fto, Po, Mo] = ved4orth300(P4, tmax, 5);
%[[1:16]',Po, Mo]

% Summing OK because mass-scaling in hitting angle variable.
for o = 1:north
    if any(D == o) % Orthant maps to double target response
         Gd = Gd + Gto(o, :);
         Igd = Igd + Fto(o, :);
         Stat(1) = Stat(1) + Po(o);
         Stat(3) = Stat(3) + Po(o) * Mo(o);
    else 
         Gs = Gs + Gto(o, :);
         Igs = Igs + Fto(o, :);
         Stat(2) = Stat(2) + Po(o);
         Stat(4) = Stat(4) + Po(o) * Mo(o);
    end;
end
Stat(3:4) = (Stat(3:4) / (Stat(1:2) + epsilon)) + Ter + st;

T = T + Ter + st / 2;
% --------------------
% Convolve with Ter.
% --------------------
if st > 2 * h
   m = round(st/h);
   n = length(T);
   fe = ones(1, m) / m;
   Gd = conv(Gd, fe); 
   Gs = conv(Gs, fe);
   Gd = Gd(1:n);
   Gs = Gs(1:n);
   Igd = conv(Igd, fe); 
   Igs = conv(Igs, fe);
   Igd = Igd(1:n);
   Igs = Igs(1:n);
end;


% Calculate quantiles.
if sw
   Qf5 = [.1, .3, .5, .7, .9]; % Summary quantiles (can be changed).

   % Normalise distributions, select subrange spanning quantiles.
   Fd = Igd / max(Igd); 
   Fs = Igs / max(Igs);
   Dummy5 = [.001,.002,.003,.004,.005]';
   
   % New code to deal with possible nonmonotonicities.
   % ------------------------------------------------
   fdmx = find(Fd == max(Fd));
   fsmx = find(Fs == max(Fs));
   Fd = Fd(1:fdmx);
   Fs = Fs(1:fsmx);
   % -------------------------------------------------     

   Id = (Fd >= .025 & Fd <= .975);  
   if min(diff(Fd(Id))) <= 0 
        Qd = Dummy5;
        disp('Cannot compute F1 quantiles.');    
   else
        Qd=interp1(Fd(Id)', T(Id)', Qf5);
   end;

   Is = (Fs >= .025 & Fs <= .975);  
   if min(diff(Fs(Is))) <= 0 
        Qs = Dummy5;
        disp('Cannot compute F2 quantiles.');    
   else
        Qs=interp1(Fs(Is)', T(Is)', Qf5);
   end;
   Qn = [Qd',Qs'];
else
   Qn = [];
end;

if sw == 2
     fprintf('ET = %10.3f \n',etime(clock, t0));
end;
