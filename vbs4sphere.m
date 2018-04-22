function [T, G1, G0, Ig1, Ig0, Stat, Qn] = vbs4sphere(P, tmax, sw)
% =====================================================================
%  Diffusion on 3-sphere, Single-target task, with orthant bias.
%  Drift variability, pools within 16 orthants.
%  Separate across-trial drift variability parameters.
%  Returns target-present/target-absent responses.
%
%  [T, G1, G0, Ig1, Ig0, Stat, Qn] = vbs4sphere(P, tmax, sw)
%        P = [v1, v2, v3, v4, eta1, eta2, eta3, eta4, sigma, a, Ter, st, c]
%  Dependencies:  (same C-code)
%       vbs4sphere --> vbd4orth300 --> ved4circle300
%  Building C-code:
%            mex vbd4orth300.c -lgsl -lgslcblas  -lm 
%  Stat = [P1, P2, M1, M2]
%  sw: 0 or omitted, no Stat or Qn
%      : 1 Calculate Stat and Qn
%      : 2 Calculate Stat Qn and execution time
%
%  Orthonant: target +, nontarget -, response rule: 1 - one; 0 - none
%  Alternating differential form property of mu4 - NOT simple factorial.
%
%  (Corrected)
%  1 + + + +  1
%  2 + + - +  1
%  3 + + - -  1
%  4 + + + -  1
%  5 + - + +  1
%  6 + - - +  1
%  7 + - - -  1
%  8 + - + -  1
%  9 - + + +  1
% 10 - + - +  1
% 11 - + - -  1
% 12 - + + -  1
% 13 - - + +  1
% 14 - - - +  1
% 15 - - - -  0
% 16 - - + -  1
% =======================================================================

name = 'VBS4SPHERE: ';
% Set me ! must agree with C-mex spacing.
% -------------------------------------- 
NS = 300;
north = 16; % Number of orthants
NP = 13;
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
G1 = zeros(1, NS);
G0 = zeros(1, NS);
Ig1 = zeros(1, NS);
Ig0 = zeros(1, NS);

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
c = P(13);


% Corrected
T1p = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,16];  % Target-present orthant.
T1a = [15]; % Target-absent orthant
P4 = [P(1:10), P(13)];
% Use biased orthants
[T, Gto, Fto, Po, Mo] = vbd4orth300x(P4, tmax, 5);
%plot(T, Gto(T1p, :), '--', T, Gto(T1a, :));
%[[1:16]',Po, Mo]
%pause

% Summing OK because mass-scaling in hitting angle variable.
for o = 1:north
    if any(T1p == o) % Orthant maps to target-present response
         G1 = G1 + Gto(o, :);
         Ig1 = Ig1 + Fto(o, :);
         Stat(1) = Stat(1) + Po(o);
         Stat(3) = Stat(3) + Po(o) * Mo(o);
    else 
         G0 = G0 + Gto(o, :);
         Ig0 = Ig0 + Fto(o, :);
         Stat(2) = Stat(2) + Po(o);
         Stat(4) = Stat(4) + Po(o) * Mo(o);
    end;
    %[o, Stat]
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
   G1 = conv(G1, fe); 
   G0 = conv(G0, fe);
   G1 = G1(1:n);
   G0 = G0(1:n);
   Ig1 = conv(Ig1, fe); 
   Ig0 = conv(Ig0, fe);
   Ig1 = Ig1(1:n);
   Ig0 = Ig0(1:n);
end;


% Calculate quantiles.
if sw
   Qf5 = [.1, .3, .5, .7, .9]; % Summary quantiles (can be changed).

   % Normalise distributions, select subrange spanning quantiles.
   Fd = Ig1 / max(Ig1); 
   Fs = Ig0 / max(Ig0);
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
