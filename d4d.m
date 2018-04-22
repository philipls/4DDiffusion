function [G2, B, Pred] = d4d(Pvar, Pfix, Sel, Dmatrix, trace)
% ===================================================================
% Fit 4D diffusion model to Corbett & Smith (2017) double-target task 
% 20 drift rates, 10 x etas.
% July 5, 2017
% Usage:
%     [G2, B, Pred] = d4d(Pvar, Pfix, Sel, Dmatrix, trace);
%
%      P = [v1...20,  a,  eta1..eta5, eta6.eta10,   Ter   st]
%          [1...20,   21     22:26      27:31     32    33]
%  Parallel-friendly trace: second element is save file suffix
% ====================================================================
t0 = clock;
name = 'D4D:';
errmg1 = 'Incorrect number of parameters for model, exiting...';
errmg2 = 'Incorrect length selector vector, exiting...';
errmg3 = 'pardiff1 requires (16,14) size dat matrix, exiting...';
if nargin < 5
    trace = 0;
end;

% -----------------------------------
% Manual setting of tmax needed here.
% SET ME PLEASE!
% -----------------------------------
NS = 300;
ncond = 10;
N = 168; % number of observations. - actually 168, keep here, but adjust below
tmax = 4.0; % corrected.
np = 33;       
lp = length(Pvar) + length(Pfix);
if lp ~= np
     [name, errmg1], length(Pvar), length(Pfix), return;
end;
if length(Sel) ~= np
     [name, errmg2], length(Sel), return;
end;    
   
% Assemble parameter vector.
P = zeros(1,np);
P(Sel==1) = Pvar;
P(Sel==0) = Pfix;

Ptemp17 = P;
save Ptemp17 Ptemp17
if trace
   P
end;

q = 1.0; % From most highly parameterized model
if ~isempty(trace)
    if sum(size(trace)) > 2
       runid = trace(2);
       Savename = ['Ptemp17', int2str(runid), '.mat'];
       if sum(size(trace)) > 3
            q = trace(3);
       end
       trace = trace(1,1);
    else 
       Savename = 'Ptemp17.mat';
    end
    if trace
       P
    end
end
save(Savename, 'Ptemp17')


if any(size(Dmatrix) ~= [10,14])
     [name, errmg3], size(Dmatrix), return;
else
     Pred = zeros(10, 14);  % Predicted values.
     Pf1 = [0,.1,.3,.5,.7,.9,1.0];
     Od = ones(10,12);
     Wt = ones(1,10);
end;

% Fixing parameters equal across constrast
% If the first is 1, rest are 0, all set to the first.

% Constraints to fix V2a as constant
if Sel(7) == 0
   P(7) = P(6);
end
if Sel(8) == 0
   P(8) = P(7);
end
if Sel(9) == 0
   P(9) = P(8);
end
if Sel(10) == 0
   P(10) = P(9);
end


% Constraints to fix V1a to V2p
if Sel(16) == 0
   P(16) = P(1);
end
if Sel(17) == 0
   P(17) = P(2);
end
if Sel(18) == 0
   P(18) = P(3);
end
if Sel(19) == 0
   P(19) = P(4);
end
if Sel(20) == 0
   P(20) = P(5);
end



% Constraints to fix V1p as constant
if Sel(12) == 0
   P(12) = P(11);
end
if Sel(13) == 0
   P(13) = P(12);
end
if Sel(14) == 0
   P(14) = P(13);
end
if Sel(15) == 0
   P(15) = P(14);
end

% Constraints to fix etas constant within target present/target absent conditions
% Present
if Sel(23) == 0
   P(23) = P(22);
end
if Sel(24) == 0
   P(24) = P(23);
end
if Sel(25) == 0
   P(25) = P(24);
end
if Sel(26) == 0
   P(26) = P(25);
end
%Absent
if Sel(28) == 0
   P(28) = P(27);
end
if Sel(29) == 0
   P(29) = P(28);
end
if Sel(30) == 0
   P(30) = P(29);
end
if Sel(31) == 0
   P(31) = P(30);
end



V = P(1:20);
a = P(21);
eta = P(22:31);
Ter = P(32);
st = P(33);
sigma = 1.0; % Fixed


V2p = [V(1),V(1), V(2),V(2), V(3),V(3), V(4),V(4), V(5),V(5)];
V2a = [V(6),V(6), V(7),V(7), V(8),V(8), V(9),V(9), V(10),V(10)]; 
V1p = [V(11),V(11), V(12),V(12), V(13),V(13), V(14),V(14), V(15),V(15)];
V1a = [V(16),V(16), V(17),V(17), V(18),V(18), V(19),V(19), V(20),V(20)];
Teri = ones(2,1) * Ter;
Teri = Teri(:)';
%Ix = [1,3,5,7,9];
%U = [V2p(Ix);V2a(Ix);V1p(Ix);V1a(Ix)]

Etap = ones(2,1) * eta(1:5);  % First five are target present.
Etaa = ones(2,1) * eta(6:10); % Second five are target absent.
Etap = Etap(:)';
Etaa = Etaa(:)';
% Constrain Ter to be smaller than the smallest first data quantile.
epsn = .001;
Q1 = min(Dmatrix(:,[3,10]));
q1min = min(Q1);

% -----------------------------------------------------------------------------
%    v1...v20,         a      eta1           Ter,      st
%     1...20          21      22:31          32        29
% -----------------------------------------------------------------------------

Ub= [ 7.0*ones(1,20), 5.0, 5.0*ones(1,10),  q1min,   0.4]; 
Lb= [-7.0*ones(1,20), 0.4,    zeros(1,10),  0.02,     0 ];
Pub=[ 6.5*ones(1,20), 4.5, 4.4*ones(1,10),  q1min,  0.35];
Plb=[-6.5*ones(1,20), 0.5,    zeros(1,10),  0.02,   0.01 ];

if any(P - Ub > 0) | any(Lb - P > 0)
   G2 = 1e5 + ...
         1e3 * (sum(max(P - Pub, 0).^2) + sum(max(Plb - P).^2));
   X2 = G2;
   BIC = 0;
   if trace
        max(P - Ub, 0)
        max(Lb - P, 0)
   end;      
else
   % ----------------------
   % Also SET ME
   % ----------------------
   h = tmax/NS;
   Qni = [3:7,10:14];    % Index the quantiles in data matrix.
   Snid = [1,2,8,9];     % Index of summary stats in data matrix for double
   Snis = [8,9,1,2];     % Exchange for single
 
   G2 = 0;
   X2 = 0;
   for i = 1:10
        teri = Ter; 
        %i
        if rem(i, 2)  % 1, 3, 5, 7, 9 
           %disp('Double')
            Pi = [-V2a(i), -V2a(i),  V2p(i), V2p(i), Etaa(i), Etaa(i), Etap(i), Etap(i), ...
                  sigma, a, teri, st];          
        else % 2, 6, 6, 8, 10
            %disp('Single')
            Pi = [-V1a(i),-V1a(i),   -V1a(i),  V1p(i), Etaa(i), Etaa(i), Etaa(i), Etap(i), ...
                  sigma, a, teri, st];
           
        end
        %Pi
        [t,g1,g2,ig1,ig2,StPr,QnPr] = ved4sphere(Pi, tmax, 1);
        %plot(t, g1, t, g2)
        %StPr
        %pause
        %StPr
        
        % Normalize response probabilities for robustness
        StPr(1) = StPr(1)/(StPr(1) + StPr(2));  
        StPr(2) = StPr(2)/(StPr(1) + StPr(2));
        Pred(i,Qni) = [QnPr(:,1)',QnPr(:,2)'];  % transpose quantiles.
        %<ved4sphere> returns double, single RESPONSES, agrees with convention in "Duncan"
        Pred(i,Snid) = StPr([1,3,2,4]);
        Qnrlo = floor((Dmatrix(i,Qni) - teri)/ h) + 1;
        Qnrhi = ceil((Dmatrix(i,Qni) - teri) / h) + 1;
        Tlo = [t(Qnrlo(1:5)),t(Qnrlo(6:10))];
        Tfrac = (Dmatrix(i,Qni) - Tlo) ./ h;
        % Predicted values weighted implicitly by probability.
        Iglo = [ig1(Qnrlo(1:5)),    ig2(Qnrlo(6:10))];
        Ighi = [ig1(Qnrlo(1:5)+1),  ig2(Qnrlo(6:10)+1)];
        Ig = [Iglo(1:5) + Tfrac(1:5) .* (Ighi(1:5) - Iglo(1:5)), ...
        Iglo(6:10) + Tfrac(6:10) .* (Ighi(6:10) - Iglo(6:10))]; %(10x1)
        % Objective function.
        % Modified version to make a proper chi-square.
        PredProb = [diff([0, Ig(1:5), ig1(length(t))]),  ...
                diff([0, Ig(6:10), ig2(length(t))])];   % (12 x 1)
        PredProb = max(PredProb, eps);
        ObsProb = [diff(Dmatrix(i,1) * Pf1), diff(Dmatrix(i,8) * Pf1)];
        Ginc = sum((ObsProb .* log(ObsProb ./(PredProb+eps)))); % - ... 
                 %(ObsProb - PredProb))./Od(i,:));
        %[i, N * Ginc] % ###
        Resi = ObsProb .* log(ObsProb ./(PredProb+eps));
        Res(i,:) = N * Resi;        
        X2inc = N * sum((ObsProb - PredProb).^2./(PredProb+eps));
        % SSE: X2inc = N * sum((ObsProb - PredProb).^2); % ## SSE ##
        X2= X2 + X2inc;
        G2 = G2 + Wt(i) * Ginc;
     end;
     G2 = 2 * N * G2; % Because it is conventional.   
      Penalty =  1e5 * (sum(sum(max(P - Pub, 0).^2))  ...
                 + sum(max(max(Plb - P, 0).^2)));
     if trace == 2
         Penalty
         max(P - Pub, 0)
         max(Plb - P, 0)
         P
    end;
  % Calculate AIC on unpenalised G2
    nparams = sum(Sel);
    dfr = 11 * ncond - nparams;
    trueG2 = G2; %  Need /2 for group when we used 168.
    %N = N / 2;  % Because real 168 rather than 336
    adjG2 = trueG2 / q;
    k = nparams;
    QAIC = adjG2  + 2 * k; 
    BIC =  adjG2 + k * log(ncond * N);
    B = [trueG2, adjG2, BIC, QAIC, dfr]; 
    G2  = G2 + Penalty;
    X2  = X2 + Penalty;

    if ~isreal(G2)
        G2 = 1e8;
    end;
    if ~isreal(X2)
        X2 = 1e8;
    end;
    %fprintf('ET = %10.3f \n',etime(clock, t0));
end


