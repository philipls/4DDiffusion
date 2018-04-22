function [G2, B, Pred] = d4s(Pvar, Pfix, Sel, Dmatrix, trace)
% ================================================================
% Fit 4D diffusion model to Corbett & Smith (2017) single-target task
% With orthant bias.
% -45: separate target-absent drifts for 0t and 1t conditions.
% -46: 10 x etas
% -47: etas follow the stimulus, not the condition.
% Data are organized by probability of target-present response.
% Contrast x (Single, None)
% July 9, 2017
% Usage:
%     [G2, B, Pred] = d4s(Pvar,Pfix, Sel, Dmatrix, trace);
%
%      P = [v1...v15, a,  eta1...eta10,    Ter, st, c]
%          [1...15,   16, 17:21, 22:26,    27   28 29]
%          1:5: 1tp, 6:10: 1ta; 11:15: 0ta
% Second element of trace is run-number for parallel fits.
% ================================================================
t0 = clock;
name = 'D4S: ';
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
N = 168; % 336; % number of observations. (check this is true)
ncond = 10;
tmax = 4.0;  % Used 2.0 for original fits.
np = 29;       
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

Ptemp47 = P;
save Ptemp47 Ptemp47
if trace
   P
end;
q = 1.00; % From most highly parameterized model
if ~isempty(trace)
    if sum(size(trace)) > 2
       runid = trace(2);
       Savename = ['Ptemp47', int2str(runid), '.mat'];
       if sum(size(trace)) > 3
            q = trace(3);
       end
       trace = trace(1,1);
    else 
       Savename = 'Ptemp47.mat';
    end
    if trace
       P
    end
end
save(Savename, 'Ptemp47')



Ptemp47 = P;
save Ptemp47 Ptemp47
if trace
   P
end;
if any(size(Dmatrix) ~= [10,14])
     [name, errmg3], size(Dmatrix), return;
else
     Pred = zeros(10, 14);  % Predicted values.
     Pf1 = [0,.1,.3,.5,.7,.9,1.0];
     Od = ones(10,12);
     Wt = ones(1,10);
end;


% Constraints to fix V1a as constant
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


%Constraints to fix etas constant within a category
% Present
if Sel(18) == 0
   P(18) = P(17);
end
if Sel(19) == 0
   P(19) = P(18);
end
if Sel(20) == 0
   P(20) = P(19);
end
if Sel(21) == 0
   P(21) = P(20);
end
%Absent
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


V = P(1:15);
a = P(16);
eta = P(17:26);
Ter = P(27);
st = P(28);
c = P(29);
sigma = 1.0; % Fixed

V1p = ones(2,1) * V(1:5);
V1a = ones(2,1) * V(6:10);
V0a = ones(2,1) * V(11:15);

V1p = V1p(:)';
V1a = V1a(:)';
V0a = V0a(:)';

Etap = ones(2,1) * eta(1:5);  % First five are target present.
Etaa = ones(2,1) * eta(6:10); % Second five are target absent.
Etap = Etap(:)';
Etaa = Etaa(:)';

Eta = eta' * ones(1,5);
Eta = Eta(:)';


% Constrain Ter to be smaller than the smallest first data quantile.
epsn = .001;
Q1 = min(Dmatrix(:,[3,10]));
q1min = min(Q1);

% --------------------------------------------------------------
%    v1...v15,       a,      eta1...eta2,     Ter,    st      c
%     1...15         16        17...26         27     28      29
% --------------------------------------------------------------

Ub= [ 7.0*ones(1,15), 4.5,    5.0*ones(1,10),  q1min,  0.5   100.0]; 
Lb= [-7.0*ones(1,15), 0.4,       zeros(1,10),  0.02,     0 -100.0];
Pub=[ 6.5*ones(1,15), 4.0,    4.4*ones(1,10),  q1min,  0.4   10.0];
Plb=[-6.5*ones(1,15), 0.5,       zeros(1,10),  0.15,   0.01  -10.0];  % changed st bounds
                                               % was .10
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
        teri = Ter; % only one at the moment.
        %i
        if rem(i, 2)  % 1, 3, 5, 7, 9 
          %disp('Single')
          Pi = [-V1a(i),-V1a(i), -V1a(i), V1p(i), ....
                Etaa(i), Etaa(i), Etaa(i), Etap(i), sigma, a, Ter, st, c];          
        else % 2, 6, 6, 8, 10
           %disp('None')
           Pi = [-V0a(i), -V0a(i), -V0a(i), -V0a(i), ...
                Etaa(i), Etaa(i),Etaa(i),Etaa(i), sigma, a, Ter, st, c];            
        end
       % Pi
        % g1 is present, g2 is absent.
        %[i,Pi]
        [t,g1,g2,ig1,ig2,StPr,QnPr] = vbs4sphere(Pi, tmax, 1);
        %StPr
        %plot(t, g1, t, g2, '--')
        %pause
        % Normalize response probabilities for robustness
        StPr(1) = StPr(1)/(StPr(1) + StPr(2));  
        StPr(2) = StPr(2)/(StPr(1) + StPr(2));
        Pred(i,Qni) = [QnPr(:,1)',QnPr(:,2)'];  % transpose quantiles.
        Pred(i,Snid) = StPr([1,3,2,4]);
        % Predicted distribution functions linearly interpolated at 
        % the quantiles.
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


