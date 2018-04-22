function duqpfs(Data, Pred, ymax)
% ==================================================================
% Plots qlpf's for Duncan single-target task. Present/absent x contrast
% Data and Pred are 10 x 14 matrix. (2 tasks x 5 contrasts)
%
% Usage:
%     duqpfs(Data, Pred, [ymax])
% ==================================================================

% ----------------------------
% Millisecond scaling
% ----------------------------
Ix = [3:7,10:14];
Data(:,Ix) = 1000.0*Data(:,Ix);
Pred(:,Ix) = 1000.0*Pred(:,Ix);

if size(Data) ~= [10,14]
   disp('Data must be a 10 x 14 matrix, exiting...');
   return;
end;

if size(Pred) ~= [10,14]
   disp('Pred must be a 10 x 14 matrix, exiting...');
   return;
end;
y1 = 0.75 * min(min(min(Data(:,[3,10]))),min(min(Pred(:,[3,10]))));
y2 = 1.1 * max(max(max(Data(:,[6,14]))),max(max(Pred(:,[6,14]))));
y1 = 0.05;

if nargin < 3
   y2 = 2200;
else
   y2 = ymax;
end;

axhandle1=setfig2;
sat = 0.6;
red =0.25;
for cond = 1:2
    axes(axhandle1(cond));
    hold on
    Cx = (cond - 1) + [1, 3, 5, 7, 9]; % alternating conditions, odd numbers are double.
    OProb1 = Data(Cx,1);
    OProb2 = Data(Cx,8);
    OQnt1 =  Data(Cx,3:7);

    OQnt2 = Data(Cx,10:14);
    PProb = [Pred(Cx,1);Pred(Cx,8)];    
    PQnt =  [Pred(Cx,3:7);Pred(Cx,10:14)];
    [Temp,Px]=sort(PProb);

    PProba = Pred(Cx,1);
    PProbb = Pred(Cx,8);    
    PQnta =  Pred(Cx,3:7);
    PQntb =  Pred(Cx,10:14);
    [Temp1,Px1]=sort(PProba);
    [Temp2,Px2]=sort(PProbb);
    first = PProba(Px1);
    second = PProbb(Px2);
    symbol = ['o', 's', 'd', 'v', '^'];
    for i = 1:5
        plot(OProb1, OQnt1(:,i), symbol(i), ...
             'MarkerSize', 6, 'MarkerEdgeColor', [0,sat,0], 'MarkerFaceColor', [0,sat,0]);
         plot(OProb2, OQnt2(:,i), symbol(i), ...
             'MarkerSize', 6, 'MarkerEdgeColor', [red,0,sat], 'MarkerFaceColor',[red,0,sat]);            
        plot(PProba(Px1), PQnta(Px1, i), 'k.-' , 'MarkerSize', 10);
        plot(PProbb(Px2), PQntb(Px2, i), 'k.-' , 'MarkerSize', 10);
    end;
    set(gca, 'XLim', [0,1.0]);
    set(gca, 'YLim', [y1,y2]); 

    hold off
end;    
label(axhandle1(1), .60, .90, 'Single')
label(axhandle1(2), .60, .90, 'None')

