%% Main optimization algorithm
clc; clear all
%% Load steel data
load('Steeldata.mat');
%% Load trained model
load('MGP20000_01.mat','subdata1','gprMdl1','MuX1','Priors1','SigmaX1','Mu1','Sigma1','subdata2','gprMdl2','MuX2','Priors2','SigmaX2','Mu2','Sigma2');
MuX1 = MuX1(1:22,:);
MuX2 = MuX2(1:24,:);
SigmaX1 = SigmaX1(1:22,1:22,:);
SigmaX2 = SigmaX2(1:24,1:24,:);
D1 = size(MuX1,1);
K1 = size(MuX1,2);
D2 = size(MuX2,1);
K2 = size(MuX2,2);
%--------------------------------------------------------------------------
%% Random parameters
% Normal random samples for DL
% DL = 20(kN/m)
no_DL = 20; mu_DL = no_DL; COV_DL = 0.1; sigma_DL = COV_DL*mu_DL;
%--------------------------------------------------------------------------
% Lognormal random samples for S1
% S1= 10(kN/m)
no_S1 = 10; mu_S1 = no_S1; COV_S1 =0.3; sigma_S1 = COV_S1*mu_S1;
%--------------------------------------------------------------------------
% Lognormal random samples for S2
% S2 = 5(kN/m)
no_S2 = 5; mu_S2 = no_S2; COV_S2 =0.3; sigma_S2 = COV_S2*mu_S2;
%--------------------------------------------------------------------------
% Lognormal random samples for L1
% L1 = 10(kN/m)
no_L1 = 10; mu_L1 = no_L1; COV_L1 =0.3; sigma_L1 = COV_L1*mu_L1;
%--------------------------------------------------------------------------
% Lognormal random samples for L2
% L2 = 5(kN/m)
no_L2 = 5; mu_L2 = no_L2; COV_L2 =0.3; sigma_L2 = COV_L2*mu_L2;
%--------------------------------------------------------------------------
% Lognormal random samples for SL
% SL = 5(kN/m)
no_SL = 5; mu_SL = no_SL; COV_SL =0.3; sigma_SL = COV_SL*mu_SL;
%--------------------------------------------------------------------------
% Lognormal random samples for WL
% WL = 8(kN)
no_WL = 8; mu_WL = no_WL; COV_WL =0.3; sigma_WL = COV_WL*mu_WL;
%--------------------------------------------------------------------------
% Normal random samples for ES
% ES = 21000(kN/cm2)
no_ES = 21000; mu_ES = no_ES; COV_ES =0.04; sigma_ES = COV_ES*mu_ES; 
%--------------------------------------------------------------------------
% Normal random samples for Fy
% Fy = 23.5(kN/cm2)
no_Fy = 23.5; mu_Fy = 1.10*no_Fy; COV_Fy =0.06; sigma_Fy = COV_Fy*mu_Fy; 
%--------------------------------------------------------------------------
% Normal random samples for Fu
% Fu = 36(kN/cm2)
no_Fu = 36; mu_Fu = 1.07*no_Fu; COV_Fu =0.08; sigma_Fu = COV_Fu*mu_Fu; 
%--------------------------------------------------------------------------
%% Main loop
max_point = 50;

new_gprMdl1 = cell(1,K1,max_point);
new_gprMdl1{1} = gprMdl1;

new_gprMdl2 = cell(1,K2,max_point);
new_gprMdl2{1} = gprMdl2;

% Lower and upper bounds of design variables
lb = ones(1,7);
ub = 12*ones(1,7);

initialPoint = randi(12,1,7); % Random
lamda = zeros(12,max_point);

tolerance= zeros(max_point,1);

store_x = zeros(max_point,7);
store_x(1,:)= initialPoint;

store_fval = zeros(max_point,1);
store_fval(1,:) = objective(initialPoint,columndata,beamdata);

store_inv = zeros(12,max_point);
store_constraint = zeros(12,max_point);

for i = 2:max_point
    % Print the current state
    formatSpec = 'Current iteration is %1.0f\n';
    fprintf(formatSpec,i)
   % ga
    % rng(1,'twister') % for reproducibility
   nvars = 7;
   bound = [lb; ub];
   populationSize = 2000;
   stallGenLimit = 1000;
   generations = 30;
   intcon = [1 2 3 4 5 6 7];
   options = optimoptions('ga','ConstraintTolerance',1e-12,'InitialPopulationMatrix',initialPoint,...
        'InitialPopulationRange',bound,'MaxGenerations',generations,...
        'PopulationSize',populationSize,'MaxStallGenerations',stallGenLimit,...
        'CrossoverFraction',0.65,...
        'EliteCount',1,...
        'FunctionTolerance',10^-6,...
        'Display','iter',... 
        'PlotFcn', @gaplotbestf);
    [x,fval] = ga(@(x)objective(x,columndata,beamdata),nvars,[],[],[],[],lb,ub,@(x)constraints_real(x,columndata,beamdata,mu_DL,mu_S1,mu_S2,mu_L1,mu_L2,mu_SL,mu_WL,mu_ES,mu_Fy,mu_Fu,new_gprMdl1{i-1},Priors1,MuX1,SigmaX1,new_gprMdl2{i-1},Priors2,MuX2,SigmaX2,lamda(:,i-1)),intcon,options);
    initialPoint = x;
%     x = round(x)
    store_x(i,:) = x;
%     fval
    store_fval(i,:) = fval;
    % Mapping
    p = MapVariables(x,columndata,beamdata);
    % Depth and web thickness of column and beam sections
    d1 = p(1); tw1 = p(2); d2 = p(4); tw2 = p(5);
    d3 = p(7); tw3 = p(8); d4 = p(10); tw4 = p(11);
    d5 = p(13); tw5 = p(14); d6 = p(16); tw6 = p(17); d7 = p(19); tw7 = p(20);
    % Mean vectors
    X_mean1 = [mu_DL mu_S1 mu_S2 mu_L1 mu_L2 mu_SL mu_WL mu_ES ...
         d1 tw1 d2 tw2 d3 tw3 d4 tw4 d5 tw5 d6 tw6 d7 tw7];
    X_mean2 = [mu_DL mu_S1 mu_S2 mu_L1 mu_L2 mu_SL mu_WL mu_ES ...
         d1 tw1 d2 tw2 d3 tw3 d4 tw4 d5 tw5 d6 tw6 d7 tw7 mu_Fy mu_Fu];
    % LSFs at mean vectors
    store_constraint(1:5,i) = MOGPE(X_mean1, new_gprMdl1{i-1}, Priors1, MuX1, SigmaX1,5);
    store_constraint(6:12,i)  = MOGPE(X_mean2, new_gprMdl2{i-1}, Priors2, MuX2, SigmaX2,7);
    
    inver_y = inverse_g(x,columndata,beamdata,new_gprMdl1{i-1},Priors1,MuX1,SigmaX1,6.68*10^-2,new_gprMdl2{i-1},Priors2,MuX2,SigmaX2,1.35*10^-3);
    store_inv(:,i) = inver_y;
    lamda(:,i) = store_constraint(:,i)-inver_y;
    
    tol = norm(lamda(:,i)- lamda(:,i-1));
    tolerance(i) = tol;
    % Check convergence
    if tol < 10^-3
        break;
    else
    if store_x(i,:)==store_x(i-1,:)
        break;
    else
    % Update MGP
    [maxDrift,peakDrift,MaxNormBeamDis5,MaxNormBeamDis6,MaxNormBeamDis7] = Deformation(mu_DL,mu_S1,mu_S2,mu_L1,mu_L2,mu_SL,mu_WL,mu_ES,...
    d1,tw1,d2,tw2,d3,tw3,d4,tw4,d5,tw5,d6,tw6,d7,tw7);

    [Max_PC1,Max_PC2,Max_PC3,Max_PC4,Max_PB5,Max_PB6,Max_PB7,~,~] = Intenalforces(mu_DL,mu_S1,mu_S2,mu_L1,mu_L2,mu_SL,mu_WL,mu_ES,...
         d1,tw1,d2,tw2,d3,tw3,d4,tw4,d5,tw5,d6,tw6,d7,tw7,mu_Fy,mu_Fu);
     
    new_point1 = [mu_DL mu_S1 mu_S2 mu_L1 mu_L2 mu_SL mu_WL mu_ES d1 tw1 d2 tw2 d3 tw3 d4 tw4 d5 tw5 d6 tw6 d7 tw7...
                 maxDrift peakDrift MaxNormBeamDis5 MaxNormBeamDis6 MaxNormBeamDis7];
    
    new_point2 = [mu_DL mu_S1 mu_S2 mu_L1 mu_L2 mu_SL mu_WL mu_ES d1 tw1 d2 tw2 d3 tw3 d4 tw4 d5 tw5 d6 tw6 d7 tw7 mu_Fy mu_Fu ...
                 Max_PC1 Max_PC2 Max_PC3 Max_PC4 Max_PB5 Max_PB6 Max_PB7];
             
    w1 = ME_weight(X_mean1,Priors1,MuX1,SigmaX1);
    w2 = ME_weight(X_mean2,Priors2,MuX2,SigmaX2) ; 
    [~,cl1]=max(w1); % Clusters that new_point1 belongs to
%     cl1
    [~,cl2]=max(w2); % Clusters that new_point2 belongs to
%     cl2
    % Model 01
    subdata1{cl1}=[subdata1{cl1};new_point1]; % Add new_point1 to the cluster
    % Initial value of theta
    theta01 = 10*ones(1,D1);
    % Upper and lower bounds of theta
    lob1 = 10^-3*ones(1,D1);
    upb1 = 20*ones(1,D1);
    % Local learning using DACE
    cluster1 = subdata1{cl1};
    gprMdl1{cl1} = dacefit(cluster1(:,1:22), cluster1(:,23:27), @regpoly2, @corrgauss, theta01, lob1, upb1);
    new_gprMdl1{i} = gprMdl1;
    gprMdl1  = new_gprMdl1{i};
    % Model 02
    subdata2{cl2}=[subdata2{cl2};new_point2]; % Add new_point1 to the cluster
    % Initial value of theta
    theta02 = 10*ones(1,D2);
    % Upper and lower bounds of theta
    lob2 = 10^-3*ones(1,D2);
    upb2 = 20*ones(1,D2);
    % Local learning using DACE
    cluster2 = subdata2{cl2};
    gprMdl2{cl2} = dacefit(cluster2(:,1:24), cluster2(:,25:31), @regpoly2, @corrgauss, theta02, lob2, upb2);
    new_gprMdl2{i} = gprMdl2;
    gprMdl2  = new_gprMdl2{i};
    end
    end
end
%--------------------------------------------------------------------------
%% Objective function
function f = objective(x,columndata,beamdata)
x = round(x);
p = MapVariables(x,columndata,beamdata);
% Total weight of the structure
f = 4*3.5*p(3) + 4*3.5*p(6) + 4*3.5*p(9) + 4*3.5*p(12) +...
    6*6*p(15) + 2*6*p(18) + 4*3*p(21);
end
%--------------------------------------------------------------------------
%% Constraint functions
function [c, ceq] = constraints_real(x,columndata,beamdata,mu_DL,mu_S1,mu_S2,mu_L1,mu_L2,mu_SL,mu_WL,mu_ES,mu_Fy,mu_Fu,gprMdl1,Priors1,MuX1,SigmaX1,gprMdl2,Priors2,MuX2,SigmaX2,lamda)
%   [c, ceq] = calculates the constraints
% Problem parameters
p = MapVariables(x,columndata,beamdata);

% Depth and web thickness of column and beam sections
d1 = p(1); tw1 = p(2); d2 = p(4); tw2 = p(5);
d3 = p(7); tw3 = p(8); d4 = p(10); tw4 = p(11);
d5 = p(13); tw5 = p(14); d6 = p(16); tw6 = p(17); d7 = p(19); tw7 = p(20);

% Flange thickness and width of column sections
tf1 = (2.091*tw1*10-3.3595)/10;
tf2 = (2.091*tw2*10-3.3595)/10;
tf3 = (2.091*tw3*10-3.3595)/10;
tf4 = (2.091*tw4*10-3.3595)/10;
bf1 = d1;
bf2 = d2;
bf3 = d3;
bf4 = d4;

% Flange thickness and width of beam sections
tf5 = (1.6522*tw5*10-0.8304)/10;
tf6 = (1.6522*tw6*10-0.8304)/10;
tf7 = (1.6522*tw7*10-0.8304)/10;
bf_b = [46;55;64;73;82;91;100;110;120;135;150;160;170;180;190;200;210;220];
d_b = [80;100;120;140;160;180;200;220;240;270;300;330;360;400;450;500;550;600];
pp = polyfit(d_b,bf_b,2);
bf5 = polyval(pp,d5*10)/10;
bf6 = polyval(pp,d6*10)/10;
bf7 = polyval(pp,d7*10)/10;

% Mean vectors
X_mean1 = [mu_DL mu_S1 mu_S2 mu_L1 mu_L2 mu_SL mu_WL mu_ES ...
         d1 tw1 d2 tw2 d3 tw3 d4 tw4 d5 tw5 d6 tw6 d7 tw7];
X_mean2 = [mu_DL mu_S1 mu_S2 mu_L1 mu_L2 mu_SL mu_WL mu_ES ...
         d1 tw1 d2 tw2 d3 tw3 d4 tw4 d5 tw5 d6 tw6 d7 tw7 mu_Fy mu_Fu];
     
% Responses at mean vector
% y_mean1:5 by 1; grad1: 22 by 5
y_mean1 = MOGPE(X_mean1, gprMdl1, Priors1, MuX1, SigmaX1,5);

% y_mean2:7 by 1; grad1: 24 by 7
y_mean2 = MOGPE(X_mean2, gprMdl2, Priors2, MuX2, SigmaX2,7);

% Inter-story drift constraint
itedrift = y_mean1(1);
% Total building drift constraint
totaldrift = y_mean1(2);

% Beam deflection constraint
beamdeflection5 = y_mean1(3);
beamdeflection6 = y_mean1(4);
beamdeflection7 = y_mean1(5);

deform = [-itedrift + lamda(1);...
          -totaldrift + lamda(2);...
          -beamdeflection5 + lamda(3);...
          -beamdeflection6 + lamda(4);...
          -beamdeflection7 + lamda(5)];
      
% Constraints on the maximum stress in columns and beams
column_stress1 = y_mean2(1);
column_stress2 = y_mean2(2);
column_stress3 = y_mean2(3);
column_stress4 = y_mean2(4);
beam_stress5 = y_mean2(5);
beam_stress6 = y_mean2(6);
beam_stress7 = y_mean2(7);

stress = [-column_stress1 + lamda(6);...
          -column_stress2 + lamda(7);...
          -column_stress3 + lamda(8);...
          -column_stress4 + lamda(9);...
          -beam_stress5 + lamda(10);...
          -beam_stress6 + lamda(11);...
          -beam_stress7 + lamda(12)];

% Inequality constraints on geometry at beam-column

geo_cb1 = bf5 - bf1; % flange column 1 > flange beam(5)
geo_cb2 = bf5 - bf3; % flange column 3 > flange beam(5)
geo_cb3 = bf7 - bf3; % flange column 3 > flange beam(7)
geo_cb4 = bf5 - bf2; % flange column 2 > flange beam(5)
geo_cb5 = bf5 - bf4; % flange column 4 > flange beam(5)
geo_cb6 = bf7 - bf4; % flange column 4 > flange beam(7)
geo_cb7 = bf6 - bf2; % flange column 2 > flange beam(6)
geo_cb8 = bf6 - bf4; % flange column 4 > flange beam(6)
geo_bb1 = d7 - d5; % height beam 5 > height beam(7)
geo_bb2 = d7 - d6; % height beam 6 > height beam(7)

%  Inequality constraints on depths of column sections of two consecutive stories
geo_cc1 = d2 - d1;      % height column 1 > height column 2
geo_cc2 = d1 - d2 - 16; % height column 1 - height column 2 < 160mm
geo_cc3 = d4 - d3;      % height column 3 > height column 4 
geo_cc4 = d3 - d4 -16;  % height column 3 - height column 4 < 160
geo_cc5 = d1 - d3;      % height column 1 < height column 3
geo_cc6 = d2 - d4;      % height column 2 < height column 4

geometry = [geo_cb1;geo_cb2;geo_cb3;geo_cb4;geo_cb5;geo_cb6;geo_cb7;geo_cb8;...
            geo_bb1;geo_bb2;...
            geo_cc1;geo_cc2;geo_cc3;geo_cc4;geo_cc5;geo_cc6];
        
% All inequality constraints
c = [deform;stress;geometry];
% No equality constraints
ceq = [];
end
%--------------------------------------------------------------------------
%% Inverse probability function
function [inv] = inverse_g(x,columndata,beamdata,gprMdl1,Priors1,MuX1,SigmaX1,Pf_expect1,gprMdl2,Priors2,MuX2,SigmaX2,Pf_expect2)
% Initialization
p = MapVariables(x,columndata,beamdata);
% Depth and web thickness of column and beam sections
d1 = p(1); tw1 = p(2); d2 = p(4); tw2 = p(5);
d3 = p(7); tw3 = p(8); d4 = p(10); tw4 = p(11);
d5 = p(13); tw5 = p(14); d6 = p(16); tw6 = p(17); d7 = p(19); tw7 = p(20);
%--------------------------------------------------------------------------
% Normal random samples for DL
% DL = 20(kN/m)
no_DL = 20; mu_DL = no_DL; COV_DL = 0.1; sigma_DL = COV_DL*mu_DL;
%--------------------------------------------------------------------------
% Lognormal random samples for S1
% S1= 10(kN/m)
no_S1 = 10; mu_S1 = no_S1; COV_S1 =0.3; sigma_S1 = COV_S1*mu_S1;
%--------------------------------------------------------------------------
% Lognormal random samples for S2
% S2 = 5(kN/m)
no_S2 = 5; mu_S2 = no_S2; COV_S2 =0.3; sigma_S2 = COV_S2*mu_S2;
%--------------------------------------------------------------------------
% Lognormal random samples for L1
% L1 = 10(kN/m)
no_L1 = 10; mu_L1 = no_L1; COV_L1 =0.3; sigma_L1 = COV_L1*mu_L1;
%--------------------------------------------------------------------------
% Lognormal random samples for L2
% L2 = 5(kN/m)
no_L2 = 5; mu_L2 = no_L2; COV_L2 =0.3; sigma_L2 = COV_L2*mu_L2;
%--------------------------------------------------------------------------
% Lognormal random samples for SL
% SL = 5(kN/m)
no_SL = 5; mu_SL = no_SL; COV_SL =0.3; sigma_SL = COV_SL*mu_SL;
%--------------------------------------------------------------------------
% Lognormal random samples for WL
% WL = 8(kN)
no_WL = 8; mu_WL = no_WL; COV_WL =0.3; sigma_WL = COV_WL*mu_WL;
%--------------------------------------------------------------------------
% Normal random samples for ES
% ES = 21000(kN/cm2)
no_ES = 21000; mu_ES = no_ES; COV_ES =0.04; sigma_ES = COV_ES*mu_ES; 
%--------------------------------------------------------------------------
% Normal random samples for Fy
% Fy = 23.5(kN/cm2)
no_Fy = 23.5; mu_Fy = 1.10*no_Fy; COV_Fy =0.06; sigma_Fy = COV_Fy*mu_Fy; 
%--------------------------------------------------------------------------
% Normal random samples for Fu
% Fu = 36(kN/cm2)
no_Fu = 36; mu_Fu = 1.07*no_Fu; COV_Fu =0.08; sigma_Fu = COV_Fu*mu_Fu; 
%--------------------------------------------------------------------------
% Mean and variance vectors
X_mean1 = [mu_DL mu_S1 mu_S2 mu_L1 mu_L2 mu_SL mu_WL mu_ES ...
         d1 tw1 d2 tw2 d3 tw3 d4 tw4 d5 tw5 d6 tw6 d7 tw7];
X_mean2 = [mu_DL mu_S1 mu_S2 mu_L1 mu_L2 mu_SL mu_WL mu_ES ...
         d1 tw1 d2 tw2 d3 tw3 d4 tw4 d5 tw5 d6 tw6 d7 tw7 mu_Fy mu_Fu];
Sigma_X1 = [sigma_DL sigma_S1 sigma_S2 sigma_L1 sigma_L2 sigma_SL sigma_WL sigma_ES zeros(1,14)]';
Sigma_X2 = [sigma_DL sigma_S1 sigma_S2 sigma_L1 sigma_L2 sigma_SL sigma_WL sigma_ES zeros(1,14) sigma_Fy sigma_Fu]';
D1 = length(X_mean1);
D2 = length(X_mean2);
% First three moments of the input variables
mux1_2 = zeros(D1,1);
mux1_3 = zeros(D1,1);
for i = 1:D1
    mux1_2(i,:) = lognormal_cem(X_mean1(i),Sigma_X1(i),2);
    mux1_3(i,:) = lognormal_cem(X_mean1(i),Sigma_X1(i),3);
end
    mux1_2(1,:) = normal_cem(X_mean1(1),Sigma_X1(1),2);
    mux1_3(1,:) = normal_cem(X_mean1(1),Sigma_X1(1),3);
    
    mux1_2(8,:) = normal_cem(X_mean1(8),Sigma_X1(8),2);
    mux1_3(8,:) = normal_cem(X_mean1(8),Sigma_X1(8),3);    
    
mux2_2 = zeros(D2,1);
mux2_3 = zeros(D2,1);
for i = 1:D2
    mux2_2(i,:) = lognormal_cem(X_mean2(i),Sigma_X2(i),2);
    mux2_3(i,:) = lognormal_cem(X_mean2(i),Sigma_X2(i),3);
end
    mux2_2(1,:) = normal_cem(X_mean2(1),Sigma_X2(1),2);
    mux2_3(1,:) = normal_cem(X_mean2(1),Sigma_X2(1),3);
    
    mux2_2(8,:) = normal_cem(X_mean2(8),Sigma_X2(8),2);
    mux2_3(8,:) = normal_cem(X_mean2(8),Sigma_X2(8),3);
    
    mux2_2(23,:) = normal_cem(X_mean2(23),Sigma_X2(23),2);
    mux2_3(23,:) = normal_cem(X_mean2(23),Sigma_X2(23),3);
    
    mux2_2(24,:) = normal_cem(X_mean2(24),Sigma_X2(24),2);
    mux2_3(24,:) = normal_cem(X_mean2(24),Sigma_X2(24),3);
    
% Responses and their gradients at mean vector
% y_mean1:5 by 1; grad1: 22 by 5
[y_mean1,grad1] = MOGPE(X_mean1, gprMdl1, Priors1, MuX1, SigmaX1,5);
grad1(9:22,:)=0;

% y_mean2:7 by 1; grad1: 24 by 7
[y_mean2,grad2] = MOGPE(X_mean2, gprMdl2, Priors2, MuX2, SigmaX2,7);
grad2(9:22,:)=0;

% First three moments of the limit state functions
% (1) For maximum inter-story drift
[muy1_1,muy1_2,muy1_3] = three_momenty(y_mean1(1,:),grad1(:,1),mux1_2,mux1_3);
% (2) For maximum total drift
[muy2_1,muy2_2,muy2_3] = three_momenty(y_mean1(2,:),grad1(:,2),mux1_2,mux1_3);
% (3) For maximum vertical deflection of beams in group 5
[muy3_1,muy3_2,muy3_3] = three_momenty(y_mean1(3,:),grad1(:,3),mux1_2,mux1_3);
% (4) For maximum vertical deflection of beams in group 6
[muy4_1,muy4_2,muy4_3] = three_momenty(y_mean1(4,:),grad1(:,4),mux1_2,mux1_3);
% (5) For maximum vertical deflection of beams in group 7
[muy5_1,muy5_2,muy5_3] = three_momenty(y_mean1(5,:),grad1(:,5),mux1_2,mux1_3);

% (6) For maximum stress of columns in group 1
[muy6_1,muy6_2,muy6_3] = three_momenty(y_mean2(1,:),grad2(:,1),mux2_2,mux2_3);
% (7) For maximum stress of columns in group 2
[muy7_1,muy7_2,muy7_3] = three_momenty(y_mean2(2,:),grad2(:,2),mux2_2,mux2_3);
% (8) For maximum stress of columns in group 3
[muy8_1,muy8_2,muy8_3] = three_momenty(y_mean2(3,:),grad2(:,3),mux2_2,mux2_3);
% (9) For maximum stress of columns in group 4
[muy9_1,muy9_2,muy9_3] = three_momenty(y_mean2(4,:),grad2(:,4),mux2_2,mux2_3);
% (10) For maximum stress of beams in group 5
[muy10_1,muy10_2,muy10_3] = three_momenty(y_mean2(5,:),grad2(:,5),mux2_2,mux2_3);
% (11) For maximum stress of beams in group 6
[muy11_1,muy11_2,muy11_3] = three_momenty(y_mean2(6,:),grad2(:,6),mux2_2,mux2_3);
% (12) For maximum stress of beams in groups 7
[muy12_1,muy12_2,muy12_3] = three_momenty(y_mean2(7,:),grad2(:,7),mux2_2,mux2_3);

% Inverse probability of failure
inv = zeros(12,1);

y2 = [-1:0.001:1];
Pf_1 = zeros(length(y2),1);
Pf_2 = zeros(length(y2),1);
Pf_3 = zeros(length(y2),1);
Pf_4 = zeros(length(y2),1);
Pf_5 = zeros(length(y2),1);
Pf_6 = zeros(length(y2),1);
Pf_7 = zeros(length(y2),1);
Pf_8 = zeros(length(y2),1);
Pf_9 = zeros(length(y2),1);
Pf_10 = zeros(length(y2),1);
Pf_11 = zeros(length(y2),1);
Pf_12 = zeros(length(y2),1);

for i = 1:length(y2)
    Pf_1(i)=third_SAB2(muy1_1,muy1_2,muy1_3,y2(i));
    if Pf_1(i)>=Pf_expect1
        inv(1) = y2(i);
        break;
    end
end

for i = 1:length(y2)
    Pf_2(i)=third_SAB2(muy2_1,muy2_2,muy2_3,y2(i));
    if Pf_2(i)>=Pf_expect1
        inv(2) = y2(i);
        break;
    end
end

for i = 1:length(y2)
    Pf_3(i)=third_SAB2(muy3_1,muy3_2,muy3_3,y2(i));
    if Pf_3(i)>=Pf_expect1
        inv(3) = y2(i);
        break;
    end
end

for i = 1:length(y2)
    Pf_4(i)=third_SAB2(muy4_1,muy4_2,muy4_3,y2(i));
    if Pf_4(i)>=Pf_expect1
        inv(4) = y2(i);
        break;
    end
end

for i = 1:length(y2)
    Pf_5(i)=third_SAB2(muy5_1,muy5_2,muy5_3,y2(i));
    if Pf_5(i)>=Pf_expect1
        inv(5) = y2(i);
        break;
    end
end

for i = 1:length(y2)
    Pf_6(i)=third_SAB2(muy6_1,muy6_2,muy6_3,y2(i));
    if Pf_6(i)>=Pf_expect2
        inv(6) = y2(i);
        break;
    end
end

for i = 1:length(y2)
    Pf_7(i)=third_SAB2(muy7_1,muy7_2,muy7_3,y2(i));
    if Pf_7(i)>=Pf_expect2
        inv(7) = y2(i);
        break;
    end
end

for i = 1:length(y2)
    Pf_8(i)=third_SAB2(muy8_1,muy8_2,muy8_3,y2(i));
    if Pf_8(i)>=Pf_expect2
        inv(8) = y2(i);
        break;
    end
end

for i = 1:length(y2)
    Pf_9(i)=third_SAB2(muy9_1,muy9_2,muy9_3,y2(i));
    if Pf_9(i)>=Pf_expect2
        inv(9) = y2(i);
        break;
    end
end

for i = 1:length(y2)
    Pf_10(i)=third_SAB2(muy10_1,muy10_2,muy10_3,y2(i));
    if Pf_10(i)>=Pf_expect2
        inv(10) = y2(i);
        break;
    end
end

for i = 1:length(y2)
    Pf_11(i)=third_SAB2(muy11_1,muy11_2,muy11_3,y2(i));
    if Pf_11(i)>=Pf_expect2
        inv(11) = y2(i);
        break;
    end
end

for i = 1:length(y2)
    Pf_12(i)=third_SAB2(muy12_1,muy12_2,muy12_3,y2(i));
    if Pf_12(i)>=Pf_expect2
        inv(12) = y2(i);
        break;
    end
end
end
%--------------------------------------------------------------------------