% Local update
clc; clear all; close all;
%--------------------------------------------------------------------------
% Load initial MGP
load('MoE20000_cop.mat','subdata','gprMdl','Mu','MuX','Priors','Sigma','SigmaX');
D = size(MuX,1);
K = size(MuX,2);
%--------------------------------------------------------------------------
% Maximum addtional points
max_point = 50; 
new_gprMdl = cell(1,K,max_point); % Store MGP
new_gprMdl{1} = gprMdl;
%------------------------------------------------------------------------
%% Random parameters
% Normal random samples for P1
% P1 = 60(kN)
mu_P1 = 60; COV_P1 =0.2; sigma_P1 = COV_P1*mu_P1; 
%------------------------------------------------------------------------
% Normal random samples for P2
% P2 = 40(kN)
mu_P2 = 40; COV_P2 =0.2; sigma_P2 = COV_P2*mu_P2; 
%------------------------------------------------------------------------
% Normal random samples for P3
% P3 = 10(kN)
mu_P3 = 10; COV_P3 =0.2; sigma_P3 = COV_P3*mu_P3; 
%------------------------------------------------------------------------
% Normal random samples for Es
% Es = 20000; % Elastic modulus (kN/cm2))
mu_Es = 20000; COV_Es =0.1; sigma_Es = COV_Es*mu_Es; 
%------------------------------------------------------------------------
% Normal random samples for L
% L = 100; % cm
mu_L = 100; COV_L =0.05; sigma_L = COV_L*mu_L; 
%------------------------------------------------------------------------
%% Main loop
% Lower and upper bounds of design variables
lb = ones(1,10);
ub = 20*ones(1,10);
initialPoint = 0.5*(lb+ub);
lamda = zeros(max_point,1); % Lamda vector
tole= zeros(max_point,1);
% Store output data
store_x = zeros(max_point,10);
store_x(1,:)= initialPoint;
store_fval = zeros(max_point,1);
store_fval(1,:) = objectives_real(initialPoint);
store_constraint = zeros(max_point,1);
%------------------------------------------------------------------------
% SDO
for i = 2:max_point
    formatSpec = 'Current iteration is %1.0f\n';
    fprintf(formatSpec,i)
  % fmincon
    options = optimoptions('fmincon','ConstraintTolerance',10^-6,'MaxIterations',50,'OptimalityTolerance',10^-6,'StepTolerance',10^-6,'PlotFcn',{'optimplotfval','optimplotx'},'Display','iter');
    [x,fval] = fmincon(@objectives_real,initialPoint,[],[],[],[],lb,ub,@(x)constraints_real(x, mu_P1, mu_P2, mu_P3, mu_Es, mu_L, new_gprMdl{i-1}, Priors, MuX, SigmaX,lamda(i-1)),options);
    store_x(i,:) = x; % store design variables
    store_fval(i,:) = fval; % store objective values
   % Mean vector
    X_mean = [mu_P1 mu_P2 mu_P3 mu_Es mu_L x];
    % Call current constraint function
    constr = MOGPE(X_mean, new_gprMdl{i-1}, Priors, MuX, SigmaX); 
    % Store constraint function values
    store_constraint(i) = constr; 
    % Inverse SA
    inver_y = inverse_y(x,new_gprMdl{i-1}, Priors, MuX, SigmaX,0.7*6.21*10^-3);
    % Update lamda
    lamda(i) = constr-inver_y;
    % stopping criterion
    tol = norm( lamda(i)- lamda(i-1));
    tole(i) = tol;
    if tol < 10^-9
        break;
    else
    % Call FEM
    g_3 = LSF3(mu_P1,mu_P2,mu_P3,mu_Es,mu_L,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10));
    % New point
    new_point = [mu_P1 mu_P2 mu_P3 mu_Es mu_L x g_3];
    % Comput weight vector
    w = ME_weight(new_point,Priors,Mu,Sigma);
    % Clusters that new_point belongs to
    [mx,cl]=max(w); 
    subdata{cl}=[subdata{cl};new_point]; % Add new_point to the cluster
    % Initial value of theta
    theta0 = 10*ones(1,D);
    % Upper and lower bounds of theta
    lob = 10^-3*ones(1,D);
    upb = 20*ones(1,D);
    % Local learning using DACE
    cluster = subdata{cl};
    gprMdl{cl} = dacefit(cluster(:,1:15), cluster(:,16), @regpoly2, @corrgauss, theta0, lob, upb);
    new_gprMdl{i} = gprMdl;
    gprMdl  = new_gprMdl{i};
    end
end
%------------------------------------------------------------------------
%% Additional functions
%------------------------------------------------------------------------
% Objective function
function y = objectives_real(x)
y = x(1) + x(2) + x(3) + x(4) + x(5) + x(6) + x(7) + x(8) + x(9) + x(10);
end
%------------------------------------------------------------------------
% Constraint function
function [c, ceq] = constraints_real(x,mu_P1, mu_P2, mu_P3, mu_Es, mu_L, gprMdl, Priors, MuX, SigmaX,delta)
% Mean vector
X_mean = [mu_P1 mu_P2 mu_P3 mu_Es mu_L x(1) x(2) x(3) x(4) x(5) x(6) x(7) x(8) x(9) x(10)];
g_3 = LSF3(mu_P1,mu_P2,mu_P3,mu_Es,mu_L,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10));
new_p = [X_mean g_3];
% Determine constraint value using MGP
g_mean = -MOGPE(X_mean, gprMdl, Priors, MuX, SigmaX)+delta;

% Inequality constraints
c = [g_mean];
% No equality constraints
ceq = [];
end
%------------------------------------------------------------------------
% Inverse SA
function [inv] = inverse_y(x,gprMdl, Priors, MuX, SigmaX,Pf_expect)
% Random parameters
% Normal random samples for P1
mu_P1 = 60; COV_P1 =0.2; sigma_P1 = COV_P1*mu_P1; 
%------------------------------------------------------------------------
% Normal random samples for P2
mu_P2 = 40; COV_P2 =0.2; sigma_P2 = COV_P2*mu_P2; 
%------------------------------------------------------------------------
% Normal random samples for P3
mu_P3 = 10; COV_P3 =0.2; sigma_P3 = COV_P3*mu_P3; 
%------------------------------------------------------------------------
% Normal random samples for Es
mu_Es = 20000; COV_Es =0.1; sigma_Es = COV_Es*mu_Es; 
%------------------------------------------------------------------------
% Normal random samples for L
mu_L = 100; COV_L =0.05; sigma_L = COV_L*mu_L; 
%------------------------------------------------------------------------
% Mean and standard variation vector
X_mean = [mu_P1 mu_P2 mu_P3 mu_Es mu_L x(1) x(2) x(3) x(4) x(5) x(6) x(7) x(8) x(9) x(10)];
Sigma_X = [sigma_P1 sigma_P2 sigma_P3 sigma_Es sigma_L zeros(1,10)]';
% First eight moments of the input variables
D = length(X_mean);
mux2 = zeros(D,1);
mux3 = zeros(D,1);
mux4 = zeros(D,1);
mux5 = zeros(D,1);
mux6 = zeros(D,1);
mux7 = zeros(D,1);
mux8 = zeros(D,1);
for i = 1:D
    mux2(i,:) = normal_cem(X_mean(i),Sigma_X(i),2);
    mux3(i,:) = normal_cem(X_mean(i),Sigma_X(i),3);
    mux4(i,:) = normal_cem(X_mean(i),Sigma_X(i),4);
    mux5(i,:) = normal_cem(X_mean(i),Sigma_X(i),5);
    mux6(i,:) = normal_cem(X_mean(i),Sigma_X(i),6);
    mux7(i,:) = normal_cem(X_mean(i),Sigma_X(i),7);
    mux8(i,:) = normal_cem(X_mean(i),Sigma_X(i),8);
end
% Gradient, hessianm, and mean value of the limit state function
func = @(u)Limit(u, gprMdl, Priors, MuX, SigmaX);
hess = hessian(func,X_mean);
hess((6:15),(6:15))=0;
[y_mean,grad] = MOGPE(X_mean, gprMdl, Priors, MuX, SigmaX);
grad((6:15),:)=0;
% y_mean
% mux2
% dia = diag(hess)
% 0.5*dot(mux2,dia)
% grad
% First four moments of the limit state function
[muy1,muy2,muy3,~] = four_ordery(y_mean,grad,hess,mux2,mux3,mux4,mux5,mux6,mux7,mux8);

% Inverse probability of failure
y2 = [-0.6:0.001:0.5];
Pf = zeros(length(y2),1);
for i = 1:length(y2)
    Pf(i)=third_SAB2(muy1,muy2,muy3,y2(i)); % Call SA
    if Pf(i)>=Pf_expect
        inv = y2(i);
        break;
    end
end
end
%------------------------------------------------------------------------
function [y] = Limit(u, gprMdl, Priors, MuX, SigmaX)
X= [u(1) u(2) u(3) u(4) u(5) u(6) u(7) u(8) u(9) u(10) u(11) u(12) u(13) u(14) u(15)];
y = MOGPE(X, gprMdl, Priors, MuX, SigmaX);
end
%------------------------------------------------------------------------

