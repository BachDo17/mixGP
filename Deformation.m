function [maxDrift,peakDrift,MaxNormBeamDis5,MaxNormBeamDis6,MaxNormBeamDis7] = Deformation(DL,S1,S2,L1,L2,SL,WL,ES,...
         d1,tw1,d2,tw2,d3,tw3,d4,tw4,d5,tw5,d6,tw6,d7,tw7)
%% clc, clear all, close all % clear screen
%--------------------------------------------------------------------------
% Model inputs
% DL = Dead load(kN/cm)
% S1 = Short term live load 1(kN/cm)
% S2 = Short term live load 2(kN/cm)
% L1 = Long term live load 1(kN/cm)
% L2 = Long term live load 2(kN/cm)
% SL = Snow load (kN/cm)
% WL = Wind load (kN)
% ES = Elastic modulus (kN/cm2)
% d,tw = Sectional depth and web thickness
%--------------------------------------------------------------------------
% Outputs
% maxDrift: Maximum inter-story drift
% peakDrift: Overall lateral displacement
% MaxNormBeamDis5,MaxNormBeamDis6,MaxNormBeamDis7: Maximum deflection of beam groups 5-7
%--------------------------------------------------------------------------
% Compute A = Crossectional area (cm2) and % and Ix = Moment inertia (cm4)for column and beam groups
% Groups 1-4: Columns tf = 2.091tw-3.3595 and bf/d = 1
tf1 = (2.091*tw1*10-3.3595)/10; 
tf2 = (2.091*tw2*10-3.3595)/10; 
tf3 = (2.091*tw3*10-3.3595)/10; 
tf4 = (2.091*tw4*10-3.3595)/10;
bf1 = d1; 
bf2 = d2; 
bf3 = d3; 
bf4 = d4;


A1 = 2*bf1*tf1+tw1*(d1-2*tf1);
A2 = 2*bf2*tf2+tw2*(d2-2*tf2);
A3 = 2*bf3*tf3+tw3*(d3-2*tf3);
A4 = 2*bf4*tf4+tw4*(d4-2*tf4);

Ix1 = (d1^3*bf1)/12-(bf1-tw1)*(d1-2*tf1)^3/12;
Ix2 = (d2^3*bf2)/12-(bf2-tw2)*(d2-2*tf2)^3/12;
Ix3 = (d3^3*bf3)/12-(bf3-tw3)*(d3-2*tf3)^3/12;
Ix4 = (d4^3*bf4)/12-(bf4-tw4)*(d4-2*tf4)^3/12;

% Group 5-7: Beams tf = 1.6522tw-0.8304 and bf = -0.0004d^2+0.6384d-6.3582
tf5 = (1.6522*tw5*10-0.8304)/10; 
tf6 = (1.6522*tw6*10-0.8304)/10; 
tf7 = (1.6522*tw7*10-0.8304)/10;

bf_b = [46;55;64;73;82;91;100;110;120;135;150;160;170;180;190;200;210;220];
d_b = [80;100;120;140;160;180;200;220;240;270;300;330;360;400;450;500;550;600];
pp = polyfit(d_b,bf_b,2);

bf5 = polyval(pp,d5*10)/10; 
bf6 = polyval(pp,d6*10)/10; 
bf7 = polyval(pp,d7*10)/10;

A5 = 2*bf5*tf5+tw5*(d5-2*tf5);
A6 = 2*bf6*tf6+tw6*(d6-2*tf6);
A7 = 2*bf7*tf7+tw7*(d7-2*tf7);

Ix5 = (d5^3*bf5)/12-(bf5-tw5)*(d5-2*tf5)^3/12;
Ix6 = (d6^3*bf6)/12-(bf6-tw6)*(d6-2*tf6)^3/12;
Ix7 = (d7^3*bf7)/12-(bf7-tw7)*(d7-2*tf7)^3/12;
%--------------------------------------------------------------------------
%% Uniformly distributed loads
Q11 =(DL + S1 + L1); Q12 = (DL + S2 + L2); Q13 = (DL + S1 + L1); 
Q21 = (DL + S1 + L1); Q22 = (DL + S2 + L2); Q23 = (DL + S1 + L1); 
Q31 = (DL + S1 + L1); Q32 = (DL + S2 + L2); Q33 = (DL + S1 + L1);
Q41 = (DL + SL); Q42 = (DL + SL); Q43 = (DL + SL); % in kN/m
% Distributed loads in local directions
Q_local = 0.01*[0 0;0 0;0 0;0 0;...
           0 -Q11;0 -Q12;0 -Q13;...
           0 0;0 0;0 0;0 0;...
           0 -Q21;0 -Q22;0 -Q23;...
           0 0;0 0;0 0;0 0;...
           0 -Q31;0 -Q32;0 -Q33;...
           0 0;0 0;0 0;0 0;...
            0 -Q41;0 -Q42;0 -Q43]; % in kN/cm
%--------------------------------------------------------------------------
%% PART I: Joint coordinate maxtrix, COORD(NJx2)
NJ = 20; % Number of joints
COORD = 100*[0 0;6 0;9 0;15 0;...
         0 3.5;6 3.5;9 3.5;15 3.5;...
         0 7;6 7;9 7;15 7;...
         0 10.5;6 10.5;9 10.5;15 10.5;...
         0 14;6 14;9 14;15 14]; % in cm
%--------------------------------------------------------------------------
%% PART II: Support data, MSUP(NSx4)
% 0, 0, 0 free
% 1, 0, 0 roller with horizontal reaction
% 0, 1, 0 roller with vertical reaction
% 1, 1, 0 hinge
% 0, 0, 1 support which prevents rotation, but not translation
% 1, 1, 1 fixed
NS = 4; % Number of joints that are attached to supports
MSUP = [1 1 1 1;2 1 1 1;3 1 1 1;4 1 1 1]; % Support data matrix, MSUP(NSx4)
%--------------------------------------------------------------------------
%% PART III: Material properties = Elastic modulus vector, EM(MNPx1)
MNP = 1;% Number of material property sets (E)
EM = [ES]; % Elastic modulus vector, EM(MNPx1)
%--------------------------------------------------------------------------
% PART IV:  Cross-sectional properties, CP(NCP,2), including area and
% moment inertia
NCP = 7;% Number of cross-sectional property sets
% Cross-sectional property vector, CP(NCPx2)[A Ix]
CP = [A1 Ix1;...% column, group 1
      A2 Ix2;...% column, group 2
      A3 Ix3;...% column, group 3
      A4 Ix4;...% column, group 4
      A5 Ix5;...% beam, group 5
      A6 Ix6;...% beam, group 6
      A7 Ix7];% beam, group 7
%--------------------------------------------------------------------------
%% PART V:   Storing member data including begining and end joints, material
% property No., and cross-section type no., MPRP(NEx4)
NE = 28; % Number of elements
MPRP = [1	5	1	1;...
        2	6	1	3;...
        3	7	1	3;...
        4	8	1	1;...
        5	6	1	5;...
        6	7	1	7;...
        7	8	1	5;...
        5	9	1	1;...
        6	10	1	3;...
        7	11	1	3;...
        8	12	1	1;...
        9	10	1	5;...
        10	11	1	7;...
        11	12	1	5;...
        9	13	1	2;...
        10	14	1	4;...
        11	15	1	4;...
        12	16	1	2;...
        13	14	1	5;...
        14	15	1	7;...
        15	16	1	5;...
        13	17	1	2;...
        14	18	1	4;...
        15	19	1	4;...
        16	20	1	2;...
        17	18	1	6;...
        18	19	1	7;...
        19	20	1	6]; % Integer member data matrix, MPRP(NEx4)
%--------------------------------------------------------------------------
%% PART VI:  Storing joint loads, PJ(NJLx3)
NJL = 4; % Number of joints that are subjected to external loads 
JP =[5;9;13;17]; % Joints where external loads applied, JP(NJLx1)
PJ = [WL 0 0;...
      WL 0 0;...
      WL 0 0;...
      WL 0 0];% Joint Load matrix, PJ(NJLx3)
%--------------------------------------------------------------------------
% Storing member load matrices, MP(NMLx2)and PM(NMLx4)
NML = 12; % Number of different load types applying to members
% Load types include:
% 1 = point load W,l1
% 2 = bending moment load M,l1
% 3 = uniformly distributed load w,l1,l2
% 4 = trapezoidal distributed load w1,w2,l1,l2
% 5 = axial point load, W,l1
% 6 = axially distributed load w,l1,l2
MP = [5 3;...
      6 3;...
      7 3;...
     12 3;...
     13 3;...
     14 3;...
     19 3;...
     20 3;...
     21 3;...
     26 3;...
     27 3;...
     28 3];% Member load matrix,MP(NMLx2)
% PM = [W,M,w,or w1|w2(if type 4),0(otherwise)|l1|l2(if load type 3,4,or 6)
% ,0(otherwise)]
PM = 0.01*[Q11 0 0 0;...
      Q12 0 0 0;...
      Q13 0 0 0;...
      Q21 0 0 0;...
      Q22 0 0 0;...
      Q23 0 0 0;...
      Q31 0 0 0;...
      Q32 0 0 0;...
      Q33 0 0 0;...
      Q41 0 0 0;...
      Q42 0 0 0;...
      Q43 0 0 0];% Member load matrix,PM(NMLx4) in kN/cm
%--------------------------------------------------------------------------
%% PART VII:Determining number of degrees of freedom NDOF of the structure
NR = COUNTRESTR(MSUP);% The number of restrained coordinates
NDOF = 3*NJ-NR; %  Number of degrees of freedom
%
% Part VIII:Forming the structure coordinate number vector, NSC(the number
% of structure coordinates per joint (NCJT) times the number of joints of
% the structure(NJ),1) 
NSC = NSCF(MSUP,NJ,NDOF); %  Sructure coordinate number vector, NSC(3*NJx1)
%--------------------------------------------------------------------------
%% PART IX:  Generating the structure stiffness matrix and the structure
% load vector due to member load
S = zeros(NDOF);% Structure stiffness matrix, S(NDOFxNDOF)
for im = 1:NE
    JB = MPRP(im,1);% Beginning joint number 
    JE = MPRP(im,2);% End joint number 
    I1 = MPRP(im,3);% Material no.
    I2 = MPRP(im,4);% Crossection type no.
    ex = [COORD(JB,1) COORD(JE,1)];
    ey = [COORD(JB,2) COORD(JE,2)];
    ep = [EM(I1,1) CP(I2,:)];
    T  = TRANS(ex,ey);% Transformation matrix
    k  = BEAM2E(ex,ey,ep);% Member stiffness matrix in local coordinate system
    K  = T'*k*T;% Member stiffness matrix in global coordinate system
    S  = STORES(JB,JE,NDOF,NSC,K,S);
end
% Storing structure load vector due to member load, Pf(NDOFx1)
Pf = zeros(NDOF,1);% Structure load vector due to member load, Pf(NDOFx1)
QfMatrix = zeros(6,NE);% Local member end force vector matrix
LMatrix = zeros(NE,1);
for im = 1:NE
    JB = MPRP(im,1);% Beginning joint number 
    JE = MPRP(im,2);% End joint number
    ex = [COORD(JB,1) COORD(JE,1)];
    ey = [COORD(JB,2) COORD(JE,2)];
    b  = [ex(2)-ex(1);ey(2)-ey(1)];
    L = sqrt(b'*b); % Member length 
    LMatrix(im,1)=L;
    T  = TRANS(ex,ey);% Transformation matrix
    Qf = zeros(6,1);% Local member end force vector,Qf(6x1)
    for iml = 1:NML
        if MP(iml,1)== im
           LTYPE = MP(iml,2);% Load type
           w1 = PM(iml,1);% Load type value 1
           w2 = PM(iml,2);% Load type value 2
           l1 = PM(iml,3);% Load loacation value 1
           l2 = PM(iml,4);% Load loacation value 2
           Qf = STOREQF(LTYPE,w1,w2,l1,l2,L,Qf);% Local member end force vector,Qf(6x1) 
        end
    end
    QfMatrix(:,im) = Qf;% Storing local member end force vector
    Ff = T'*Qf;% Global member end force vector,Ff(6x1)
    Pf = STOREPF(JB,JE,NDOF,NSC,Ff,Pf);% Structure load vector due to member load, Pf(NDOFx1)
end
%--------------------------------------------------------------------------
%% PART X:   Storing joint loads in the structure load vector
P = zeros(NDOF,1);
P = STOREP(JP,PJ,NDOF,NSC,P);% Joint load vector, P(NDOFx1)
Pt = P - Pf;% Joint load vector minus Structure load vector due to member load
%--------------------------------------------------------------------------
%% PART XI:  Calculatiing the structure's joint displacements by solving the
% stiffness relation, KU = P, using Gauss-Jordan elimination.
d = GaussPivot(S,Pt); % Solve for joint displacement vector, d(NDOFx1), by Gauss Pivot
%--------------------------------------------------------------------------
%% PART XII: Determining the member end force vectors and the support
% reaction vector
disMatrix = zeros(6,NE);% Storing displacements of joints of each member,disMatrix(6xNE) 
for im = 1:NE
    JB = MPRP(im,1);% Beginning joint number 
    JE = MPRP(im,2);% End joint number 
    V = DETERDISG(JB,JE,NDOF,NSC,d);% Displacement vector of each member
    disMatrix(:,im) = V;
end
disMatrix; % Global displacement matrix (6xNE)
%--------------------------------------------------------------------------
%% Postprocessing: determining local member end force vectors,Qlocal(6,1),the support
% Reaction vectors, ReacMatrix
QlocalMatrix = zeros(6,NE);% End force matrix
ReacMatrix = zeros(6,NE);% Reaction matrix
BeamDis = zeros(9,NE); % Beam displacement maxtrix
for im = 1:NE
    JB = MPRP(im,1);% Beginning joint number 
    JE = MPRP(im,2);% End joint number
    ex = [COORD(JB,1) COORD(JE,1)];
    ey = [COORD(JB,2) COORD(JE,2)];
    I1 = MPRP(im,3);% Material no.
    I2 = MPRP(im,4);% Crossection type no.
    ep = [EM(I1,1) CP(I2,:)];
    ed = disMatrix(:,im)';
    eq = Q_local(im,:);
    T  = TRANS(ex,ey);% Transformation matrix
    k  = BEAM2E(ex,ey,ep);% Member stiffness matrix in local coordinate system
    [es,edi,eci] = BEAM2S(ex,ey,ep,ed,eq,9);
    BeamDis(:,im) = edi(:,2); % in cm
    Qlocal = LOCALEF(k,T,disMatrix(:,im),QfMatrix(:,im));% Local member end force vector
    QlocalMatrix(:,im) = Qlocal;% Storing calculated member end force vector to end force matrix
    Qglobal = GLOBALEF(T,Qlocal);% Global member end force vector
    RE = REACTION(JB,JE,NDOF,NSC,Qglobal);% Reaction vector
    ReacMatrix(:,im)= RE;% Reaction matrix
end
%--------------------------------------------------------------------------
%% (1) MGP - Maximum interstory drift in cm
disMatrixDrift = JODISM(disMatrix,MPRP,NJ);
Drift = [disMatrixDrift(8,2)- disMatrixDrift(4,2);...
        disMatrixDrift(12,2)- disMatrixDrift(8,2);...
        disMatrixDrift(16,2)- disMatrixDrift(12,2);...
        disMatrixDrift(20,2)- disMatrixDrift(16,2)];
maxDrift = 1-max(abs(Drift))/(350/300);% constraint check: if <= 0 failure
%--------------------------------------------------------------------------
%% (2) MGP - Total drift
peakDrift = 1-disMatrixDrift(20,2)/(1400/400);% constraint check: if <= 0 failure
%--------------------------------------------------------------------------
%% Long-term beam deflections
% (3) MGP - Beam group 5
NormBeamDis5 = zeros(6,1);
NormBeamDis5(1,:) = 1.5*abs(min(BeamDis(:,5)))/(LMatrix(5,:)/360);
NormBeamDis5(2,:) = 1.5*abs(min(BeamDis(:,7)))/(LMatrix(7,:)/360);
NormBeamDis5(3,:) = 1.5*abs(min(BeamDis(:,12)))/(LMatrix(12,:)/360);
NormBeamDis5(4,:) = 1.5*abs(min(BeamDis(:,14)))/(LMatrix(14,:)/360);
NormBeamDis5(5,:) = 1.5*abs(min(BeamDis(:,19)))/(LMatrix(19,:)/360);
NormBeamDis5(6,:) = 1.5*abs(min(BeamDis(:,21)))/(LMatrix(21,:)/360);
MaxNormBeamDis5 = 1-max(NormBeamDis5);
%--------------------------------------------------------------------------
% (4) MGP - Beam group 6
NormBeamDis6 = zeros(2,1);
NormBeamDis6(1,:) = 1.5*abs(min(BeamDis(:,26)))/(LMatrix(26,:)/360);
NormBeamDis6(2,:) = 1.5*abs(min(BeamDis(:,28)))/(LMatrix(28,:)/360);
MaxNormBeamDis6 = 1-max(NormBeamDis6);
%--------------------------------------------------------------------------
% (5) MGP -  Beam group 7
NormBeamDis7 = zeros(4,1);
NormBeamDis7(1,:) = 1.5*abs(min(BeamDis(:,6)))/(LMatrix(6,:)/360);
NormBeamDis7(2,:) = 1.5*abs(min(BeamDis(:,13)))/(LMatrix(13,:)/360);
NormBeamDis7(3,:) = 1.5*abs(min(BeamDis(:,20)))/(LMatrix(20,:)/360);
NormBeamDis7(4,:) = 1.5*abs(min(BeamDis(:,27)))/(LMatrix(27,:)/360);
MaxNormBeamDis7 = 1-max(NormBeamDis7);
end
