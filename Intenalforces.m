function [Max_PC1,Max_PC2,Max_PC3,Max_PC4,Max_PB5,Max_PB6,Max_PB7,AxialMaxtrix,BendingMatrix] = Intenalforces(DL,S1,S2,L1,L2,SL,WL,ES,...
         d1,tw1,d2,tw2,d3,tw3,d4,tw4,d5,tw5,d6,tw6,d7,tw7,Fy,Fu)
% Calculate axial forces and beanding moments at sections along beams
% (0,...,L) with factored loads applied
%--------------------------------------------------------------------------
% Inputs
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
% Max_PC1,Max_PC2,Max_PC3,Max_PC4: Maximum performance of column groups 1-4
% Max_PB5,Max_PB6,Max_PB7: Maximum performance of beam groups 5-7
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
Q11 =(1.2*DL + 1.6*S1 + 1.6*L1); Q12 = (1.2*DL + 1.6*S2 + 1.6*L2); Q13 = (1.2*DL + 1.6*S1 + 1.6*L1); 
Q21 = (1.2*DL + 1.6*S1 + 1.6*L1); Q22 = (1.2*DL + 1.6*S2 + 1.6*L2); Q23 = (1.2*DL + 1.6*S1 + 1.6*L1); 
Q31 = (1.2*DL + 1.6*S1 + 1.6*L1); Q32 = (1.2*DL + 1.6*S2 + 1.6*L2); Q33 = (1.2*DL + 1.6*S1 + 1.6*L1);
Q41 = (1.2*DL + 1.6*SL); Q42 = (1.2*DL + 1.6*SL); Q43 = (1.2*DL + 1.6*SL); % in kN/m
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
AxialMaxtrix = zeros(9,NE);
BendingMatrix = zeros(9,NE);
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
    AxialMaxtrix(:,im) = es(:,1); % in kN
    BendingMatrix(:,im) = es(:,3); % in kN.cm
    Qlocal = LOCALEF(k,T,disMatrix(:,im),QfMatrix(:,im));% Local member end force vector
    QlocalMatrix(:,im) = Qlocal;% Storing calculated member end force vector to end force matrix
    Qglobal = GLOBALEF(T,Qlocal);% Global member end force vector
    RE = REACTION(JB,JE,NDOF,NSC,Qglobal);% Reaction vector
    ReacMatrix(:,im)= RE;% Reaction matrix
end
%--------------------------------------------------------------------------
%% Axial and moment at 9 sections along columns
Axial_colunm = zeros(9,16);
Moment_colunm = zeros(9,16);
for i = 1:16
    if i <=4
        Axial_colunm(:,i)= AxialMaxtrix(:,i);
        Moment_colunm(:,i)= BendingMatrix(:,i);
    end
    if i>4 && i <=8
        Axial_colunm(:,i)= AxialMaxtrix(:,i+3);
        Moment_colunm(:,i)= BendingMatrix(:,i+3);
    end
    if i>8 && i <=12
        Axial_colunm(:,i)= AxialMaxtrix(:,i+6);
        Moment_colunm(:,i)= BendingMatrix(:,i+6);
    end
    if i>12 && i <=16
        Axial_colunm(:,i)= AxialMaxtrix(:,i+9);
        Moment_colunm(:,i)= BendingMatrix(:,i+9);
    end    
end
%--------------------------------------------------------------------------
%% Moment at 9 sections along beams
Moment_beam = zeros(9,12);
for i = 1:12
    if i <=3
        Moment_beam(:,i) = BendingMatrix(:,i+4);
    end
    if i >3 && i <=6
        Moment_beam(:,i) = BendingMatrix(:,i+8);
    end
    if i >6 && i<=9
        Moment_beam(:,i) = BendingMatrix(:,i+12);
    end
    if i >9 && i<=12
        Moment_beam(:,i) = BendingMatrix(:,i+16);
    end
end
%--------------------------------------------------------------------------
%% Calculate available combined axial and flexural strengths for columns 
Ga = zeros(16,1);
Gb = zeros(16,1);
for i = 1:4
Ga(i,:) = 1;
end
Ga(5,:) = 2*(Ix1/350)/(Ix5/600);
Ga(6,:) = 2*(Ix3/350)/(Ix5/600+Ix7/300);
Ga(7,:) = Ga(6,:);
Ga(8,:) = Ga(5,:);
Ga(9,:) = (Ix1/350+Ix2/350)/(Ix5/600);
Ga(10,:) = (Ix3/350+Ix4/350)/(Ix5/600+Ix7/300);
Ga(11,:) = Ga(10,:);
Ga(12,:) = Ga(9,:);
Ga(13,:) = 2*(Ix2/350)/(Ix5/600);
Ga(14,:) = 2*(Ix4/350)/(Ix5/600+Ix7/300);
Ga(15,:) = Ga(14,:);
Ga(16,:) = Ga(13,:);

for i = 1:12
Gb(i,:) = Ga(i+4,:);
end
Gb(13,:) =  (Ix2/350)/(Ix6/600);
Gb(14,:) =  (Ix4/350)/(Ix6/600+Ix7/300);
Gb(15,:) =  Gb(14,:);
Gb(16,:) =  Gb(13,:);

Ke_column = zeros(16,1);
for i = 1:16
    Ke_column(i,:) = Unbraced_K(Ga(i,:),Gb(i,:));
end
% Column compressive strength
Pc_c = zeros(16,1);
Pc_c(1,:) = CompStr_AISC(d1,tw1,bf1,tf1,ES,Fy,Ke_column(1,:),350,0,0);  % Column group 1
Pc_c(2,:) = CompStr_AISC(d3,tw3,bf3,tf3,ES,Fy,Ke_column(2,:),350,0,0);  % Column group 3
Pc_c(3,:) = CompStr_AISC(d3,tw3,bf3,tf3,ES,Fy,Ke_column(3,:),350,0,0);  % Column group 3
Pc_c(4,:) = CompStr_AISC(d1,tw1,bf1,tf1,ES,Fy,Ke_column(4,:),350,0,0);  % Column group 1
Pc_c(5,:) = CompStr_AISC(d1,tw1,bf1,tf1,ES,Fy,Ke_column(5,:),350,0,0);  % Column group 1
Pc_c(6,:) = CompStr_AISC(d3,tw3,bf3,tf3,ES,Fy,Ke_column(6,:),350,0,0);  % Column group 3
Pc_c(7,:) = CompStr_AISC(d3,tw3,bf3,tf3,ES,Fy,Ke_column(7,:),350,0,0);  % Column group 3
Pc_c(8,:) = CompStr_AISC(d1,tw1,bf1,tf1,ES,Fy,Ke_column(8,:),350,0,0);  % Column group 1
Pc_c(9,:) = CompStr_AISC(d2,tw2,bf2,tf2,ES,Fy,Ke_column(9,:),350,0,0);  % Column group 2
Pc_c(10,:) = CompStr_AISC(d4,tw4,bf4,tf4,ES,Fy,Ke_column(10,:),350,0,0);% Column group 4
Pc_c(11,:) = CompStr_AISC(d4,tw4,bf4,tf4,ES,Fy,Ke_column(11,:),350,0,0);% Column group 4
Pc_c(12,:) = CompStr_AISC(d2,tw2,bf2,tf2,ES,Fy,Ke_column(12,:),350,0,0);% Column group 2
Pc_c(13,:) = CompStr_AISC(d2,tw2,bf2,tf2,ES,Fy,Ke_column(13,:),350,0,0);% Column group 2
Pc_c(14,:) = CompStr_AISC(d4,tw4,bf4,tf4,ES,Fy,Ke_column(14,:),350,0,0);% Column group 4
Pc_c(15,:) = CompStr_AISC(d4,tw4,bf4,tf4,ES,Fy,Ke_column(15,:),350,0,0);% Column group 4
Pc_c(16,:) = CompStr_AISC(d2,tw2,bf2,tf2,ES,Fy,Ke_column(16,:),350,0,0);% Column group 2

% Column tensile strength
Pc_t = zeros(16,1);
Pc_t(1,:) = TensStr_AISC(d1,tw1,bf1,tf1,Fy,Fu);% Column group 1
Pc_t(4,:) = Pc_t(1,:); Pc_t(5,:) = Pc_t(1,:); Pc_t(8,:) = Pc_t(1,:); % Column group 1

Pc_t(2,:) = TensStr_AISC(d3,tw3,bf3,tf3,Fy,Fu);% Column group 3
Pc_t(3,:) = Pc_t(2,:);Pc_t(6,:) = Pc_t(2,:);Pc_t(7,:) = Pc_t(2,:);% Column group 3

Pc_t(9,:) = TensStr_AISC(d2,tw2,bf2,tf2,Fy,Fu);% Column group 2
Pc_t(12,:) = Pc_t(9,:);Pc_t(13,:) = Pc_t(9,:);Pc_t(16,:) = Pc_t(9,:); % Column group 2

Pc_t(10,:) = TensStr_AISC(d4,tw4,bf4,tf4,Fy,Fu);% Column group 4
Pc_t(11,:) = Pc_t(10,:);Pc_t(14,:) = Pc_t(10,:);Pc_t(15,:) = Pc_t(10,:);% Column group 4

% Column flexural strength
Mcx_c = zeros(16,1);
Mcx_c(1,:) = FlexStrXc_AISC(d1,tw1,bf1,tf1,ES,Fy,350,1);
Mcx_c(4,:) = Mcx_c(1,:); Mcx_c(5,:) = Mcx_c(1,:); Mcx_c(8,:) = Mcx_c(1,:);

Mcx_c(2,:) = FlexStrXc_AISC(d3,tw3,bf3,tf3,ES,Fy,350,1);
Mcx_c(3,:) = Mcx_c(2,:);Mcx_c(6,:) = Mcx_c(2,:);Mcx_c(7,:) = Mcx_c(2,:);

Mcx_c(9,:) = FlexStrXc_AISC(d2,tw2,bf2,tf2,ES,Fy,350,1);
Mcx_c(12,:) = Mcx_c(9,:);Mcx_c(13,:) = Mcx_c(9,:);Mcx_c(16,:) = Mcx_c(9,:);

Mcx_c(10,:) = FlexStrXc_AISC(d4,tw4,bf4,tf4,ES,Fy,350,1);
Mcx_c(11,:) = Mcx_c(10,:);Mcx_c(14,:) = Mcx_c(10,:);Mcx_c(15,:) = Mcx_c(10,:);
% Mcy_c
% = FlexStrY_AISC(d,tw,bf,tf,E,Fy)
%--------------------------------------------------------------------------
%% Calculate available flexural strength for beams
Mcx_b = zeros(12,1);
Mcx_b(1,:) = FlexStrXb_AISC(d5,tw5,bf5,tf5,ES,Fy,600/2,1.3);
Mcx_b(3,:) = Mcx_b(1,:);
Mcx_b(4,:) = Mcx_b(1,:);Mcx_b(6,:) = Mcx_b(1,:);
Mcx_b(7,:) = Mcx_b(1,:); Mcx_b(9,:) = Mcx_b(1,:);

Mcx_b(2,:) = FlexStrXb_AISC(d7,tw7,bf7,tf7,ES,Fy,300,1.2);
Mcx_b(5,:) = Mcx_b(2,:);Mcx_b(8,:) = Mcx_b(2,:);Mcx_b(11,:) = Mcx_b(2,:);

Mcx_b(10,:) = FlexStrXb_AISC(d6,tw6,bf6,tf6,ES,Fy,600/2,1.3);
Mcx_b(12,:) = Mcx_b(10,:);
%--------------------------------------------------------------------------
%% Column performace
Axial_ratio = zeros(9,16);
for i = 1:16
    Axial_ratio(:,i) = Axial_colunm(:,i)./Pc_c(i,:);
end
PerformanceColumn= zeros(9,16);
Max_performaceColumn = zeros(16,1);
for i = 1:16
if max(abs(Axial_ratio(:,i)))>= 0.2
    PerformanceColumn(:,i) = abs(Axial_colunm(:,i)./Pc_c(i,:))+8*abs(Moment_colunm(:,i)./Mcx_c(i,:))/9;
else
    PerformanceColumn(:,i) = abs(Axial_colunm(:,i)./(2*Pc_c(i,:)))+abs(Moment_colunm(:,i)./Mcx_c(i,:));
end
Max_performaceColumn(i,:) = max(PerformanceColumn(:,i));
end
% (6) MGP -  Column group 1
GrC1 = [Max_performaceColumn(1),Max_performaceColumn(4),Max_performaceColumn(5),Max_performaceColumn(8)];
Max_PC1 = 1-max(GrC1); 

% (7) MGP -  Column group 2
GrC2 = [Max_performaceColumn(9),Max_performaceColumn(12),Max_performaceColumn(13),Max_performaceColumn(16)];
Max_PC2 = 1-max(GrC2); 

% (8) MGP -  Column group 3
GrC3 = [Max_performaceColumn(2),Max_performaceColumn(3),Max_performaceColumn(6),Max_performaceColumn(7)];
Max_PC3 = 1-max(GrC3);

% (9) MGP -  Column group 4
GrC4 = [Max_performaceColumn(10),Max_performaceColumn(11),Max_performaceColumn(14),Max_performaceColumn(15)];
Max_PC4 = 1-max(GrC4);
%--------------------------------------------------------------------------
%% Beam performace

PerformanceBeam = zeros(9,12);
Max_performaceBeam = zeros(12,1);
for i = 1:12
    PerformanceBeam(:,i) = Moment_beam(:,i)./Mcx_b(i,:);
    Max_performaceBeam(i,:) = max(abs(PerformanceBeam(:,i)));
end
% (10) MGP - Beam group 5
GrB5 = [Max_performaceBeam(1),Max_performaceBeam(3),Max_performaceBeam(4),Max_performaceBeam(6),Max_performaceBeam(7),Max_performaceBeam(9)];
Max_PB5 = 1-max(GrB5);

% (11) MGP - Beam group 6
GrB6 = [Max_performaceBeam(10),Max_performaceBeam(12)];
Max_PB6 = 1-max(GrB6);

% (12) MGP - Beam group 7
GrB7 = [Max_performaceBeam(2),Max_performaceBeam(5),Max_performaceBeam(8),Max_performaceBeam(11)];
Max_PB7 = 1-max(GrB7);
end

