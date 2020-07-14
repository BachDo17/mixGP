function [g_3] = Truss_FEM(P1,P2,P3,Es,L,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10)
%----- Topology matrix Edof -------------------------------------
NJ = 6;
NR = 4;
NE = 10;
NDOF_T = 2*NJ;
NDOF = NDOF_T - NR;
% Edof = [element_ID dof1_ID dof2_ID dof3_ID dof4_ID] 
 Edof=[1   1  2  3  4;
       2   3  4  5  6;
       3   7  8  9  10;
       4   9  10 11 12;
       5   3  4  9  10;
       6   5  6  11 12;
       7   3  4  7  8;
       8   1  2  9  10;
       9   3  4  11 12;
      10   5  6  9  10];
  
%----- Global stiffness matrix K and load vector F --------------------
 K=zeros(NDOF_T);
 % Self-weight loading
 f_s = zeros(NDOF_T,1);
 f_s(2) = -[0.5*0.01*7850*L/100*A1/(10000)+0.5*sqrt(2)*0.01*7850*L/100*A8/(10000)];
 f_s(4) = -[3*0.5*0.01*7850*L/100*A1/(10000)+2*0.5*sqrt(2)*0.01*7850*L/100*A8/(10000)];
 f_s(6) = -[2*0.5*0.01*7850*L/100*A1/(10000)+0.5*sqrt(2)*0.01*7850*L/100*A8/(10000)];
 f_s(8) = -[0.5*0.01*7850*L/100*A1/(10000)+0.5*sqrt(2)*0.01*7850*L/100*A8/(10000)];
 f_s(10) = -[3*0.5*0.01*7850*L/100*A1/(10000)+2*0.5*sqrt(2)*0.01*7850*L/100*A8/(10000)];
 f_s(12) = -[2*0.5*0.01*7850*L/100*A1/(10000)+0.5*sqrt(2)*0.01*7850*L/100*A8/(10000)];
 % External point loading
 F_E=zeros(NDOF_T,1);
 F_E(4)=-P1;
 F_E(5)=P3;
 F_E(6)=-P2;
 F = f_s+F_E;
%----- Element properties ---------------------------------------
Ep = zeros(NE,2);
Ep(1,:) = [Es A1]; Ep(2,:) = [Es A2]; 
Ep(3,:) = [Es A3]; Ep(4,:) = [Es A4]; 
Ep(5,:) = [Es A5]; Ep(6,:) = [Es A6];
Ep(7,:) = [Es A7]; Ep(8,:) = [Es A8];
Ep(9,:) = [Es A9]; Ep(10,:) = [Es A10];

%----- Global coordinates and topology --------------------------
% Node coordinates [x y]
 Coord=[0 0;
        L 0;
        2*L 0;
        0 L;
        L L;
        2*L L];
 % Corresponing DOFs
 Dof=[ 1  2;
       3  4;
       5  6;
       7  8;
       9 10;
      11 12];
 
%----- Element coordinates -------------------------------------- 

 [Ex,Ey]=coordxtr(Edof,Coord,Dof,2);

%----- Create element stiffness matrices Ke and assemble into K -
% All the element matrices are computed and assembled in the loop
 for i=1:10
    Ke=bar2e(Ex(i,:),Ey(i,:),Ep(i,:));
    K=assem(Edof(i,:),K,Ke);
 end
 
%----- Solve the system of equations ----------------------------
% Boundary condition: which DOFs are fixed
 bc= [1 0;2 0;7 0;8 0];   
 % Dis: Displacement solution including boundary values
 % Reac: Reaction force vector
 [Dis,Reac] = solveq(K,F,bc);
 g_3 = 0.4 - abs(Dis(6,1));
%----- Element forces -------------------------------------------
% Ed: Element displacement matrix
 Ed=extract(Edof,Dis);
 N = zeros(NE,1);
 for i=1:NE
    N(i,:)=bar2s(Ex(i,:),Ey(i,:),Ep(i,:),Ed(i,:));
 end
%---------------------------- end -------------------------------
end
