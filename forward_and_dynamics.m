clear all
close all

syms g m1 m2 G1 G2 assume positive;
syms l1 l2 lc1 lc2 assume positive;
syms q1 q2 q dq1 dq2 dq ddq1 ddq2 assume real;
syms x_p y_p;
syms x1 y1 x2 y2;
syms I1 Ixx1 Iyy1 Izz1 Ixy1 Iyz1 Ixz1;
syms I2 Ixx2 Iyy2 Izz2 Ixy2 Iyz2 Ixz2;

% L(1)=Link('d', 0, 'a', 0.25, 'alpha', 0);
% L(2)=Link('d', 0, 'a', 0.25, 'alpha', 0);


%%
%%%%%% kinematics

%%%forward kinematics
% lc1=0.198;
% lc2=0.109;

DH_table=...
    [0,      l1,     0,     q1
     0,      l2,     0,     q2
    ];
DHc_table=...
    [0,      lc1,     0,     q1
     0,      lc2,     0,     q2
    ];

%assign the parameters
alpha=DH_table(:,1); 
a=DH_table(:,2);
d=DH_table(:,3);
q=DH_table(:,4);
 
A01 =  [cos(q(1))  -sin(q(1))    0   a(1)*cos(q(1));
        sin(q(1))   cos(q(1))     0   a(1)*sin(q(1));
            0            0        1         0       ;
            0            0        0         1      ];
        
A12 = [cos(q(2))   -sin(q(2))     0   a(2)*cos(q(2));
       sin(q(2))    cos(q(2))     0   a(2)*sin(q(2));
           0            0         1         0       ;
           0            0         0         1      ];
%%%% center DH
alpha=DHc_table(:,1); 
a_c=DHc_table(:,2);
d_c=DHc_table(:,3);
q=DHc_table(:,4);

A0c1 = [cos(q(1))  -sin(q(1))     0   a_c(1)*cos(q(1));
        sin(q(1))   cos(q(1))     0   a_c(1)*sin(q(1));
            0            0        1         0       ;
            0            0        0         1      ];
A1c2 = [cos(q(2))   -sin(q(2))    0   a_c(2)*cos(q(2));
       sin(q(2))    cos(q(2))     0   a_c(2)*sin(q(2));
           0            0         1         0       ;
           0            0         0         1      ];
       
%transfer matrix
T02 = simplify(A01*A12);
Tc02=simplify(A01*A1c2);
%%
%%% Inverse kinematics
% P = [x_p, y_p, 0];
%  D = (x_p.^2+y_p.^2-lc1^2-lc2^2)/(2*lc1*lc2);
%  q2= -abs(acos(D));
%  q1= abs(atan(y_p./x_p)-atan(sin(q2)*lc2./(lc1+cos(q2)*lc2)));

%%
%%% Jacobian matrix
o0=[0;0;0];
oc1=A0c1(1:3,4);
oc2=Tc02(1:3,4);
z=[0;0;1];
j1v=cross(z,(oc1-o0));
j2v=cross(z,(oc1-oc1));
j1=vertcat(j1v,z);
j2=vertcat(j2v,z);
J1=[j1 j2];

o1= A01(1:3,4);
j1_v=cross(z,(oc2-o0));
j2_v=cross(z,(oc2-o1));
%j1w = j2w = z;
j_1=vertcat(j1_v,z);
j_2=vertcat(j2_v,z);
J2=[j_1 j_2];

%J2_endeffector
o2=T02(1:3,4);
j_1_v=cross(z,(o2-o0));
j_2_v=cross(z,(o2-o1));
j__1=vertcat(j_1_v,z);
j__2=vertcat(j_2_v,z);
J2ee=[j__1 j__2]


%%
%%%%%%%% Dynamics
Jvc1 = J1(1:3,1:2);
Jvc2 = J2(1:3,1:2);

%%%% G %%%%
G1 = (m1*lc1+m2*l1)*g*cos(q1)+m2*lc2*g*cos(q1+q2);
G2 = m2*lc2*g*cos(q1+q2);
G=[G1;G2];

%%%% M %%%%
K=m1*(Jvc1)'*Jvc1+m2*(Jvc2)'*Jvc2;
I=[I1+I2,I2;I2,I2];
M=simplify(K+I);

%%%% C %%%%
q=[q1;q2];
dq=[dq1;dq2];
ddq=[ddq1;ddq2];

C=C_analyse(M,q1,q2,dq2,dq1);

%%%% tau %%%%
tau=simplify(M*ddq+C*dq+G);

%%%% substitute the value
% I1=[Ixx1,Ixy1,Ixz1;Ixy1,Iyy1,Iyz1;Ixz1,Iyz1,Izz1];
% I2=[Ixx2,Ixy2,Ixz2;Ixy2,Iyy2,Iyz2;Ixz2,Iyz2,Izz2];
% m1=0.374; 
% m2=0.129;
% g=9.81;
% lc1=0.198;
% lc2=0.109;

%%%%substitute the values
K_s=subs(K,{l1,lc1,lc2,m1,m2},{0.25,0.198,0.109,0.374,0.129});
%vpa(simplify(K_s),3)
M_s=simplify(K_s+subs(I, {I1, I2},{0.0027, 0.00113}));
%define precision
M_matrix=vpa(simplify(K_s));
G_matrix=vpa(subs(G,{g,l1,lc1,lc2,m1,m2},{9.81,0.25,0.198,0.109,0.374,0.129}));
C_matrix=vpa(subs(C,{l1,lc1,lc2,m1,m2},{0.25,0.198,0.109,0.374,0.129}));

tau_matrix = vpa(M_matrix*ddq+G_matrix+C_matrix*dq,5);

tau1=vpa(simplify(tau_matrix(1,1)),5);
tau2=tau_matrix(2,1);

