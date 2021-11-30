clear all
close all
clc

% Trajectory for Point-to-Point Motion with Polynomials (p. 188 ff.)
% Matlab code: (c) Uwe Mettin, 2010
% Modified: Szabolcs Fodor, 2015

%% Joint Space

% query: initial/target position, velocities and accelerations of joints
q_0T = [   0, 2; 
        -0.5, 1];
v_0T = [   0, 0;
           0, 0];
a_0T = [   0, 0;
           0, 0];

% time
t_int = [0 3];
t = linspace(t_int(1),t_int(2),100);

% constraints matrix
A_poly = [1, t_int(1), t_int(1)^2, t_int(1)^3, t_int(1)^4, t_int(1)^5;...
          0, 1, 2*t_int(1), 3*t_int(1)^2, 4*t_int(1)^3, 5*t_int(1)^4; ...
          0, 0, 2, 6*t_int(1), 12*t_int(1)^2, 20*t_int(1)^3; ...
          1, t_int(2), t_int(2)^2, t_int(2)^3, t_int(2)^4, t_int(2)^5;...
          0, 1, 2*t_int(2), 3*t_int(2)^2, 4*t_int(2)^3, 5*t_int(2)^4; ...
          0, 0, 2, 6*t_int(2), 12*t_int(2)^2, 20*t_int(2)^3];       

% coefficients joint 1
poly_1 = A_poly \ [q_0T(1,1); v_0T(1,1); a_0T(1,1); q_0T(1,2); v_0T(1,2); a_0T(1,2);];
poly_1 = [poly_1(6), poly_1(5), poly_1(4), poly_1(3), poly_1(2), poly_1(1)];
      
% coefficients joint 2
poly_2 = A_poly \ [q_0T(2,1); v_0T(2,1); a_0T(2,1); q_0T(2,2); v_0T(2,2); a_0T(2,2)];
poly_2 = [poly_2(6), poly_2(5), poly_2(4), poly_2(3), poly_2(2), poly_2(1)];

q1 = polyval(poly_1,t);
q2 = polyval(poly_2,t);
Dq1 = polyval(polyder(poly_1),t);
Dq2 = polyval(polyder(poly_2),t);
Ddq1ref = polyval(polyder(polyder(poly_1)),t);
Ddq2ref = polyval(polyder(polyder(poly_2)),t);
l1 = 0.25; 
l2 = 0.25; 

%%
%%%%%% generate matrix with corresponding time stamp 
%%%%%% to use as input of 'from worksapce' block
sim_time = 3;    
sim_step = 100;    

for i=1:sim_step 
    time_stamp(i,1) = i * sim_time / sim_step;
end
%%% time stamp
theta1 = [time_stamp,q1'];
theta2 = [time_stamp,q2'];
Dtheta1 = [time_stamp,Dq1'];
Dtheta2 = [time_stamp,Dq2'];
Ddtheta1 = [time_stamp,Ddq1ref'];
Ddtheta2 = [time_stamp,Ddq2ref'];

