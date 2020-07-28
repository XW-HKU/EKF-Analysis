function [X_prid,X_next,RotM_next,P_prid,P_next,F] = System_Mode_all_inputs(X_in,RotM_last,Z_in,P_in,dt,mode)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%X_def=[v;q;omega;a_s;tau;ba;bg;bm;g];
%Z=[m_v;m_a;m_omega;m_m];
global COMP_T;
X2(:,1)=zeros(42,1);X2(9,1)=-9.8;
P2(:,:,1)=eye(length(X2(:,1)));

Q_p=0.0000*ones(1,3);
Q_v=0.002^2*0.1*ones(1,3);
Q_q=0.0000*ones(1,4);
Q_a_s=0.00000*ones(1,3);
Q_a_s1=0.000*ones(1,3);
Q_a_s2=0*ones(1,3);
Q_a_s3=0*ones(1,3);
Q_omega=0.0*ones(1,3);
Q_tau=0.00*ones(1,3);
Q_tau1=0*ones(1,3);
Q_tau2=0*ones(1,3);
Q_ba=0.002^2*0.01*ones(1,3);
Q_bg=0.002^2*0.01*ones(1,3);
Q_l_ic=0.000000*ones(1,3);
Q_ptheta=0.002^2*0.5*ones(1,3);

Q=diag([Q_p,Q_v,Q_ba,Q_bg,Q_l_ic,Q_ptheta]);
R=diag([0.0001*ones(1,3),0.002*ones(1,3)]);

% get inputs and assignments current states varibles
[ p1, p2, p3, v1, v2, v3, a_s1, a_s2, a_s3, a_s_11, a_s_12, a_s_13, a_s_21, a_s_22, a_s_23, a_s_31, a_s_32, a_s_33, omega1, omega2, omega3, tau1, tau2, tau3, tau_11, tau_12, tau_13, tau_21, tau_22, tau_23, ba1, ba2, ba3, bg1, bg2, bg3, l_ic1, l_ic2, l_ic3] = Assignment_MEKF_HighOrder(X_in);
[ R1_1, R1_2, R1_3, R2_1, R2_2, R2_3, R3_1, R3_2, R3_3 ] = RotMatAssi(RotM_last);
m_p1=Z_in(1);m_p2=Z_in(2);m_p3=Z_in(3);
a_s_in1=Z_in(4);a_s_in2=Z_in(5);a_s_in3=Z_in(6);
omega_in1=Z_in(7);omega_in2=Z_in(8);omega_in3=Z_in(9);
m_m1=Z_in(10);m_m2=Z_in(11);m_m3=Z_in(12);

% calculate states pridiction
tic
X_prid_p=[p1 + dt*v1
                                                                          p2 + dt*v2
                                                                          p3 + dt*v3
        v1 + dt*(R1_1*(a_s_in1 - ba1) + R1_2*(a_s_in2 - ba2) + R1_3*(a_s_in3 - ba3))
        v2 + dt*(R2_1*(a_s_in1 - ba1) + R2_2*(a_s_in2 - ba2) + R2_3*(a_s_in3 - ba3))
 v3 + dt*(R3_1*(a_s_in1 - ba1) + R3_2*(a_s_in2 - ba2) + R3_3*(a_s_in3 - ba3) + 49/5)
                                                                                 ba1
                                                                                 ba2
                                                                                 ba3
                                                                                 bg1
                                                                                 bg2
                                                                                 bg3
                                                                               l_ic1
                                                                               l_ic2
                                                                               l_ic3
                                                  0
                                                  0
                                                  0];
                                              
RotM_prid=RotM_last*expm(cross_mat(dt*(Z_in(7:9)-[bg1 bg2 bg3]')));

% calculate outputs pridiction
[ p1, p2, p3, v1, v2, v3, a_s1, a_s2, a_s3, a_s_11, a_s_12, a_s_13, a_s_21, a_s_22, a_s_23, a_s_31, a_s_32, a_s_33, omega1, omega2, omega3, tau1, tau2, tau3, tau_11, tau_12, tau_13, tau_21, tau_22, tau_23, ba1, ba2, ba3, bg1, bg2, bg3, l_ic1, l_ic2, l_ic3] = Assignment_MEKF_HighOrder(X_in);
[ R1_1, R1_2, R1_3, R2_1, R2_2, R2_3, R3_1, R3_2, R3_3 ] = RotMatAssi(RotM_prid);
Z_prid_p=[p1 + R1_1*l_ic1 + R1_2*l_ic2 + R1_3*l_ic3
 p2 + R2_1*l_ic1 + R2_2*l_ic2 + R2_3*l_ic3
 p3 + R3_1*l_ic1 + R3_2*l_ic2 + R3_3*l_ic3
                                conj(R1_1)
                                conj(R1_2)
                                conj(R1_3)];

% calculate F and H matrixes
F_p =[...
[ 1, 0, 0, dt,  0,  0,              0,              0,              0,   0,   0,   0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 1, 0,  0, dt,  0,              0,              0,              0,   0,   0,   0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 1,  0,  0, dt,              0,              0,              0,   0,   0,   0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  1,  0,  0, -dt*conj(R1_1), -dt*conj(R1_2), -dt*conj(R1_3),   0,   0,   0, 0, 0, 0, R1_3*a_s2*dt - R1_2*a_s3*dt, R1_1*a_s3*dt - R1_3*a_s1*dt, R1_2*a_s1*dt - R1_1*a_s2*dt]
[ 0, 0, 0,  0,  1,  0, -dt*conj(R2_1), -dt*conj(R2_2), -dt*conj(R2_3),   0,   0,   0, 0, 0, 0, R2_3*a_s2*dt - R2_2*a_s3*dt, R2_1*a_s3*dt - R2_3*a_s1*dt, R2_2*a_s1*dt - R2_1*a_s2*dt]
[ 0, 0, 0,  0,  0,  1, -dt*conj(R3_1), -dt*conj(R3_2), -dt*conj(R3_3),   0,   0,   0, 0, 0, 0, R3_3*a_s2*dt - R3_2*a_s3*dt, R3_1*a_s3*dt - R3_3*a_s1*dt, R3_2*a_s1*dt - R3_1*a_s2*dt]
[ 0, 0, 0,  0,  0,  0,              1,              0,              0,   0,   0,   0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,              0,              1,              0,   0,   0,   0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,              0,              0,              1,   0,   0,   0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,              0,              0,              0,   1,   0,   0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,              0,              0,              0,   0,   1,   0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,              0,              0,              0,   0,   0,   1, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,              0,              0,              0,   0,   0,   0, 1, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,              0,              0,              0,   0,   0,   0, 0, 1, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,              0,              0,              0,   0,   0,   0, 0, 0, 1,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,              0,              0,              0, -dt,   0,   0, 0, 0, 0,                           1,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,              0,              0,              0,   0, -dt,   0, 0, 0, 0,                           0,                           1,                           0]
[ 0, 0, 0,  0,  0,  0,              0,              0,              0,   0,   0, -dt, 0, 0, 0,                           0,                           0,                           1]
];
F_p(end-2:end,end-2:end)=expm(cross_mat(-dt*[omega1 omega2 omega3]));

H =[[ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, conj(R1_1), conj(R1_2), conj(R1_3), R1_3*l_ic2 - R1_2*l_ic3, R1_1*l_ic3 - R1_3*l_ic1, R1_2*l_ic1 - R1_1*l_ic2]
[ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, conj(R2_1), conj(R2_2), conj(R2_3), R2_3*l_ic2 - R2_2*l_ic3, R2_1*l_ic3 - R2_3*l_ic1, R2_2*l_ic1 - R2_1*l_ic2]
[ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, conj(R3_1), conj(R3_2), conj(R3_3), R3_3*l_ic2 - R3_2*l_ic3, R3_1*l_ic3 - R3_3*l_ic1, R3_2*l_ic1 - R3_1*l_ic2]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,          0,          0,          0,                       0,             -conj(R1_3),              conj(R1_2)]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,          0,          0,          0,              conj(R1_3),                       0,             -conj(R1_1)]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,          0,          0,          0,             -conj(R1_2),              conj(R1_1),                       0]
];

% Pridict
P_prid_p=F_p*P_in([1:6 end-11:end],[1:6 end-11:end])*F_p'+Q;

% update
K_next=P_prid_p*H'*inv(H*P_prid_p*H'+R);
X_next_p=X_prid_p+K_next*(Z_in([1:3 end-2:end])-Z_prid_p);
RotM_next=RotM_prid*expm(cross_mat(X_next_p(end-2:end)));
P_next_p=(eye(length(X_prid_p))-K_next*H)*P_prid_p;
COMP_T=COMP_T+toc;

X_prid=X_in;
P_prid=P_in;
X_next=X_prid;
P_next=P_prid;
F =eye(length(X_in));
X_prid([1:6 end-11:end])=X_prid_p;
P_prid([1:6 end-11:end],[1:6 end-11:end])=P_prid_p;
P_next([1:6 end-11:end],[1:6 end-11:end])=P_next_p;
F([1:6 end-11:end],[1:6 end-11:end])=F_p;
X_next([1:6 end-11:end])=X_next_p;
end

