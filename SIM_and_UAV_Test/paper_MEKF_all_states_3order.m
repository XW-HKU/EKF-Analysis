%% simulation data generation
clear;clc
syms a b c d g dt real;
[X_gt, pos_in, acc_gryo_in, mag_in, Q_t, Q_gt, Z_gt] = Simulation(0.002, 15);

Z_in=zeros(12,1);
bias(:,1)=zeros(length(Z_in),1);
noise(:,1)=zeros(length(Z_in),1);
for i=2:length(X_gt)
    dt=acc_gryo_in(i,1)-acc_gryo_in(i-1,1);
    dbias=[0.0*randn(3,1);1*randn(3,1);1*randn(3,1);0.0*randn(3,1)];
    bias(:,i)=bias(:,i-1)+dbias*dt;
    %bias(4:9,i)=0.5;
    noise(:,i)=[0.03*randn(3,1);1*randn(3,1);1*randn(3,1);0.02*randn(3,1)];
end
%%
% I-CM estimation
% load noise
% load bias
X2(:,1)=zeros(42,1);X2(9,1)=-9.8;X2(37:39,1)=-0.2;
P2(:,:,1)=eye(length(X2(:,1)));

bias_gy(:,1)=zeros(3,1);

Q_p=0.0000*ones(1,3);
Q_v=0.002^2*0*ones(1,3);
Q_q=0.0000*ones(1,4);
Q_a_s=0*ones(1,3);%q1=100000;q2=100000000;q3=0.00001;q4=0.0000001;
Q_a_s1=10*ones(1,3);
Q_a_s2=8*ones(1,3);
Q_a_s3=0.01*ones(1,3);
Q_omega=0.000*ones(1,3);
Q_tau=1*ones(1,3);
Q_tau1=0.5*ones(1,3);
Q_tau2=0.0000001*ones(1,3);
Q_ba=0.002^2*[0.3 0.5 0.1];
Q_bg=0.002^2*1*ones(1,3);
Q_l_ic=0.000000*ones(1,3);
Q_ptheta=0.002^2*0.*ones(1,3);

Q=diag([Q_p,Q_v,Q_a_s,Q_a_s1,Q_a_s2,Q_a_s3,Q_omega,Q_tau,Q_tau1,Q_tau2,Q_ba,Q_bg,Q_l_ic,Q_ptheta]);
R=diag([0.0005*ones(1,3),1*ones(1,3),1*ones(1,3),0.0004*ones(1,3)]);
RotMatrix(:,:,1)=eye(3);
i=2;

for i=2:length(X_gt)
    dt=acc_gryo_in(i,1)-acc_gryo_in(i-1,1);
    
    % get measurements
    Z_in=[Z_gt([1:3],i);Z_gt([11:13],i);Z_gt([14:16],i);Z_gt([17:19],i)]+bias(:,i)+noise(:,i);
    
    a_s_in1=Z_in(4);a_s_in2=Z_in(5);a_s_in3=Z_in(6);
    omega_in1=Z_in(7);omega_in2=Z_in(8);omega_in3=Z_in(9);
    omega_in=[omega_in1 omega_in2 omega_in3];
    
    % Assignment and input F
    X_in=X2(:,i-1);
    [ p1, p2, p3, v1, v2, v3, a_s1, a_s2, a_s3, a_s_11, a_s_12, a_s_13, a_s_21, a_s_22, a_s_23, a_s_31, a_s_32, a_s_33, omega1, omega2, omega3, tau1, tau2, tau3, tau_11, tau_12, tau_13, tau_21, tau_22, tau_23, ba1, ba2, ba3, bg1, bg2, bg3, l_ic1, l_ic2, l_ic3] = Assignment_MEKF_HighOrder(X_in);
    [ R1_1, R1_2, R1_3, R2_1, R2_2, R2_3, R3_1, R3_2, R3_3 ] = RotMatAssi( RotMatrix(:,:,i-1) );

    F(:,:,i) =[...
[ 1, 0, 0, dt,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 1, 0,  0, dt,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 1,  0,  0, dt,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  1,  0,  0, dt*conj(R1_1), dt*conj(R1_2), dt*conj(R1_3),  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, R1_3*a_s2*dt - R1_2*a_s3*dt, R1_1*a_s3*dt - R1_3*a_s1*dt, R1_2*a_s1*dt - R1_1*a_s2*dt]
[ 0, 0, 0,  0,  1,  0, dt*conj(R2_1), dt*conj(R2_2), dt*conj(R2_3),  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, R2_3*a_s2*dt - R2_2*a_s3*dt, R2_1*a_s3*dt - R2_3*a_s1*dt, R2_2*a_s1*dt - R2_1*a_s2*dt]
[ 0, 0, 0,  0,  0,  1, dt*conj(R3_1), dt*conj(R3_2), dt*conj(R3_3),  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, R3_3*a_s2*dt - R3_2*a_s3*dt, R3_1*a_s3*dt - R3_3*a_s1*dt, R3_2*a_s1*dt - R3_1*a_s2*dt]
[ 0, 0, 0,  0,  0,  0,             1,             0,             0, dt,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             1,             0,  0, dt,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             1,  0,  0, dt,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  1,  0,  0, dt,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  1,  0,  0, dt,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  1,  0,  0, dt,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  1,  0,  0, dt,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  1,  0,  0, dt,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  1,  0,  0, dt,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, dt,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, dt,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, dt,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, dt,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, dt,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, dt,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, dt,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, dt,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, dt, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 1, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 1, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 1, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 1, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 1, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0, dt,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           1,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, dt,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           1,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, dt,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           1]
];

F(end-2:end,end-2:end,i)=expm(cross_mat(-dt*[omega1 omega2 omega3]));%+eye(3);

% predict states and coveriance matrix
    X2_prid(:,i)=[p1 + dt*v1
                                         p2 + dt*v2
                                         p3 + dt*v3
        v1 + dt*(R1_1*a_s1 + R1_2*a_s2 + R1_3*a_s3)
        v2 + dt*(R2_1*a_s1 + R2_2*a_s2 + R2_3*a_s3)
 v3 + dt*(R3_1*a_s1 + R3_2*a_s2 + R3_3*a_s3 + 49/5)
                                   a_s1 + a_s_11*dt
                                   a_s2 + a_s_12*dt
                                   a_s3 + a_s_13*dt
                                 a_s_11 + a_s_21*dt
                                 a_s_12 + a_s_22*dt
                                 a_s_13 + a_s_23*dt
                                 a_s_21 + a_s_31*dt
                                 a_s_22 + a_s_32*dt
                                 a_s_23 + a_s_33*dt
                                             a_s_31
                                             a_s_32
                                             a_s_33
                                   omega1 + dt*tau1
                                   omega2 + dt*tau2
                                   omega3 + dt*tau3
                                   tau1 + dt*tau_11
                                   tau2 + dt*tau_12
                                   tau3 + dt*tau_13
                                 tau_11 + dt*tau_21
                                 tau_12 + dt*tau_22
                                 tau_13 + dt*tau_23
                                             tau_21
                                             tau_22
                                             tau_23
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
    RotMatrix(:,:,i)=RotMatrix(:,:,i-1)*expm(cross_mat(dt*[omega1 omega2 omega3]));
    % Predict output
    [ p1, p2, p3, v1, v2, v3, a_s1, a_s2, a_s3, a_s_11, a_s_12, a_s_13, a_s_21, a_s_22, a_s_23, a_s_31, a_s_32, a_s_33, omega1, omega2, omega3, tau1, tau2, tau3, tau_11, tau_12, tau_13, tau_21, tau_22, tau_23, ba1, ba2, ba3, bg1, bg2, bg3, l_ic1, l_ic2, l_ic3] = Assignment_MEKF_HighOrder(X_in);
    [ R1_1, R1_2, R1_3, R2_1, R2_2, R2_3, R3_1, R3_2, R3_3 ] = RotMatAssi(RotMatrix(:,:,i));
    Z_prid(:,i)=[p1 + R1_1*l_ic1 + R1_2*l_ic2 + R1_3*l_ic3
 p2 + R2_1*l_ic1 + R2_2*l_ic2 + R2_3*l_ic3
 p3 + R3_1*l_ic1 + R3_2*l_ic2 + R3_3*l_ic3
                                a_s1 + ba1
                                a_s2 + ba2
                                a_s3 + ba3
                              bg1 + omega1
                              bg2 + omega2
                              bg3 + omega3
                                      R1_1
                                      R2_1
                                      R3_1];
 
    H2(:,:,i) = [...
[ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, conj(R1_1), conj(R1_2), conj(R1_3), R1_3*l_ic2 - R1_2*l_ic3, R1_1*l_ic3 - R1_3*l_ic1, R1_2*l_ic1 - R1_1*l_ic2]
[ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, conj(R2_1), conj(R2_2), conj(R2_3), R2_3*l_ic2 - R2_2*l_ic3, R2_1*l_ic3 - R2_3*l_ic1, R2_2*l_ic1 - R2_1*l_ic2]
[ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, conj(R3_1), conj(R3_2), conj(R3_3), R3_3*l_ic2 - R3_2*l_ic3, R3_1*l_ic3 - R3_3*l_ic1, R3_2*l_ic1 - R3_1*l_ic2]
[ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,          0,          0,          0,                       0,                       0,                       0]
[ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,          0,          0,          0,                       0,                       0,                       0]
[ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,          0,          0,          0,                       0,                       0,                       0]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,          0,          0,          0,                       0,                       0,                       0]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,          0,          0,          0,                       0,                       0,                       0]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,          0,          0,          0,                       0,                       0,                       0]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,          0,          0,          0,                       0,                   -R1_3,                    R1_2]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,          0,          0,          0,                       0,                   -R2_3,                    R2_2]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,          0,          0,          0,                       0,                   -R3_3,                    R3_2]
];

    P_in=P2(:,:,i-1);
    F_w=eye(length(Q));
    %F_w(4:6,4:6)=-RotMatrix(:,:,i);
    F_w(end-2:end,end-2:end)=-eye(3);
    P2_prid(:,:,i)=F(:,:,i)*P_in*F(:,:,i)'+F_w*Q*F_w';
    
% update
    K_next=P2_prid(:,:,i)*H2(:,:,i)'*inv(H2(:,:,i)*P_in*H2(:,:,i)'+R);
    X_next=X2_prid(:,i)+K_next*(Z_in-Z_prid(:,i));
    RotMatrix(:,:,i)=RotMatrix(:,:,i)*expm(cross_mat(X_next(end-2:end)));
    det(RotMatrix(:,:,i));
    L=eye(length(X2_prid(:,i)));
    L(end-2:end,end-2:end)=A(X_next(end-2:end))';
    P_next=L*(eye(length(X_in))-K_next*H2(:,:,i))*P2_prid(:,:,i)*L';
    
    as_B_f(:,i)=X_next(7:9);%as_B_f(:,i)=RotMatrix(:,:,i)'*X_next(7:9);
    m_p_f(:,i)=X_next(1:3)+RotMatrix(:,:,i)*X_next(end-5:end-3);
    Euler_sim(:,i)=rotm2eul(RotMatrix(:,:,i));
    
    X2(:,i)=X_next;
    P2(:,:,i)=P_next;
end

 N=length(X2);
latex_font=12

%%
latex_font=12;
[b,a] = butter(2,40/850);
d1 = designfilt('lowpassiir','FilterOrder',1, ...
    'HalfPowerFrequency',40/850,'DesignMethod','butter');
%freqz(b,a);
POS_raw=Z_gt(1:3,:)+noise(1:3,:);
IMU_raw=Z_gt(11:16,:)+bias(4:9,:)+noise(4:9,:);
imu_butfilt(1,:) = filter(b,a,IMU_raw(1,:));
imu_butfilt(2,:) = filter(b,a,IMU_raw(2,:));
imu_butfilt(3,:) = filter(b,a,IMU_raw(3,:));

imu_zerofilt(1,:)=filtfilt(d1,IMU_raw(1,:));
imu_zerofilt(2,:)=filtfilt(d1,IMU_raw(2,:));
imu_zerofilt(3,:)=filtfilt(d1,IMU_raw(3,:));

[b,a] = butter(2,20/850);
d1 = designfilt('lowpassiir','FilterOrder',1, ...
    'HalfPowerFrequency',20/850,'DesignMethod','butter');
imu_butfilt(4,:) = filter(b,a,IMU_raw(4,:));
imu_butfilt(5,:) = filter(b,a,IMU_raw(5,:));
imu_butfilt(6,:) = filter(b,a,IMU_raw(6,:));
imu_zerofilt(4,:)=filtfilt(d1,IMU_raw(4,:));
imu_zerofilt(5,:)=filtfilt(d1,IMU_raw(5,:));
imu_zerofilt(6,:)=filtfilt(d1,IMU_raw(6,:));

f1=figure(11)

subplot(2,1,1)
plot(acc_gryo_in(:,1),IMU_raw(2,:),'color',[150 150 150]/255,'linewidth',0.3);hold on;
plot(acc_gryo_in(:,1),Z_gt(12,:),'k','linewidth',1.5);hold on;
plot(acc_gryo_in(:,1),imu_butfilt(2,:),'r:','linewidth',1.5);hold on;
plot(acc_gryo_in(:,1),imu_zerofilt(2,:),'g-.','linewidth',1.5);hold on;
plot(acc_gryo_in(:,1),as_B_f(2,:),'b:','linewidth',1.5);hold on;
%xlim([10.8,12.2]);
title('Body-X Acceleration','interpreter','latex','FontSize', latex_font);
ylabel('$m^2/s$','interpreter','latex','FontSize', latex_font);
legend('Measurement','Ground Truth','Normal Low-Pass','Zero Phase Low-Pss','Our Method');%title('Pos Y');

subplot(2,1,2)
plot(acc_gryo_in(:,1),IMU_raw(6,:),'color',[150 150 150]/255,'linewidth',0.3);hold on;
plot(acc_gryo_in(:,1),Z_gt(16,:),'k','linewidth',1.5);hold on;
plot(acc_gryo_in(:,1),imu_butfilt(6,:),'r:','linewidth',1.5);hold on;
plot(acc_gryo_in(:,1),imu_zerofilt(6,:),'g-.','linewidth',1.5);hold on;
plot(acc_gryo_in(:,1),X2(21,:),'b:','linewidth',1.5);hold on;
%xlim([10.8,12.2]);
%plot(acc_gryo_in(:,1),X2_smooth(20,:));
title('Body-X Angular Velocity','interpreter','latex','FontSize', latex_font);
xlabel('time (s)','interpreter','latex','FontSize', latex_font);
ylabel('$rad/s$','interpreter','latex','FontSize', latex_font);
legend('Measurement','Ground Truth','Normal Low-Pass','Zero Phase Low-Pss','Our Method');%title('Pos Y');

%%
f1=figure(1)
%set(gcf,'outerposition',get(0,'screensize'));
latex_font=12
index=[round(0.5*length(Z_gt)):length(Z_gt)];

ha = tight_subplot(3,2,[.07 .08],[0.05 .05],[.09 .008])
% p
axes(ha(1));
plot(acc_gryo_in(:,1),Z_gt(1,:),'k','linewidth',1);hold on;
plot(acc_gryo_in(:,1),m_p_f(1,:),'r:','linewidth',1);hold on;
title('Position Estimation Results','interpreter','latex','FontSize', latex_font)
ylabel('Position X($m$)','interpreter','latex','FontSize', latex_font);
%ylim([-0.25,0.25]);
legend('Ground Truth','Our Method','Zero Phase Low-Pss','Our Method');

axes(ha(2));
[errorNum,xaxis] = hist(m_p_f(1,index)-Z_gt(1,index),[-0.1:0.001:0.1]); 
bar(xaxis,errorNum);hold on;


axes(ha(3));
plot(acc_gryo_in(:,1),Z_gt(2,:),'k','linewidth',1);hold on;
plot(acc_gryo_in(:,1),m_p_f(2,:),'r:','linewidth',1);hold on;
%title('Position Estimation Results','interpreter','latex','FontSize', latex_font)
ylabel('Position Y($m$)','interpreter','latex','FontSize', latex_font);
%ylim([-0.25,0.25]);
legend('Ground Truth','Our Method','Zero Phase Low-Pss','Our Method');
axes(ha(4));
[errorNum,xaxis] = hist(m_p_f(2,index)-Z_gt(2,index),[-0.1:0.001:0.1]); 
bar(xaxis,errorNum);hold on;

axes(ha(5));
plot(acc_gryo_in(:,1),-Z_gt(4,:),'k','linewidth',1);hold on;
plot(acc_gryo_in(:,1),-m_v_f(1,:),'r:','linewidth',1);hold on;
%title('Position Estimation Results','interpreter','latex','FontSize', latex_font)
ylabel('Position Z($m$)','interpreter','latex','FontSize', latex_font);
%ylim([-0.25,0.25]);
legend('Ground Truth','Our Method','Zero Phase Low-Pss','Our Method');
axes(ha(6));
[errorNum,xaxis] = hist(m_p_f(3,index)-Z_gt(3,index),[-0.1:0.001:0.1]); 
bar(xaxis,errorNum);hold on;


%%

subplot(6,2,3)
plot(acc_gryo_in(:,1),Z_gt(2,:),'k','linewidth',1);hold on;
plot(acc_gryo_in(:,1),m_p_f(2,:),'r:','linewidth',1);hold on;
title('Position Estimation Results','interpreter','latex','FontSize', latex_font)
ylabel('Position Y($m$)','interpreter','latex','FontSize', latex_font);
%ylim([-0.25,0.25]);
legend('Ground Truth','Our Method','Zero Phase Low-Pss','Our Method');
subplot(6,2,4)
[errorNum,xaxis] = hist(m_p_f(2,index)-Z_gt(2,index),[-0.1:0.001:0.1]); 
bar(xaxis,errorNum);hold on;

subplot(6,2,5)
%plot(acc_gryo_in(:,1),-POS_raw(3,:),'color',[150 150 150]/255,'linewidth',0.3);hold on;
plot(acc_gryo_in(:,1),-Z_gt(3,:),'k','linewidth',1);hold on;
plot(acc_gryo_in(:,1),-m_p_f(3,:),'r:','linewidth',1);hold on;
%title('Position Estimation Results','interpreter','latex','FontSize', latex_font)
ylabel('Position Z($m$)','interpreter','latex','FontSize', latex_font);
%ylim([-0.25,0.25]);
legend('Ground Truth','Our Method','Zero Phase Low-Pss','Our Method');
subplot(6,2,6)
[errorNum,xaxis] = hist(m_p_f(3,index)-Z_gt(3,index),[-0.1:0.001:0.1]); 
bar(xaxis,errorNum);hold on;

subplot(6,2,7)
%plot(acc_gryo_in(:,1),IMU_raw(1,:),'color',[150 150 150]/255,'linewidth',0.3);hold on;
plot(acc_gryo_in(:,1),Z_gt(11,:),'k','linewidth',1);hold on;
plot(acc_gryo_in(:,1),as_B_f(1,:),'r:','linewidth',1);hold on;
% plot(acc_gryo_in(:,1),as_B_f(1,:)-Z_gt(11,:));hold on;
% plot(acc_gryo_in(:,1),as_B_f(2,:)-Z_gt(12,:));hold on;
% plot(acc_gryo_in(:,1),as_B_f(3,:)-Z_gt(13,:));
ylim([-3,3]);
title('Acceleration Estimation Error','interpreter','latex','FontSize', latex_font);
ylabel('$m^2/s$','interpreter','latex','FontSize', latex_font);
legend('Ground Truth','Our Method','Zero Phase Low-Pss','Our Method');
subplot(6,2,8)
[errorNum,xaxis] = hist(as_B_f(1,index)-Z_gt(11,index),[-2:0.02:2]);
bar(xaxis,errorNum);hold on;

subplot(6,1,4)
plot(acc_gryo_in(:,1),X2(19,:)-X_gt(11,:));hold on;
plot(acc_gryo_in(:,1),X2(20,:)-X_gt(12,:));hold on;
plot(acc_gryo_in(:,1),X2(21,:)-X_gt(13,:));hold on;
ylim([-3,3]);
%plot(acc_gryo_in(:,1),X2_smooth(20,:));
title('Angular Velocity Estimation Error','interpreter','latex','FontSize', latex_font);
ylabel('$rad/s$','interpreter','latex','FontSize', latex_font);
legend('Body-X','Body-Y','Body-Z');%title('omega Y');

subplot(6,1,5)
eul_err(:,1)=Euler_gt(:,1)-interp1(acc_gryo_in(:,1),Euler_sim(1,:)',Q_gt(:,1));
eul_err(:,2)=Euler_gt(:,2)-interp1(acc_gryo_in(:,1),Euler_sim(2,:)',Q_gt(:,1));
eul_err(:,3)=Euler_gt(:,3)-interp1(acc_gryo_in(:,1),Euler_sim(3,:)',Q_gt(:,1));
plot(Q_gt(:,1),eul_err(:,1));hold on;
plot(Q_gt(:,1),eul_err(:,2));hold on;
plot(Q_gt(:,1),eul_err(:,3));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(20,:));
title('Attitude Estimation Error','interpreter','latex','FontSize', latex_font);
ylabel('$rad$','interpreter','latex','FontSize', latex_font);
legend('Yaw','Pitch','Roll');%title('omega Y');
ylim([-0.3,0.3]);

subplot(6,1,6)
plot(acc_gryo_in(:,1),X2(37,:)+0.5);hold on;
plot(acc_gryo_in(:,1),-X2(38,:)-0.5);hold on;
plot(acc_gryo_in(:,1),X2(39,:)+0.5);hold on;
%plot(acc_gryo_in(:,1),X2_smooth(20,:));
title('Transitional Offset Estimation Error','interpreter','latex','FontSize', latex_font);
ylabel('$m$','interpreter','latex','FontSize', latex_font);
legend('Body-X','Body-Y','Body-Z');%title('omega Y');
ylim([-1,1]);

%%
f1=figure(2)
%a_s
subplot(3,3,4)
plot(acc_gryo_in(:,1),Z_gt(11,:)+bias(4,:)+noise(4,:));hold on;
plot(acc_gryo_in(:,1),as_B_f(1,:));hold on;
plot(acc_gryo_in(:,1),Z_gt(11,:));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(7,:));
legend('measure','filtered','gt-IMU');title('acc X');
subplot(3,3,5)
plot(acc_gryo_in(:,1),Z_gt(12,:)+bias(5,:)+noise(5,:));hold on;
plot(acc_gryo_in(:,1),as_B_f(2,:));hold on;
plot(acc_gryo_in(:,1),Z_gt(12,:));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(9,:));
legend('measure','filtered','gt-IMU');title('acc Y');
subplot(3,3,6)
plot(acc_gryo_in(:,1),Z_gt(13,:)+bias(6,:)+noise(6,:));hold on;
plot(acc_gryo_in(:,1),as_B_f(3,:));hold on;
plot(acc_gryo_in(:,1),Z_gt(13,:));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(9,:));
legend('measure','filtered','gt-IMU' );title('acc Z');

%omega
subplot(3,3,7)
plot(acc_gryo_in(:,1),X_gt(11,:)+bias(7,:)+noise(7,:));hold on;
plot(acc_gryo_in(:,1),X2(19,:));hold on;
plot(acc_gryo_in(:,1),X_gt(11,:));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(19,:));
legend('measure','filtered','gt');title('omega X');
subplot(3,3,8)
plot(acc_gryo_in(:,1),X_gt(12,:)+bias(8,:)+noise(8,:));hold on;
plot(acc_gryo_in(:,1),X2(20,:));hold on;
plot(acc_gryo_in(:,1),X_gt(12,:));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(20,:));
legend('measure','filtered','gt' );title('omega Y');
subplot(3,3,9)
plot(acc_gryo_in(:,1),X_gt(13,:)+bias(9,:)+noise(9,:));hold on;
plot(acc_gryo_in(:,1),X2(21,:));hold on;
plot(acc_gryo_in(:,1),X_gt(13,:));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(21,:));
legend('measure','filtered','gt' );title('omega Z');

f2=figure(3)
set(gcf,'outerposition',get(0,'screensize'));

%euler angle
subplot(4,3,1)
plot(Q_gt(:,1),Euler_gt(:,1).*57.3);hold on;
plot(acc_gryo_in(:,1),Euler_sim(1,:).*57.3);hold on;
%plot(acc_gryo_in(:,1),Euler_smooth(1,:).*57.3);
legend('gt','filtered' );title('yaw');
subplot(4,3,2)
plot(Q_gt(:,1),Euler_gt(:,2).*57.3);hold on;
plot(acc_gryo_in(:,1),Euler_sim(2,:).*57.3);hold on;
%plot(acc_gryo_in(:,1),Euler_smooth(2,:).*57.3);
legend('gt','filtered' );title('pitch');
subplot(4,3,3)
plot(Q_gt(:,1),Euler_gt(:,3).*57.3);hold on;
plot(acc_gryo_in(:,1),Euler_sim(3,:).*57.3);hold on;
%plot(acc_gryo_in(:,1),Euler_smooth(3,:).*57.3);
legend('gt','filtered' );title('roll');

% bias_acc
subplot(4,3,4)
plot(acc_gryo_in(:,1),bias(4,:));hold on;
plot(acc_gryo_in(:,1),X2(31,:));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(22,:));
legend('gt','filtered' );title('Acc Bias X');
subplot(4,3,5)
plot(acc_gryo_in(:,1),bias(5,:));hold on;
plot(acc_gryo_in(:,1),X2(32,:));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(23,:));
legend('gt','filtered' );title('Acc Bias Y');
subplot(4,3,6)
plot(acc_gryo_in(:,1),bias(6,:));hold on;
plot(acc_gryo_in(:,1),X2(33,:));
%plot(acc_gryo_in(:,1),X2_smooth(24,:));
legend('gt','filtered' );title('Acc Bias Z');

% bias_gyro
subplot(4,3,7)
plot(acc_gryo_in(:,1),bias(7,:));hold on;
plot(acc_gryo_in(:,1),X2(34,:));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(25,:));
legend('gt','filtered' );title('Gyro Bias X');
subplot(4,3,8)
plot(acc_gryo_in(:,1),bias(8,:));hold on;
plot(acc_gryo_in(:,1),X2(35,:));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(26,:));
legend('gt','filtered' );title('Gyro Bias Y');
subplot(4,3,9)
plot(acc_gryo_in(:,1),bias(9,:));hold on;
plot(acc_gryo_in(:,1),X2(36,:));
%plot(acc_gryo_in(:,1),X2_smooth(27,:));
legend('gt','filtered' );title('Gyro Bias Z');

% L_ic
subplot(4,3,10)
plot(acc_gryo_in(:,1),-X_gt(47,:));hold on;
plot(acc_gryo_in(:,1),X2(37,:));hold on;
plot(acc_gryo_in(:,1),X2_smooth(28,:));
legend('gt','filtered' );title('L IC X');
subplot(4,3,11)
plot(acc_gryo_in(:,1),-X_gt(48,:));hold on;
plot(acc_gryo_in(:,1),X2(38,:));hold on;
plot(acc_gryo_in(:,1),X2_smooth(29,:));
legend('gt','filtered' );title('L IC Y');
subplot(4,3,12)
plot(acc_gryo_in(:,1),-X_gt(49,:));hold on;
plot(acc_gryo_in(:,1),X2(39,:));
plot(acc_gryo_in(:,1),X2_smooth(30,:));
legend('gt','filtered' );title('L IC Z');
%%
f1=figure(2)
subplot(2,1,1)
plot(acc_gryo_in(:,1),Z_gt(1,:));hold on;
plot(acc_gryo_in(:,1),Z_gt(2,:));hold on;
plot(acc_gryo_in(:,1),-Z_gt(3,:));hold on;
title('Position','interpreter','latex','FontSize', latex_font)
ylabel('$m$','interpreter','latex','FontSize', latex_font);
ylim([-2,10]);
legend('Inertial-X','Inertial-Y','Inertial-Z');

subplot(2,1,2)
plot(Q_gt(:,1),2*Euler_gt(:,1));hold on;
plot(Q_gt(:,1),2*Euler_gt(:,2));hold on;
plot(Q_gt(:,1),2*Euler_gt(:,3));
ylim([-1,1]);
title('Attitude','interpreter','latex','FontSize', latex_font);
ylabel('$rad$','interpreter','latex','FontSize', latex_font);
legend('Yaw','Pitch','Roll');%title('Pos Y');

