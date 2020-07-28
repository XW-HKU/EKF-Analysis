%% simulation data generation
clear;clc
syms a b c d g dt real;
[X_gt, pos_in, acc_gryo_in, mag_in, Q_t, Q_gt, Z_gt] = Simulation(0.002, 25);

Z_in=zeros(12,1);
bias(:,1)=zeros(length(Z_in),1);
noise(:,1)=zeros(length(Z_in),1);
for i=2:length(X_gt)
    dt=acc_gryo_in(i,1)-acc_gryo_in(i-1,1);
    dbias=[0.0*randn(3,1);2*randn(3,1);1*randn(3,1);0.0*randn(3,1)];
    bias(:,i)=bias(:,i-1)+dbias*dt;
    %bias(4:9,i)=0.5;
    noise(:,i)=[0.05*randn(3,1);2*randn(3,1);2*randn(3,1);0.2*randn(3,1)];
end

%% I-CM estimation
% load noise
% load bias
X2(:,1)=zeros(33,1);X2(9,1)=-9.8;
P2(:,:,1)=eye(length(X2(:,1)));

bias_gy(:,1)=zeros(3,1);

Q_p=0.0000*ones(1,3);
Q_v=0.002^2*0*ones(1,3);
Q_q=0.0000*ones(1,4);
Q_a_s=0.02*ones(1,3);
Q_a_s1=0*ones(1,3);
Q_a_s2=0*ones(1,3);
Q_a_s3=0.0*ones(1,3);
Q_omega=0.005*ones(1,3);
Q_ba=0.002^2*5*ones(1,3);
Q_bg=0.002^2*2*ones(1,3);
Q_l_ic=0.000000*ones(1,3);
Q_ptheta=0.002^2*0.*ones(1,3);

Q=diag([Q_p,Q_v,Q_a_s,Q_a_s1,Q_a_s2,Q_a_s3,Q_omega,Q_ba,Q_bg,Q_l_ic,Q_ptheta]);
R=diag([0.05*ones(1,3),2*ones(1,3),2*ones(1,3),0.2*ones(1,3)]);
RotMatrix(:,:,1)=eye(3);
i=2;

for i=2:length(X_gt)
    dt=acc_gryo_in(i,1)-acc_gryo_in(i-1,1);
    
    % generate bias and noise
    Z_in=[Z_gt([1:3],i);Z_gt([11:13],i);Z_gt([14:16],i);Z_gt([17:19],i)]+bias(:,i)+noise(:,i);
    
    %dbias_gy=0*randn(3,1);
    %bias_gy(:,i)=bias_gy(:,i-1)+dbias_gy*dt;
    %acc_gryo_in(i,[5:7])=acc_gryo_in(i,[5:7])+bias_gy(:,i)';
    
    a_s_in1=Z_in(4);a_s_in2=Z_in(5);a_s_in3=Z_in(6);
    omega_in1=Z_in(7);omega_in2=Z_in(8);omega_in3=Z_in(9);
    omega_in=[omega_in1 omega_in2 omega_in3];
    
    % Assignment and calculate F
    X_in=X2(:,i-1);
    [p1, p2, p3, v1, v2, v3, omega1, omega2, omega3, a_s1, a_s2, a_s3, a_s_11, a_s_12, a_s_13, a_s_21, a_s_22, a_s_23, a_s_31, a_s_32, a_s_33, tau1, tau2, tau3, ba1, ba2, ba3, bg1, bg2, bg3, l_ic1, l_ic2, l_ic3] = Assignment_MEKF(X_in,X_in,acc_gryo_in(i-1,[5:7]));
    [ R1_1, R1_2, R1_3, R2_1, R2_2, R2_3, R3_1, R3_2, R3_3 ] = RotMatAssi( RotMatrix(:,:,i-1) );
    bg=[bg1 bg2 bg3];
    F(:,:,i) =[...
[ 1, 0, 0, dt,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 1, 0,  0, dt,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 1,  0,  0, dt,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  1,  0,  0, dt*conj(R1_1), dt*conj(R1_2), dt*conj(R1_3), 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, R1_3*a_s2*dt - R1_2*a_s3*dt, R1_1*a_s3*dt - R1_3*a_s1*dt, R1_2*a_s1*dt - R1_1*a_s2*dt]
[ 0, 0, 0,  0,  1,  0, dt*conj(R2_1), dt*conj(R2_2), dt*conj(R2_3), 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, R2_3*a_s2*dt - R2_2*a_s3*dt, R2_1*a_s3*dt - R2_3*a_s1*dt, R2_2*a_s1*dt - R2_1*a_s2*dt]
[ 0, 0, 0,  0,  0,  1, dt*conj(R3_1), dt*conj(R3_2), dt*conj(R3_3), 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, R3_3*a_s2*dt - R3_2*a_s3*dt, R3_1*a_s3*dt - R3_3*a_s1*dt, R3_2*a_s1*dt - R3_1*a_s2*dt]
[ 0, 0, 0,  0,  0,  0,             1,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             1,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  1,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 1, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 1, 0, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 1, 0, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 1, 0, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 1, 0, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1,                           0,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, dt,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           1,                           0,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, dt,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           1,                           0]
[ 0, 0, 0,  0,  0,  0,             0,             0,             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, dt, 0, 0, 0, 0, 0, 0, 0, 0, 0,                           0,                           0,                           1]
];

F(end-2:end,end-2:end,i)=expm(cross_mat(-dt*[omega1 omega2 omega3]));%+eye(3);
F(4:6,end-2:end)=-dt*RotMatrix(:,:,i-1)*cross_mat([a_s1 a_s2 a_s3]);
 %F(end-2:end,19:21,i)=dt*eye(3);

% predict states and coveriance matrix
    X2_prid(:,i)=[p1 + dt*v1
                                         p2 + dt*v2
                                         p3 + dt*v3
        v1 + dt*(R1_1*a_s1 + R1_2*a_s2 + R1_3*a_s3)
        v2 + dt*(R2_1*a_s1 + R2_2*a_s2 + R2_3*a_s3)
 v3 + dt*(R3_1*a_s1 + R3_2*a_s2 + R3_3*a_s3 + 49/5)
                                               a_s1
                                               a_s2
                                               a_s3
                                             a_s_11
                                             a_s_12
                                             a_s_13
                                             a_s_21
                                             a_s_22
                                             a_s_23
                                             a_s_31
                                             a_s_32
                                             a_s_33
                                             omega1
                                             omega2
                                             omega3
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
    [p1, p2, p3, v1, v2, v3, omega1, omega2, omega3, a_s1, a_s2, a_s3, a_s_11, a_s_12, a_s_13, a_s_21, a_s_22, a_s_23, a_s_31, a_s_32, a_s_33, tau1, tau2, tau3, ba1, ba2, ba3, bg1, bg2, bg3, l_ic1, l_ic2, l_ic3] = Assignment_MEKF(X2_prid(:,i),X_in,acc_gryo_in(i,[5:7]));
    [ R1_1, R1_2, R1_3, R2_1, R2_2, R2_3, R3_1, R3_2, R3_3 ] = RotMatAssi( RotMatrix(:,:,i) );
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
[ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, conj(R1_1), conj(R1_2), conj(R1_3), R1_3*l_ic2 - R1_2*l_ic3, R1_1*l_ic3 - R1_3*l_ic1, R1_2*l_ic1 - R1_1*l_ic2]
[ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, conj(R2_1), conj(R2_2), conj(R2_3), R2_3*l_ic2 - R2_2*l_ic3, R2_1*l_ic3 - R2_3*l_ic1, R2_2*l_ic1 - R2_1*l_ic2]
[ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, conj(R3_1), conj(R3_2), conj(R3_3), R3_3*l_ic2 - R3_2*l_ic3, R3_1*l_ic3 - R3_3*l_ic1, R3_2*l_ic1 - R3_1*l_ic2]
[ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,          0,          0,          0,                       0,                       0,                       0]
[ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,          0,          0,          0,                       0,                       0,                       0]
[ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,          0,          0,          0,                       0,                       0,                       0]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0,          0,          0,          0,                       0,                       0,                       0]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,          0,          0,          0,                       0,                       0,                       0]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,          0,          0,          0,                       0,                       0,                       0]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,          0,          0,          0,                       0,                   -R1_3,                    R1_2]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,          0,          0,          0,                       0,                   -R2_3,                    R2_2]
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,          0,          0,          0,                       0,                   -R3_3,                    R3_2]
];

%H2(:,22:27,i)=0;
    P_in=P2(:,:,i-1);
    F_w=eye(length(Q));
    %F_w(4:6,4:6)=-RotMatrix(:,:,i);
    F_w(end-2:end,end-2:end)=-eye(3);
    P2_prid(:,:,i)=F(:,:,i)*P_in*F(:,:,i)'+F_w*Q*F_w';  
% update
    K_next=P2_prid(:,:,i)*H2(:,:,i)'*inv(H2(:,:,i)*P_in*H2(:,:,i)'+R);
    X_next=X2_prid(:,i)+K_next*([Z_in(1:3);Z_in(4:6);Z_in(7:9);Z_in(10:12)]-Z_prid(:,i));
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
 X2_smooth(:,[N-1 N])=X2(:,[N-1,N]);
% P2_smooth(:,:,[N-1 N])=P2(:,:,[N-1,N]);
for j=1:length(Q_gt)
    Euler_gt(j,:)=quat2eul(Q_gt(j,[2:5]));
end
Euler_smooth(:,[N-1 N])=quat2eul(X2_smooth([7:10],[N-1,N])')';

f1=figure(2)
set(gcf,'outerposition',get(0,'screensize'));
% p
subplot(3,3,1)
plot(acc_gryo_in(:,1),Z_gt(1,:)+noise(1,:));hold on;
plot(acc_gryo_in(:,1),m_p_f(1,:));hold on;
plot(acc_gryo_in(:,1),Z_gt(1,:));
legend('measure','filtered','gt');title('Pos X');
subplot(3,3,2)
plot(acc_gryo_in(:,1),Z_gt(2,:)+noise(2,:));hold on;
plot(acc_gryo_in(:,1),m_p_f(2,:));hold on;
plot(acc_gryo_in(:,1),Z_gt(2,:));
legend('measure','filtered','gt');title('Pos Y');
subplot(3,3,3)
plot(acc_gryo_in(:,1),Z_gt(3,:)+noise(3,:));hold on;
plot(acc_gryo_in(:,1),m_p_f(3,:));
plot(acc_gryo_in(:,1),Z_gt(3,:));
legend('measure','filtered','gt');title('Pos Z');

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

% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(gcf,'1a','-dpdf','-r600')

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
plot(acc_gryo_in(:,1),X2(22,:));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(22,:));
legend('gt','filtered' );title('Acc Bias X');
subplot(4,3,5)
plot(acc_gryo_in(:,1),bias(5,:));hold on;
plot(acc_gryo_in(:,1),X2(23,:));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(23,:));
legend('gt','filtered' );title('Acc Bias Y');
subplot(4,3,6)
plot(acc_gryo_in(:,1),bias(6,:));hold on;
plot(acc_gryo_in(:,1),X2(24,:));
%plot(acc_gryo_in(:,1),X2_smooth(24,:));
legend('gt','filtered' );title('Acc Bias Z');

% bias_gyro
subplot(4,3,7)
plot(acc_gryo_in(:,1),bias(7,:));hold on;
plot(acc_gryo_in(:,1),X2(25,:));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(25,:));
legend('gt','filtered' );title('Gyro Bias X');
subplot(4,3,8)
plot(acc_gryo_in(:,1),bias(8,:));hold on;
plot(acc_gryo_in(:,1),X2(26,:));hold on;
%plot(acc_gryo_in(:,1),X2_smooth(26,:));
legend('gt','filtered' );title('Gyro Bias Y');
subplot(4,3,9)
plot(acc_gryo_in(:,1),bias(9,:));hold on;
plot(acc_gryo_in(:,1),X2(27,:));
%plot(acc_gryo_in(:,1),X2_smooth(27,:));
legend('gt','filtered' );title('Gyro Bias Z');

% L_ic
subplot(4,3,10)
plot(acc_gryo_in(:,1),-X_gt(47,:));hold on;
plot(acc_gryo_in(:,1),X2(28,:));hold on;
plot(acc_gryo_in(:,1),X2_smooth(28,:));
legend('gt','filtered' );title('L IC X');
subplot(4,3,11)
plot(acc_gryo_in(:,1),-X_gt(48,:));hold on;
plot(acc_gryo_in(:,1),X2(29,:));hold on;
plot(acc_gryo_in(:,1),X2_smooth(29,:));
legend('gt','filtered' );title('L IC Y');
subplot(4,3,12)
plot(acc_gryo_in(:,1),-X_gt(49,:));hold on;
plot(acc_gryo_in(:,1),X2(30,:));
plot(acc_gryo_in(:,1),X2_smooth(30,:));
legend('gt','filtered' );title('L IC Z');

%suptitle('use MEKF acc and angular rate as states')
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(gcf,'1b','-dpdf','-r600')