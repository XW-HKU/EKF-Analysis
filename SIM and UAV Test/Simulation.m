function [X_gt,pos_in,acc_gryo_in,mag_in,Q_t,Q_gt,Z_gt]=Simulation(dt, T)
%clc;clear;dt=0.002;T=20;
X_gt(:,1)=[zeros(6,1);[1 0 0 0]';zeros(3,1);[0 0 -9.8]';zeros(30,1);[1 1 1]';9.8];
Z_gt(:,1)=zeros(19,1);Z_gt(7,1)=1;
Euler_sim(:,1)=[0 0 0]';
Q_p=0.00*ones(1,3);
Q_v=0.0*ones(1,3);
Q_q=0.00*ones(1,4);
Q_omega=0.002*ones(1,3);
Q_a_s=0.001*ones(1,3);
Q_a_s_1=0.0001*ones(1,3);
Q_a_s_2=0.0001*ones(1,3);
Q_a_s_3=0.001*ones(1,3);
Q_tau=0.0001*ones(1,3);
Q_tau_1=0.001*ones(1,3);
Q_tau_2=0.001*ones(1,3);
Q_tau_3=0.001*ones(1,3);
Q_ba=0*ones(1,3);
Q_bg=0*ones(1,3);
Q_bm=0*ones(1,3);
Q_l_ic=0.00*ones(1,3);
Q_g=0;

Q=diag([Q_p,Q_v,Q_q,Q_omega,Q_a_s,Q_tau,Q_l_ic]);
R=diag([0.1*ones(1,3),0.0*ones(1,3),0.5*ones(1,4),0.5*ones(1,3),0.5*ones(1,3),0.5*ones(1,3)]);

i=2;

for t=[dt:dt:T]
    [ p1, p2, p3, v1, v2, v3, a, b, c, d, omega1, omega2, omega3, a_s1, a_s2, a_s3, a_s_11, a_s_12, a_s_13, a_s_21, a_s_22, a_s_23, a_s_31, a_s_32, a_s_33, tau1, tau2, tau3, tau_11, tau_12, tau_13, tau_21, tau_22, tau_23, tau_31, tau_32, tau_33, ba1, ba2, ba3, bg1, bg2, bg3, bm1, bm2, bm3, l_ic1, l_ic2, l_ic3, g] = Assignment(X_gt(:,i-1));

    % controller
    if t<1
        z_sp=0;
        acc_bz_in=-9.8;
        tau_1_in=0;
        tau_2_in=0;
        tau_3_in=0;
    else if t<=3
        z_sp=-5/2*(t-1);
        else if t<=12
                z_sp=-5;
            else if t<=14
                    z_sp=-5+5/2*(t-12);
                else
                    z_sp=0;
                end
            end
        end
    end
    if t>=1 && t<=14
        tt=t-1;
        acc_bz_in=-9.8+1*sin(5.5*tt);
        tau_1_in=9*sin(6*tt)+3*sin(17*tt);
        tau_2_in=9.2*sin(5*tt)+3*sin(19.4*tt);
        tau_3_in=8.3*sin(6.5*tt)+4*sin(14*tt);
    end
    
    %euler setpoint
	euler_sp(1,i)=-v2*0.03-p2*0.03;
	euler_sp(2,i)=v1*0.03+p1*0.03;
	euler_sp(3,i)=0;
    
    q_sp=angle2quat(euler_sp(3,i),euler_sp(2,i),euler_sp(1,i));
    
    ang_axis=quatmultiply(quatinv(q_sp),X_gt([7:10],i-1)');
    omega_sp(:,i)=-ang_axis([2:4])*8;
    
    %acc_bz_ctrl_int=acc_bz_ctrl_int+dt*Euler_sim(1,i-1)
%     tau_1_ctrl_int=acc_bz_ctrl_int+dt*Euler_sim(1,i-1)
%     tau_2_ctrl_int=acc_bz_ctrl_int+dt*Euler_sim(1,i-1)
%     tau_3_ctrl_int=acc_bz_ctrl_int+dt*Euler_sim(1,i-1)
    v_sp=(z_sp-p3)*2;
    
    acc_bz_ctrl=(v_sp-v3)*5;
    tau_ctrl(:,i)=(omega_sp(:,i)-[omega1 omega2 omega3]').*2;
    
    %dynamics
    X_gt(:,i)=                                                                                        [p1 + dt*v1
                                                                                                              p2 + dt*v2
                                                                                                              p3 + dt*v3
                                    v1 + dt*(a_s3*(2*a*c + 2*b*d) - a_s2*(2*a*d - 2*b*c) + a_s1*(a^2 + b^2 - c^2 - d^2))
                                    v2 + dt*(a_s1*(2*a*d + 2*b*c) - a_s3*(2*a*b - 2*c*d) + a_s2*(a^2 - b^2 + c^2 - d^2))
 v3 + dt*(a_s2*(2*a*b + 2*c*d) - a_s1*(2*a*c - 2*b*d) + a_s3*(a^2 - b^2 - c^2 + d^2) + 5569248142366845/562949953421312)
                                                                     a - dt*((b*omega1)/2 + (c*omega2)/2 + (d*omega3)/2)
                                                                     b + dt*((a*omega1)/2 + (c*omega3)/2 - (d*omega2)/2)
                                                                     c + dt*((a*omega2)/2 - (b*omega3)/2 + (d*omega1)/2)
                                                                     d + dt*((a*omega3)/2 + (b*omega2)/2 - (c*omega1)/2)
                                                                                                        omega1 + dt*tau1
                                                                                                        omega2 + dt*tau2
                                                                                                        omega3 + dt*tau3
                                                                                                                    0
                                                                                                                    0
                                                                                                                    acc_bz_in+acc_bz_ctrl
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                    tau_1_in+tau_ctrl(1,i)
                                                                                                                    tau_2_in+tau_ctrl(2,i)
                                                                                                                    tau_3_in+tau_ctrl(3,i)
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0
                                                                                                                       0.5
                                                                                                                       0.5
                                                                                                                       0.5
                                                                                                                       g];
    
    [ p1, p2, p3, v1, v2, v3, a, b, c, d, omega1, omega2, omega3, a_s1, a_s2, a_s3, a_s_11, a_s_12, a_s_13, a_s_21, a_s_22, a_s_23, a_s_31, a_s_32, a_s_33, tau1, tau2, tau3, tau_11, tau_12, tau_13, tau_21, tau_22, tau_23, tau_31, tau_32, tau_33, ba1, ba2, ba3, bg1, bg2, bg3, bm1, bm2, bm3, l_ic1, l_ic2, l_ic3, g] = Assignment(X_gt(:,i));

    q_tmp=X_gt([7:1:10],i);
    q_tmp=q_tmp./norm(q_tmp);
    Euler_sim(:,i)=quat2eul(q_tmp');
    X_gt([7:1:10],i)=q_tmp;
    Z_gt(:,i)=[ p1, p2, p3, 0, 0, 0,a,b,c,d, a_s1 - l_ic2*tau3 + l_ic3*tau2 - omega2*(l_ic1*omega2 - l_ic2*omega1) - omega3*(l_ic1*omega3 - l_ic3*omega1), a_s2 + l_ic1*tau3 - l_ic3*tau1 + omega1*(l_ic1*omega2 - l_ic2*omega1) - omega3*(l_ic2*omega3 - l_ic3*omega2), a_s3 - l_ic1*tau2 + l_ic2*tau1 + omega1*(l_ic1*omega3 - l_ic3*omega1) + omega2*(l_ic2*omega3 - l_ic3*omega2), omega1, omega2, omega3, a^2 + b^2 - c^2 - d^2 + bm1, bm2 + 2*a*d + 2*b*c, bm3 - 2*a*c + 2*b*d];
    noise=R*randn(length(R),1);
    Z_t(:,i)=Z_gt(:,i)+noise;
    i=i+1;
end
time_h=[0:dt:T];
time_l=[0.0001:7*dt:T];
acc_gryo_in=[time_h',Z_t([11:16],:)'];
pos_in=[time_l',interp1(time_h,Z_t([1 2 3],:)',time_l,'linear')];
mag_in=[time_l',interp1(time_h,Z_t([17 18 19],:)',time_l,'linear')];
Q_gt=[time_h', X_gt([7 8 9 10],:)'];
Q_t=[time_l', interp1(time_h,Z_t([7 8 9 10],:)',time_l,'linear')];
end
% time=[0:dt:T];
% figure(1)
% plot(time,X_gt([1 2 3],:))
% figure(2)
% plot(time,X_gt([4 5 6],:))
% figure(3)
% plot(time,X_gt([7 8 9 10],:))
% %plot(time,X_gt(16,:))
