clc;clear
syms a b c d g dt real
global X_def Z_def P Q R;
global COMP_T;
COMP_T=0;

%1 read data
start_time = 59.7;
end_time   = 100;
%[pos_in, acc_gryo_in, Q_gt] = bagprocess('subset_2020-01-03-22-00-28.bag');
acc_gryo_in=csvread([pwd '/03_07_10_sensor_raw_test_0.csv'],1,0);
acc_gryo_in(:,1)=acc_gryo_in(:,1)./1E6 ;
ind1=find((acc_gryo_in(:,1)>=start_time)&(acc_gryo_in(:,1)<=end_time));
acc_gryo_in=acc_gryo_in(ind1,[1 4 5 6 7 8 9]);
acc_gryo_in(:,1)=acc_gryo_in(:,1);

pos_in=csvread([pwd '/03_07_10_vehicle_visual_odometry_0.csv'],1,0);
pos_in(:,1)=pos_in(:,1)./1E6;
ind2=find((pos_in(:,1)>=start_time)&(pos_in(:,1)<=end_time));
Q_ODO=pos_in(ind2,[5 6 7 8]);
pos_full=pos_in(:,[1 2 3 4]);
pos_in=pos_in(ind2,[1 2 3 4]);
pos_in(:,1)=pos_in(:,1);

%%
Q_in=csvread([pwd '/03_07_10_vehicle_attitude_0.csv'],1,0);
Q_in(:,1)=Q_in(:,1)./1E6;
ind3=find((Q_in(:,1)>=start_time)&(Q_in(:,1)<=end_time));
Q_in=Q_in(ind3,[1 5 6 7 8]);
Q_gt(:,1)=pos_in(:,1);
Q_gt(:,[2:5])=interp1(Q_in(:,1),Q_in(:,[2:5]),pos_in(:,1),'linear');
Q_gt(1,[2:5])=Q_gt(2,[2:5]);
for i=1:length(Q_gt)
    mag_in(i,1)=Q_gt(i,1);
    mag_in(i,[2 3 4])=quat2rotm(Q_gt(i,[2:5]))'*[1 0 0]';
    R_in(:,:,i)=quat2rotm(Q_gt(i,[2:5]));
    a=Q_gt(i,2);b=Q_gt(i,3);c=Q_gt(i,4);d=Q_gt(i,5);
    R_in2(:,:,i)=[a*a+b*b-c*c-d*d 2*(b*c-a*d) 2*(a*c+b*d);
    2*(b*c+a*d)     a*a-b*b+c*c-d*d 2*(c*d-a*b);
    2*(b*d-a*c)       2*(a*b+c*d)   a*a-b*b-c*c+d*d];
    Euler_gt(:,i)=[Q_gt(i,1),quat2eul(Q_gt(i,[2:5]))];
end
G_norm=9.78;
pos_start=[pos_in(1,2) pos_in(1,3) pos_in(1,4)]';
Q_start=[Q_gt(1,2) Q_gt(1,3) Q_gt(1,4) Q_gt(1,5)]';

X_NUM=42;
X_out(:,1)=zeros(X_NUM,1);X_out(9,1)=-9.8;
RotMatrix(:,:,1)=quat2rotm(Q_start');
X_prid(:,1)=X_out(:,1);
Cur_T(1)=min(acc_gryo_in(1,1),mag_in(1,1));

%%
i=2;
mag_index=1;
iter_index=2;
Dt(iter_index)=0;
P_in=eye(X_NUM);
Euler(:,1)=quat2eul(Q_start');
Euler_gt(:,1)=[0;quat2eul(Q_start')'];
P_prid=zeros(X_NUM,X_NUM,length(acc_gryo_in)-1);
P_next=zeros(X_NUM,X_NUM,length(acc_gryo_in)-1);
F_now=zeros(X_NUM,X_NUM,length(acc_gryo_in)-1);
P_prid(:,:,1)=eye(X_NUM);
P_next(:,:,1)=eye(X_NUM);
F_now(:,:,1) =eye(X_NUM);

while i<length(acc_gryo_in)-1
    if mag_index<=1
        Z_in(:,iter_index)=[pos_start;[acc_gryo_in(i-1,2) acc_gryo_in(i-1,3) acc_gryo_in(i-1,4)]';[acc_gryo_in(i-1,5) acc_gryo_in(i-1,6) acc_gryo_in(i-1,7)]';[mag_in(1,2) mag_in(1,3) mag_in(1,4)]'];
    else
        Z_in(:,iter_index)=[[pos_in(mag_index,2) pos_in(mag_index,3) pos_in(mag_index,4)]';[acc_gryo_in(i-1,2) acc_gryo_in(i-1,3) acc_gryo_in(i-1,4)]';[acc_gryo_in(i-1,5) acc_gryo_in(i-1,6) acc_gryo_in(i-1,7)]';[mag_in(mag_index,2) mag_in(mag_index,3) mag_in(mag_index,4)]'];
    end
    
    if(mag_index==length(mag_in))
        Cur_T(iter_index)=acc_gryo_in(i,1);
        flag='acc_gy';
        i=i+1;
    elseif(acc_gryo_in(i,1)<=mag_in(mag_index,1))
        Cur_T(iter_index)=acc_gryo_in(i,1);
        flag='acc_gy';
        i=i+1;
    else
        flag='mag';
        Cur_T(iter_index)=mag_in(mag_index,1);
        mag_index=mag_index+1;
    end
    flag='None';
    
    Dt=Cur_T(iter_index)-Cur_T(iter_index-1);
    
    Use_states=1;
    if (Use_states==1)
        [X_prid(:,iter_index),X_out(:,iter_index),RotMatrix(:,:,iter_index),P_prid(:,:,iter_index),P_next(:,:,iter_index),F_now(:,:,iter_index)]=System_Mode(X_out(:,iter_index-1),RotMatrix(:,:,iter_index-1),Z_in(:,iter_index),P_in,Dt,flag);
        name='Filtered';
    else
        [X_prid(:,iter_index),X_out(:,iter_index),RotMatrix(:,:,iter_index),P_prid(:,:,iter_index),P_next(:,:,iter_index),F_now(:,:,iter_index)]=System_Mode_all_inputs(X_out(:,iter_index-1),RotMatrix(:,:,iter_index-1),Z_in(:,iter_index),P_in,Dt,flag);
        name='Filtered';
    end
    P_in=P_next(:,:,iter_index);
    Euler(:,iter_index)=rotm2eul(RotMatrix(:,:,iter_index));
    Euler_gt(:,mag_index)=[Cur_T(iter_index);quat2eul(Q_gt(mag_index,[2 3 4 5]))'];
    m_p_f(:,iter_index)=X_out(1:3,iter_index)+RotMatrix(:,:,iter_index)*X_out(end-5:end-3,iter_index);
    iter_index=iter_index+1;
end
%% Inputs
i=2;
mag_index=1;
iter_index=2;
Dt(iter_index)=0;
X_out_I(:,1)=zeros(X_NUM,1);X_out_I(9,1)=-9.8;X_out([1:3],1)=pos_start;
RotMatrix_I(:,:,1)=eye(3);
P_in=eye(X_NUM);
Euler(:,1)=quat2eul(Q_start');
P_prid(:,:,1)=eye(X_NUM);
P_next(:,:,1)=eye(X_NUM);
F_now(:,:,1) =eye(X_NUM);

while i<length(acc_gryo_in)-1
    if(mag_index==length(mag_in))
        Cur_T(iter_index)=acc_gryo_in(i,1);
        flag='acc_gy';
        i=i+1;
    elseif(acc_gryo_in(i,1)<=mag_in(mag_index,1))
        Cur_T(iter_index)=acc_gryo_in(i,1);
        flag='acc_gy';
        i=i+1;
    else
        flag='mag';
        Cur_T(iter_index)=mag_in(mag_index,1);
        mag_index=mag_index+1;
    end
    
    Dt=Cur_T(iter_index)-Cur_T(iter_index-1);
    'input';
    [X_prid(:,iter_index),X_out_I(:,iter_index),RotMatrix_I(:,:,iter_index),P_prid(:,:,iter_index),P_next(:,:,iter_index),F_now(:,:,iter_index)]=System_Mode_all_inputs(X_out_I(:,iter_index-1),RotMatrix_I(:,:,iter_index-1),Z_in(:,iter_index),P_in,Dt,flag);
    P_in=P_next(:,:,iter_index);
    Euler_I(:,iter_index)=rotm2eul(RotMatrix_I(:,:,iter_index));
    %Euler_gt(:,mag_index)=[Cur_T(iter_index);quat2eul(Q_gt(mag_index,[2 3 4 5]))'];
    m_p_f_I(:,iter_index)=X_out_I(1:3,iter_index)+RotMatrix_I(:,:,iter_index)*X_out_I(end-5:end-3,iter_index);
    iter_index=iter_index+1;
end
%%
latex_font=11
figure(111)
[b,a] = butter(2,20/850);
d1 = designfilt('lowpassiir','FilterOrder',1, ...
    'HalfPowerFrequency',40/850,'DesignMethod','butter');
%freqz(b,a);
IMU_raw=Z_in(4:9,:);
imu_butfilt(1,:) = filter(b,a,IMU_raw(1,:));
imu_butfilt(2,:) = filter(b,a,IMU_raw(2,:));
imu_butfilt(3,:) = filter(b,a,IMU_raw(3,:));

imu_zerofilt(1,:)=filtfilt(d1,IMU_raw(1,:));
imu_zerofilt(2,:)=filtfilt(d1,IMU_raw(2,:));
imu_zerofilt(3,:)=filtfilt(d1,IMU_raw(3,:));

[b,a] = butter(2,40/850);
d1 = designfilt('lowpassiir','FilterOrder',1, ...
    'HalfPowerFrequency',60/850,'DesignMethod','butter');
imu_butfilt(4,:) = filter(b,a,IMU_raw(4,:));
imu_butfilt(5,:) = filter(b,a,IMU_raw(5,:));
imu_butfilt(6,:) = filter(b,a,IMU_raw(6,:));
imu_zerofilt(4,:)=filtfilt(d1,IMU_raw(4,:));
imu_zerofilt(5,:)=filtfilt(d1,IMU_raw(5,:));
imu_zerofilt(6,:)=filtfilt(d1,IMU_raw(6,:));

% Compare filter
latex_font=12
figure(11)
set(gcf,'Position',[600,200,900,500]);
T=tiledlayout(3,2);
T.Padding = 'none';
T.TileSpacing = 'none';
t_start=21.5;
t_end=22.3;
% Body X
nexttile;
plot(Cur_T,IMU_raw(1,:),'color',[96 96 96]/150,'linewidth',0.3);hold on;
plot(Cur_T,imu_butfilt(1,:),'r:','linewidth',1.5);hold on;
plot(Cur_T,imu_zerofilt(1,:),'g-.','linewidth',1.5);hold on;
plot(Cur_T,X_out(7,:),'b:','linewidth',1.5);
xlim([t_start,t_end]);
ylim([-5 15]);
title('Acceleration','interpreter','latex','FontSize', latex_font);
ylabel('Body-X $(m/s^2)$','interpreter','latex','FontSize', latex_font);
legend('Measurement','Normal Low-Pass','Zero Phase Low-Pss','Our Method','Location','Northwest','interpreter','latex','FontSize', latex_font-2);

nexttile;
plot(Cur_T,IMU_raw(4,:),'color',[96 96 96]/150,'linewidth',0.3);hold on;
plot(Cur_T,imu_butfilt(4,:),'r:','linewidth',1.5);hold on;
plot(Cur_T,imu_zerofilt(4,:),'g-.','linewidth',1.5);hold on;
plot(Cur_T,X_out(19,:),'b:','linewidth',1.5);hold on;
xlim([t_start,t_end]);
%ylim([-10 30]);
title('Angular Velocity','interpreter','latex','FontSize', latex_font);
ylabel('Body-X $(rad/s)$','interpreter','latex','FontSize', latex_font);

% Body Y
nexttile;
plot(Cur_T,IMU_raw(2,:),'color',[96 96 96]/150,'linewidth',0.3);hold on;
plot(Cur_T,imu_butfilt(2,:),'r:','linewidth',1.5);hold on;
plot(Cur_T,imu_zerofilt(2,:),'g-.','linewidth',1.5);hold on;
plot(Cur_T,X_out(8,:),'b:','linewidth',1.5);hold on;
xlim([t_start,t_end]);
ylim([-8 8]);
ylabel('Body-Y $ (m/s^2)$','interpreter','latex','FontSize', latex_font);

nexttile;
plot(Cur_T,IMU_raw(5,:),'color',[96 96 96]/150,'linewidth',0.3);hold on;
plot(Cur_T,imu_butfilt(5,:),'r:','linewidth',1.5);hold on;
plot(Cur_T,imu_zerofilt(5,:),'g-.','linewidth',1.5);hold on;
plot(Cur_T,X_out(20,:),'b:','linewidth',1.5);hold on;
xlim([t_start,t_end]);
ylim([-1 3]);
ylabel('Body-Y $ (rad/s)$','interpreter','latex','FontSize', latex_font);
legend('Measurement','Normal Low-Pass','Zero Phase Low-Pss','Our Method','Location','Northeast','interpreter','latex','FontSize', latex_font-2);

% Body Z
nexttile;
plot(Cur_T,IMU_raw(3,:),'color',[96 96 96]/150,'linewidth',0.3);hold on;
plot(Cur_T,imu_butfilt(3,:),'r:','linewidth',1.5);hold on;
plot(Cur_T,imu_zerofilt(3,:),'g-.','linewidth',1.5);hold on;
plot(Cur_T,X_out(9,:),'b:','linewidth',1.5);hold on;
xlim([t_start,t_end]);
xlabel('Time $(s)$','interpreter','latex','FontSize', latex_font);
ylabel('Body-Z $ (m/s^2)$','interpreter','latex','FontSize', latex_font);

nexttile;
plot(Cur_T,IMU_raw(6,:),'color',[96 96 96]/150,'linewidth',0.3);hold on;
plot(Cur_T,imu_butfilt(6,:),'r:','linewidth',1.5);hold on;
plot(Cur_T,imu_zerofilt(6,:),'g-.','linewidth',1.5);hold on;
plot(Cur_T,X_out(21,:),'b:','linewidth',1.5);hold on;
xlim([t_start,t_end]);
xlabel('Time $(s)$','interpreter','latex','FontSize', latex_font);
ylabel('Body-Z $ (rad/s)$','interpreter','latex','FontSize', latex_font);

% smoother
Euler_ODO=(interp1(Euler_gt(1,:)',Euler_gt([2 3 4],:)',Cur_T','linear'))';
%%
'step4: figure'
latex_font=12;
index=[round(0.7*length(Cur_T)):round(1*length(Cur_T))];

Error_pos_Input=m_p_f_I(1:3,index)-Z_in(1:3,index);
Error_atti_Input=(Euler_I(1:3,index)-Euler_ODO(1:3,index))*0.5;
Err_c_Input=X_out_I(37:39,index)-[-0.02 -0.01 -0.070]'*ones(1,length(index));
Error_Input_UAV=[Error_pos_Input*1000;Error_atti_Input*57.3;Err_c_Input*1000];
ErrMean_Input= mean(Error_Input_UAV,2);
ErrMid_Input=median(Error_Input_UAV,2);
ErrSec_Input= [min(Error_Input_UAV,[],2),max(Error_Input_UAV,[],2)];

Error_pos_our=m_p_f(1:3,index)-Z_in(1:3,index);
Error_atti_our=(Euler(1:3,index)-Euler_ODO(1:3,index))*0.5;
Err_c_our=X_out(37:39,index)-[-0.02 -0.01 -0.07]'*ones(1,length(index));
Error_our_UAV=[Error_pos_our*1000;Error_atti_our*57.3;Err_c_our*1000];
ErrMean_our=mean(Error_our_UAV,2);
ErrMid_our=median(Error_our_UAV,2);
ErrSec_our= [min(Error_our_UAV,[],2),max(Error_our_UAV,[],2)]%-ErrMean_our;

N_err=length(ErrMean_Input);

figure(1)
plot3(Z_in(1,:),Z_in(2,:),-Z_in(3,:),'k.')
hold on; plot3(m_p_f(1,:),m_p_f(2,:),-m_p_f(3,:),'b.');
hold on; plot3(m_p_f_I(1,:),m_p_f_I(2,:),-m_p_f_I(3,:),'r.');
legend('Measurements','Our Method','Normal EKF','interpreter','latex','FontSize', latex_font-2);

figure(2)
set(gcf,'Position',[100,100,700,300]);
T=tiledlayout(3,3);
T.Padding = 'compact';
%T.TileSpacing = 'compact';
%ind_i=[1 4 7 2 5 8 3 6 9];
ind_i=1:9;
i=1;
ylabel_array=["$\textbf{p}_{x}$ $(mm)$" "$\textbf{p}_{y}$ $(mm)$" "$\textbf{p}_{z}$ $(mm)$" "Yaw $(deg)$" "Pitch $(deg)$" "Roll $(deg)$" "$\textbf{c}_{x}$ $(mm)$" "$\bf{c}_{y}$ $(mm)$" "$\bf{c}_{z}$ $(mm)$"];
title_array=["Postion $(mm)$" "Attitude $(deg)$" "Translational Offset $(mm)$"];
while i<=N_err
    nexttile;
    k=ind_i(i);
    boxplot([Error_our_UAV(ind_i(i),:)' Error_Input_UAV(ind_i(i),:)'],'symbol','');hold on;
    set(gca, 'XTickLabel',[]);
    set(gca, 'TickLabelInterpreter','latex');
    set(gca, 'FontSize', latex_font-1);
    if(i>=7)
        %title(title_array(i),'interpreter','latex','FontSize', latex_font-1);
        set(gca, 'XTickLabel',{'Our Method','Normal EKF'});
        set(gca, 'TickLabelInterpreter','latex');
        set(gca, 'FontSize', latex_font-1);
    end
    ylabel(ylabel_array(ind_i(i)),'Units', 'Normalized', 'Position', [-0.15, 0.5, 0],'interpreter','latex','FontSize', latex_font);
    i=i+1;
    grid on;
end

save('Error_our_UAV.mat','Error_our_UAV')
save('Error_Input_UAV.mat','Error_Input_UAV')

grid on;
%% p
f1=figure(71)
set(gcf,'Position',[100,100,900,900]);
T=tiledlayout(6,3);
T.Padding = 'none';
T.TileSpacing = 'none';

nexttile([1 2])
plot(Cur_T,Z_in(1,:),'k.','linewidth',1);hold on;
plot(Cur_T,m_p_f(1,:),'b.','linewidth',1);hold on;
plot(Cur_T,m_p_f_I(1,:),'r.','linewidth',1);hold on;
ylabel('Position X($m$)','interpreter','latex','FontSize', latex_font);
legend('From Mocap','Our Method','Normal EKF','interpreter','latex','FontSize',latex_font-2);

nexttile
[errorNum,xaxis] = hist(m_p_f(1,index)-Z_in(1,index),[-0.03:0.001:0.03]);
h1=bar(xaxis,errorNum);h1.FaceColor='none';h1.EdgeColor=[0 0 1];hold on;
[errorNum,xaxis] = hist(m_p_f_I(1,index)-Z_in(1,index),[-0.03:0.001:0.03]); 
h2=bar(xaxis,errorNum);h2.FaceColor='none';h2.EdgeColor=[1 0 0];hold on;
legend('Our Method','Normal EKF','interpreter','latex','FontSize',latex_font-2);
xlim([-0.02,0.02]);
ylim([0,8000]);

nexttile([1 2])
plot(Cur_T,Z_in(2,:),'k.','linewidth',1);hold on;
plot(Cur_T,m_p_f(2,:),'b.','linewidth',1);hold on;
plot(Cur_T,m_p_f_I(2,:),'r.','linewidth',1);hold on;
ylabel('Position Y($m$)','interpreter','latex','FontSize', latex_font);
ylim([-2,1]);
legend('From Mocap','Our Method','Normal EKF','interpreter','latex','FontSize',latex_font-2);

nexttile
[errorNum,xaxis] = hist(m_p_f(2,index)-Z_in(2,index),[-0.03:0.001:0.03]); 
h1=bar(xaxis,errorNum);h1.FaceColor='none';h1.EdgeColor=[0 0 1];hold on;
[errorNum,xaxis] = hist(m_p_f_I(2,index)-Z_in(2,index),[-0.03:0.001:0.03]); 
h2=bar(xaxis,errorNum);h2.FaceColor='none';h2.EdgeColor=[1 0 0];hold on;
legend('Our Method','Normal EKF','interpreter','latex','FontSize',latex_font-2);
xlim([-0.02,0.02]);
ylim([0,6000]);

nexttile([1 2])
plot(Cur_T,-Z_in(3,:),'k','linewidth',1);hold on;
plot(Cur_T,-m_p_f(3,:),'b--','linewidth',1);hold on;
plot(Cur_T,-m_p_f_I(3,:),'r-.','linewidth',1);hold on;
ylabel('Position Z($m$)','interpreter','latex','FontSize', latex_font);
ylim([0,2.5]);
legend('From Mocap','Our Method','Normal EKF','interpreter','latex','FontSize',latex_font-2);

nexttile
[errorNum,xaxis] = hist(m_p_f(3,index)-Z_in(3,index),[-0.03:0.001:0.03]); 
h1=bar(xaxis,errorNum);h1.FaceColor='none';h1.EdgeColor=[0 0 1];hold on;
[errorNum,xaxis] = hist(m_p_f_I(3,index)-Z_in(3,index),[-0.03:0.001:0.03]); 
h2=bar(xaxis,errorNum);h2.FaceColor='none';h2.EdgeColor=[1 0 0];hold on;
legend('Our Method','Normal EKF','interpreter','latex','FontSize',latex_font-2);
xlim([-0.01,0.01])
ylim([0,8000]);

nexttile([1 2])
plot(Cur_T,Euler_ODO(1,:),'k','linewidth',1);hold on;
plot(Cur_T,Euler(1,:),'b--','linewidth',1);hold on;
plot(Cur_T,Euler_I(1,:),'r-.','linewidth',1);hold on;
%title('Attitude Estimation Results','interpreter','latex','FontSize', latex_font)
ylabel('Yaw($rad$)','interpreter','latex','FontSize', latex_font);
%ylim([-0.25,0.25]);
legend('From Mocap','Our Method','Normal EKF','interpreter','latex','FontSize',latex_font-2);

nexttile
[errorNum,xaxis] = hist(Euler(1,index)-Euler_ODO(1,index),[-0.05:0.001:0.05]); 
h1=bar(xaxis,errorNum);h1.FaceColor='none';h1.EdgeColor=[0 0 1];hold on;
[errorNum,xaxis] = hist(Euler_I(1,index)-Euler_ODO(1,index),[-0.05:0.001:0.05]); 
h2=bar(xaxis,errorNum);h2.FaceColor='none';h2.EdgeColor=[1 0 0];hold on;
legend('Our Method','Normal EKF','interpreter','latex','FontSize',latex_font-2);
xlim([-0.02,0.02])
ylim([0,4500]);

nexttile([1 2])
% [errorNum,xaxis] = hist(m_p_f(2,index)-Z_in(2,index),[-0.03:0.001:0.03]); 
% bar(xaxis,errorNum);hold on;
plot(Cur_T,Euler_ODO(2,:),'k','linewidth',1);hold on;
plot(Cur_T,Euler(2,:),'b--','linewidth',1);hold on;
plot(Cur_T,Euler_I(2,:),'r-.','linewidth',1);hold on;
ylabel('Pitch($rad$)','interpreter','latex','FontSize', latex_font);
%ylim([-0.25,0.25]);
legend('From Mocap','Our Method','Normal EKF','interpreter','latex','FontSize',latex_font-2);

nexttile
[errorNum,xaxis] = hist(Euler(2,index)-Euler_ODO(2,index),[-0.3:0.001:0.3]); 
h1=bar(xaxis,errorNum);h1.FaceColor='none';h1.EdgeColor=[0 0 1];hold on;
[errorNum,xaxis] = hist(Euler_I(2,index)-Euler_ODO(2,index),[-0.25:0.001:0.25]); 
h2=bar(xaxis,errorNum);h2.FaceColor='none';h2.EdgeColor=[1 0 0];hold on;
legend('Our Method','Normal EKF','interpreter','latex','FontSize',latex_font-2);
xlim([-0.03,0.03]);
ylim([0,4000]);

nexttile([1 2])
plot(Cur_T,Euler_ODO(3,:),'k','linewidth',1);hold on;
plot(Cur_T,Euler(3,:),'b--','linewidth',1);hold on;
plot(Cur_T,Euler_I(3,:),'r-.','linewidth',1);hold on;
xlabel('Time (s)','interpreter','latex','FontSize', latex_font);
ylabel('Roll($rad$)','interpreter','latex','FontSize', latex_font);
%ylim([-0.25,0.25]);
legend('From Mocap','Our Method','Normal EKF','interpreter','latex','FontSize',latex_font-2);

nexttile
[errorNum,xaxis] = hist(Euler(3,index)-Euler_ODO(3,index),[-0.2:0.001:0.1]); 
h1=bar(xaxis,errorNum);h1.FaceColor='none';h1.EdgeColor=[0 0 1];hold on;
[errorNum,xaxis] = hist(Euler_I(3,index)-Euler_ODO(3,index),[-0.1:0.001:0.1]); 
h2=bar(xaxis,errorNum);h2.FaceColor='none';h2.EdgeColor=[1 0 0];hold on;
legend('Our Method','Normal EKF','interpreter','latex','FontSize',latex_font-2);
xlim([-0.02,0.08]);
ylim([0,3000]);
xlabel('Error distribution after converges','interpreter','latex','FontSize', latex_font);

%% L-IC
f2=figure(4);
latex_font=12;
set(gcf,'Position',[600,600,600,400]);

Cur_T=Cur_T-Cur_T(1);

T=tiledlayout(3,1);
T.Padding = 'none';
T.TileSpacing = 'none';

nexttile;
plot(Cur_T,Cur_T.*0+1,'k','linewidth',1);hold on;
plot(Cur_T,-X_out(37,:)*100+1,'b--','linewidth',1);hold on;
plot(Cur_T,-X_out_I(37,:)*100,'r-.','linewidth',1);hold on;
title('Translational Offset','interpreter','latex','FontSize', latex_font);
ylabel('Body-X (cm)','interpreter','latex','FontSize', latex_font);
legend('From Mocap','Our Method','Normal EKF','interpreter','latex','FontSize',latex_font-2);
%ylim([-5,5]);

nexttile;
plot(Cur_T,Cur_T.*0,'k','linewidth',1);hold on;
plot(Cur_T,X_out(38,:)*100,'b--','linewidth',1);hold on;
plot(Cur_T,X_out_I(38,:)*100,'r-.','linewidth',1);hold on;
ylabel('Body-Y (cm)','interpreter','latex','FontSize', latex_font);
%ylim([-5,5]);

nexttile;
plot(Cur_T,Cur_T.*0-7.0,'k','linewidth',1);hold on;
plot(Cur_T,X_out(39,:)*100,'b--','linewidth',1);hold on;
plot(Cur_T,X_out_I(39,:)*100,'r-.','linewidth',1);hold on;
%ylim([-0.2,0.1]);
ylabel('Body-Z (cm)','interpreter','latex','FontSize', latex_font);
xlabel('Time (s)','interpreter','latex','FontSize', latex_font);
%% plot
% figure(22)
% % omega
% subplot(3,1,1)
% plot(Cur_T,Z_in(7,:));hold on;plot(Cur_T,X_out(19,:));legend('Measured',name );title('Angular Velocity X','interpreter','latex','FontSize', latex_font);
% ylabel('$rad/s$','interpreter','latex','FontSize', latex_font);
% set(gca,'LooseInset',get(gca,'TightInset'));
% subplot(3,1,2)
% plot(Cur_T,Z_in(8,:));hold on;plot(Cur_T,X_out(20,:));legend('Measured',name );title('Angular Velocity Y','interpreter','latex','FontSize', latex_font);
% ylabel('$rad/s$','interpreter','latex','FontSize', latex_font);
% set(gca,'LooseInset',get(gca,'TightInset'));
% subplot(3,1,3)
% plot(Cur_T,Z_in(9,:));hold on;plot(Cur_T,X_out(21,:));legend('Measured',name );title('Angular Velocity Z','interpreter','latex','FontSize', latex_font);
% xlabel('time (s)','interpreter','latex','FontSize', latex_font); ylabel('$rad/s$','interpreter','latex','FontSize', latex_font);
% set(gca,'LooseInset',get(gca,'TightInset'));
% 
% figure(33)
% % a_s
% subplot(3,1,1)
% plot(Cur_T,Z_in(4,:));hold on;plot(Cur_T,X_out(7,:));legend('Measured',name );title('Acceleration X','interpreter','latex','FontSize', latex_font);
% ylabel('$m/s^2$','interpreter','latex','FontSize', latex_font);
% set(gca,'LooseInset',get(gca,'TightInset'));
% subplot(3,1,2)
% plot(Cur_T,Z_in(5,:));hold on;plot(Cur_T,X_out(8,:));legend('Measured',name );title('Acceleration Y','interpreter','latex','FontSize', latex_font);
% ylabel('$m/s^2$','interpreter','latex','FontSize', latex_font);
% set(gca,'LooseInset',get(gca,'TightInset'));
% subplot(3,1,3)
% plot(Cur_T,Z_in(6,:));hold on;plot(Cur_T,X_out(9,:));legend('Measured',name );title('Acceleration Z','interpreter','latex','FontSize', latex_font);
% xlabel('time (s)','interpreter','latex','FontSize', latex_font); ylabel('$m/s^2$','interpreter','latex','FontSize', latex_font);
% set(gca,'LooseInset',get(gca,'TightInset'));