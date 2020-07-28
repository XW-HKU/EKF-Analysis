clc;clear
syms a b c d g dt real
global X_def Z_def P Q R;
global COMP_T;
COMP_T=0;

%1 read data
start_time = 59.7;
end_time   = 100;
%[pos_in, acc_gryo_in, Q_gt] = bagprocess('subset_2020-01-03-22-00-28.bag');
acc_gryo_in=csvread('03_07_10_sensor_raw_test_0.csv',1,0);
acc_gryo_in(:,1)=acc_gryo_in(:,1)./1E6 ;
ind1=find((acc_gryo_in(:,1)>=start_time)&(acc_gryo_in(:,1)<=end_time));
acc_gryo_in=acc_gryo_in(ind1,[1 4 5 6 7 8 9]);
acc_gryo_in(:,1)=acc_gryo_in(:,1);

%% normal filter
[b,a] = butter(1,20/800);
freqz(b,a);
acc_butfilt(:,1) = filter(b,a,acc_gryo_in(:,2));
acc_butfilt(:,2) = filter(b,a,acc_gryo_in(:,3));
acc_butfilt(:,3) = filter(b,a,acc_gryo_in(:,4));
%plot(acc_gryo_in(:,2),)


%% spectrum
figure(11)
subplot(1,3,1);
IMU_Fr=length(acc_gryo_in(:,1))/(acc_gryo_in(end,1)-acc_gryo_in(1,1));

[pxx,f] = pwelch(acc_gryo_in(:,2),[],[],[],IMU_Fr);
loglog(f,pxx);
xlim([0.1 100]);
ylim([10^(-5) 1]);
title('IMU measurements');
subplot(1,3,2);
[pxx,f] = pwelch(acc_butfilt(:,1),[],[],[],IMU_Fr);
loglog(f,pxx);
xlim([0.1 100]);
ylim([10^(-6) 1]);
title('IMU measurements by filter');
subplot(1,3,3);
EKF_Fr=length(Cur_T)/(Cur_T(end)-Cur_T(1));
% periodogram(X_out(7,:),rectwin(length(X_out(7,:))),length(X_out(7,:)),EKF_Fr);
[pxx,f] = pwelch(X_out(7,:),[],[],[],EKF_Fr);
loglog(f,pxx);
xlim([0.1 100]);
ylim([10^(-6) 1]);
title('EKF output');