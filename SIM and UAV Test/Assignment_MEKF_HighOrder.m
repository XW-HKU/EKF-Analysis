function [ p1, p2, p3, v1, v2, v3, a_s1, a_s2, a_s3, a_s_11, a_s_12, a_s_13, a_s_21, a_s_22, a_s_23, a_s_31, a_s_32, a_s_33, omega1, omega2, omega3, tau1, tau2, tau3, tau_11, tau_12, tau_13, tau_21, tau_22, tau_23, ba1, ba2, ba3, bg1, bg2, bg3, l_ic1, l_ic2, l_ic3] = Assignment_MEKF_HighOrder(X_in)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% X=[p;v;a_s;a_s_1;a_s_2;a_s_3;omega;ba;bg;l_ic;ptheta];

p1=X_in(1);p2=X_in(2);p3=X_in(3);
v1=X_in(4);v2=X_in(5);v3=X_in(6);
a_s1=X_in(7);a_s2=X_in(8);a_s3=X_in(9);
a_s_11=X_in(10);a_s_12=X_in(11);a_s_13=X_in(12);
a_s_21=X_in(13);a_s_22=X_in(14);a_s_23=X_in(15);
a_s_31=X_in(16);a_s_32=X_in(17);a_s_33=X_in(18);
omega1=X_in(19);omega2=X_in(20);omega3=X_in(21);
tau1=X_in(22);tau2=X_in(23);tau3=X_in(24);
tau_11=X_in(25);tau_12=X_in(26);tau_13=X_in(27);
tau_21=X_in(28);tau_22=X_in(29);tau_23=X_in(30);
ba1=X_in(31);ba2=X_in(32);ba3=X_in(33);
bg1=X_in(34);bg2=X_in(35);bg3=X_in(36);
l_ic1=X_in(37);l_ic2=X_in(38);l_ic3=X_in(39);

%omega1=X_smooth(13);omega2=X_smooth(14);omega3=X_smooth(15);
%tau1=X_smooth(16);tau2=X_smooth(17);tau3=X_smooth(18);

 %omega1=X_gt(1);omega2=X_gt(2);omega3=X_gt(3);
% tau1=X_gt(26);tau2=X_gt(27);tau3=X_gt(28);
end
