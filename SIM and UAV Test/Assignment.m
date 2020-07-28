function [ p1, p2, p3, v1, v2, v3, a, b, c, d, omega1, omega2, omega3, a_s1, a_s2, a_s3, a_s_11, a_s_12, a_s_13, a_s_21, a_s_22, a_s_23, a_s_31, a_s_32, a_s_33, tau1, tau2, tau3, tau_11, tau_12, tau_13, tau_21, tau_22, tau_23, tau_31, tau_32, tau_33, ba1, ba2, ba3, bg1, bg2, bg3, bm1, bm2, bm3, l_ic1, l_ic2, l_ic3, g] = Assignment(X_in)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
p1=X_in(1);p2=X_in(2);p3=X_in(3);
v1=X_in(4);v2=X_in(5);v3=X_in(6);
a=X_in(7);b=X_in(8);c=X_in(9);d=X_in(10);
omega1=X_in(11);omega2=X_in(12);omega3=X_in(13);
a_s1=X_in(14);a_s2=X_in(15);a_s3=X_in(16);
a_s_11=X_in(17);a_s_12=X_in(18);a_s_13=X_in(19);
a_s_21=X_in(20);a_s_22=X_in(21);a_s_23=X_in(22);
a_s_31=X_in(23);a_s_32=X_in(24);a_s_33=X_in(25);
tau1=X_in(26);tau2=X_in(27);tau3=X_in(28);
tau_11=X_in(29);tau_12=X_in(30);tau_13=X_in(31);
tau_21=X_in(32);tau_22=X_in(33);tau_23=X_in(34);
tau_31=X_in(35);tau_32=X_in(36);tau_33=X_in(37);
ba1=X_in(38);ba2=X_in(39);ba3=X_in(40);
bg1=X_in(41);bg2=X_in(42);bg3=X_in(43);
bm1=X_in(44);bm2=X_in(45);bm3=X_in(46);
l_ic1=X_in(47);l_ic2=X_in(48);l_ic3=X_in(49);
g=X_in(50);
end

