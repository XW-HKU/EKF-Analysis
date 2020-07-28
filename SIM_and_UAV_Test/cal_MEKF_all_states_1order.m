clear;clc
syms p1 p2 p3 v1 v2 v3 a b c d a_s1 a_s2 a_s3 omega1 omega2 omega3 a_s1 a_s2 a_s3 a_s_11 a_s_12 a_s_13 a_s_21 a_s_22 a_s_23 a_s_31 a_s_32 a_s_33 tau1 tau2 tau3 ba1 ba2 ba3 bg1 bg2 bg3 l_ic1 l_ic2 l_ic3 g dt real;
syms ptheta1 ptheta2 ptheta3 real;
%1-calculate model
p=[p1 p2 p3]';
v=[v1 v2 v3]';
q=[a b c d]';
a_s=[a_s1 a_s2 a_s3]';
a_s_1=[a_s_11 a_s_12 a_s_13]';
a_s_2=[a_s_21 a_s_22 a_s_23]';
a_s_3=[a_s_31 a_s_32 a_s_33]';
omega=[omega1 omega2 omega3]';
tau=[tau1 tau2 tau3]';
tau_1=sym('tau_1',[3,1]);
tau_2=sym('tau_2',[3,1]);
tau_3=sym('tau_3',[3,1]);
ba=[ba1 ba2 ba3]';
bg=[bg1 bg2 bg3]';
bm=sym('bm',[3,1]);
l_ic=[l_ic1 l_ic2 l_ic3]';
R=sym('R',[3,3]);
ptheta=[ptheta1 ptheta2 ptheta3]';

syms a_s_in1 a_s_in2 a_s_in3 omega_in1 omega_in2 omega_in3 real;
a_s_in=[a_s_in1 a_s_in2 a_s_in3]';
omega_in=[omega_in1 omega_in2 omega_in3]';

dp=v;
dv=R*a_s+[0;0;9.8];
dq=0.5*[-b -c -d;a -d c;d a -b;-c b a]*(omega);
domega=zeros(3,1);%tau;
da_s=zeros(3,1);%a_s_1;
da_s_1=zeros(3,1);%a_s_2;
da_s_2=zeros(3,1);%a_s_3;
da_s_3=zeros(3,1);
dtau=zeros(3,1);%tau_1;
dtau_1=zeros(3,1);%tau_2;
dtau_2=zeros(3,1);%tau_3;
dtau_3=zeros(3,1);
dba=zeros(3,1);
dbg=zeros(3,1);
dbm=zeros(3,1);
dl_ic=zeros(3,1);
dg=0;

syms X real;

%% distance estimation
X=[p;v;a_s;a_s_1;a_s_2;a_s_3;omega;ba;bg;l_ic;ptheta];
dX=[dp;dv;zeros(3,1);zeros(3,1);zeros(3,1);zeros(3,1);zeros(3,1);zeros(3,1);zeros(3,1);zeros(3,1);omega]

X_pri = X + dt * dX;
%X_pri(end-2:end)=0

for i=1:1:length(X)
    F(:,i)=diff(X_pri, X(i))';
end

%F(end-2:end,end-2:end)=expm(cross_mat(-(omega_in-bg)*dt));
F(4:6,end-2:end)=-dt*R*cross_mat(a_s);
%F(4:6,22:24)=-dt*R

% Measurement Function
m_p=p+R*l_ic;
m_v=zeros(3,1);
m_q=q;
m_a=a_s+ba;
m_omega=omega+bg;
m_m=R*[1;0;0];

syms Z real;
Z=[m_p;m_a;m_omega;m_m];

for i=1:1:length(X)
    H(:,i)=diff(Z, X(i))';
end
H(1:3,end-2:end)=-R*cross_mat(l_ic);%[cross(l_ic,[1 0 0]'),cross(l_ic,[0 1 0]'),cross(l_ic,[0 0 1]')];
%H(4:6,end-2:end)= cross_mat(R'*a_s);%[cross(R'*a_s,[1 0 0]'),cross(R'*a_s,[0 1 0]'),cross(R'*a_s,[0 0 1]')];
H(end-2:end,end-2:end)=-R*cross_mat([1;0;0]);%[cross([1;0;0],[1 0 0]'),cross([1;0;0],[0 1 0]'),cross([1;0;0],[0 0 1]')]

X_pri
Z