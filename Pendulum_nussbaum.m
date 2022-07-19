
clc
clear all
close all
%%  %%Adaptive Terminal Sliding Mode Control a Howitzer Shell Transfer Arm with Friction Compensation  ( ANFTSM)  %%
n=1;t=0;
Dt=0.0001;

x1=0.0;
x2=0;x3=0;
e0=0;e11=0;E12=0;
b=50;
z=00;iae=0;itae=0;ce=0;
uu=0;u_1=0;a1=0;a2=0;a3=0;
ydr1=0;ydr2=0;

k1=80;k2=0.1;k=100;
r=3;m=5/3;


dx1=0.0;ei=0;e12=0;
delta=0;E12=0;
for i=1:200000
%     d=0.2*sin(1*t)*exp(-0.01*t)+1.0*sin(1*t)*cos(1*t);  %正弦误差
    d=0.1*sin(1*t);  %正弦误差
    f=(9.81*sin(x1)-0.492*0.21*0.65*x2^2*sin(2*x1))/(0.65*2/3-0.492*0.21*0.65/2*(cos(x1))^2);
    g=(-0.492*1*cos(x1))/(0.65*2/3-0.492*0.21*0.65/2*(cos(x1))^2);

    yd=1*cos(6*0.2*t)*cos(0.2*t); %输入指令
    dyd=-1*6*0.2*sin(6*0.2*t)*cos(0.2*t)-1*0.2*cos(6*0.2*t)*sin(0.2*t); %输入信号的一阶导数
    ddyd=-1*6*0.2*6*0.2*cos(6*0.2*t)*cos(0.2*t)+1*6*0.2*0.2*sin(6*0.2*t)*sin(0.2*t)-...
    1*6*0.2*0.2*sin(6*0.2*t)*sin(0.2*t)+1*0.2*0.2*cos(6*0.2*t)*cos(0.2*t); %输入信号的二阶导数
   
    y=x1; %实际输出
    dy=x2; %实际输出的一阶导数
    e1=y-yd;%误差
    e2=dy-dyd;
    
    c1=0.8;c2=0.1;
    r=2;p14=0.05477;
    B=3.001/3;
m=5/3.0;
k1=0.1;e11=0;
k2=0.1;

   s1=e1+k1*(abs(e1))^3*sign(e1)+k2*(abs(e2))^m*sign(e2)+0.8*(tanh(e1/0.08))*0;%滑模变量
    s=s1;
    uu=(1/g*(-500*s-1*(a1+11.520)*sign(s)+...
    ddyd-1/k2/m*(abs(e2*(1+k1*3*(abs(e1))^(3-1))))^(2-m)*1*sign(e2*(1+k1*3*(abs(e1))^(3-1))))); %u
    F2=k2*B*(abs(e2))^(m-1);
    a1=a1+10.05*1*abs(s)*F2*Dt;

%%%倒立摆模型
    f=(9.81*sin(x1)-0.492*0.21*0.65*x2^2*sin(2*x1))/(0.65*2/3-0.492*0.21*0.65/2*(cos(x1))^2);
    g=(-0.492*1*cos(x1))/(0.65*2/3-0.492*0.21*0.65/2*(cos(x1))^2);
   
    Dx1=x2;
    Dx2=f+g*uu+1*d;

    if n>=4880
        iae=iae+abs(e1)*Dt;
        itae=itae+t*abs(e1)*Dt;
        ce=ce+uu^2*Dt;
    end
    
    x1=x1+Dx1*Dt;
    x2=x2+Dx2*Dt;
    E12=E12+e1^2;      
    e12=e12+e1^2*Dt;
    yr_store(:,n)=[yd;ydr1];
  
    yy1_store(:,n)=[yd;y];
    uu1_store(n)=uu;s1_store(n)=s;
    ee1_store(n)=e1;
    u_1=uu;
    n=n+1;
    t=t+Dt;
    
end
iae
itae
ce
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%  Proposed method   (modified NFTSM)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wko=1*rand(1,3);%10x4
mjk=1*rand(3,1);%3x10
bjk=1*rand(3,1);
djk=1*rand(3,1);
a=1*rand(1,1);
na=0.1;nb=0.1;nw=0.1;nd=0.1;nm=0.1;h_1=zeros(3,1);

xi=0;
n=1;iae=0;itae=0;ce=0;
t=0;
Dt=0.0001;

x1=0;
x2=0;x3=0;
e0=0;tk=[];Ts=[];jko=1;Ts(1)=0;tk(1)=0;
b=50;
z=00;a3=0;
u1=0;u_1=0;
e13=0;e1=0;a1=0;a2=0;a3=0;
c=2;%10
m=5/3;fp=0;ut=0;vT=01.5;uT=0.0;
k1=0.1;tt=0;E14=0;
k2=0.01;k=0;kt=0;m11=18/19;n11=32/33;%33,22,18
k3=350;%350
k4=400;%400
u1=2;u2=0.5;l1=0.5;l2=0.1;l3=0.1;
dx1=0.0;v1=0;v2=0;
delta=0;TEm3=0;RR3=0;

    u1=1;u2=0.1;l1=0.1;l2=0.1;l3=0.1;
for i=1:200000

    d=0.1*sin(t);  %正弦误差    
    f=(9.81*sin(x1)-0.492*0.21*0.65*x2^2*sin(2*x1))/(0.65*2/3-0.492*0.21*0.65/2*(cos(x1))^2);
    g=(-0.492*1*cos(x1))/(0.65*2/3-0.492*0.21*0.65/2*(cos(x1))^2);
    p=b*x2;
    dz=-b*z+b*(-b*x2-f-g*(u_1));
    z=z+dz*Dt;
    d_est=z+p;

    yd=1*cos(6*0.2*t)*cos(0.2*t); %输入指令 期望信号20°
    dyd=-1*6*0.2*sin(6*0.2*t)*cos(0.2*t)-1*0.2*cos(6*0.2*t)*sin(0.2*t); %输入信号的一阶导数
    ddyd=-1*6*0.2*6*0.2*cos(6*0.2*t)*cos(0.2*t)+1*6*0.2*0.2*sin(6*0.2*t)*sin(0.2*t)-1*6*0.2*0.2*sin(6*0.2*t)*sin(0.2*t)+1*0.2*0.2*cos(6*0.2*t)*cos(0.2*t); %输入信号的二阶导数
    y=x1; %实际输出
    dy=x2; %实际输出的一阶导数
    e1=y-yd;
    e2=dy-dyd;
  if t<0.753
    p1=(1.15^(2/13)+1.3-1.30*exp(5*2/13*t))^(13/2)+0.05;%0.002
  else
      p1=0.05;
  end
%      p1=(1.2-0.05)*exp(-11*t)+0.05;%0.002 
%     T=9;
%      if t<T
%         p1=(0.02^1-0.0011*t)^(7)+0.0005;
%      else
%          p1=0.0005;
%      end
%      dp1=-1.4*(0.2^7-0.2*t)^(6);
%      dp1=-11*(0.02-0.001)*exp(-11*t);
%      nn=e1/p1;
%      uu=tan(p1*nn/2);
%      re=pi*(1+uu^2)/2/p1;
%      uu1=re*(e2-2*dp1*atan(uu)/pi);
%      ddp1=121*(0.02-0.001)*exp(-11*t);
% %      ddp1=1.68*(0.2^7-0.2*t)^(5);
%      F=pi*(2*uu*uu1*p1-dp1*(1+uu^2))*(e2-2*dp1*atan(uu)/pi)/2/p1^2-re*(2*ddp1*atan(uu)/pi+2*dp1*uu1/(1+uu^2)/pi);

    nn=e1/p1;
    uu=0.5*log((1+nn)/(1-nn));
    re=1/(1-nn^2)/p1;
    dp1=-11*(1.2-0.05)*exp(-11*t);%转换误差一阶导z;0.0021,8
%     dp1=-7.5*exp(5*2/13*t)*(1.15^(2/13)+1.3-1.30*exp(5*2/13*t))^(13/2-1);
    uu1=re*(e2-e1*dp1/p1);%转换误差一阶导dz
    ddp1=121*(1.2-0.05)*exp(-11*t);
%     ddp1=7.5^2*(11/13)*exp(5*4/13*t)*(1.15^(2/13)+1.3-1.30*exp(5*2/13*t))^(13/2-2)-7.5*10/13*exp(5*2/13*t)*(1.15^(2/13)+1.3-1.30*exp(5*2/13*t))^(13/2-1);
        dre=-2*nn*dp1/(1-nn^2)/p1^2;
    F=dre*(e2-dp1*e1/p1)+re*(-e2*dp1/p1-e1*ddp1/p1+e1*dp1^2/p1^2);

    c1=0.8510;c2=0.051;k1=2.0;k2=30;k3=350;k4=500;%1,200
    r=10;
    r1=0.5;B=33.000001/33;k=8;E=0.1;BB=2;%ks;p
%     uu=e1;uu1=e2;re=0;
    s=uu+c1*log(r*abs(uu)+1)*sign(uu)+c2*(abs(uu1))^(1*B)*sign(uu1)+0.8*abs(uu)^m11*(tanh((uu^n11)/0.08))*1;
    
%     k1=0.8;k2=1/0.08;m1=17/19;n1=19/33;
%         s=uu+c1*(abs(uu))^r*sign(uu)+c2*(abs(uu1))^(1*B)*sign(uu1)+0.8*1*(tanh(uu^1*k2))*1;
% %     s=uu+c1*(abs(uu))^r*sign(uu)+c2*(abs(uu1))^(1*B)*sign(uu1)+k1*(abs(uu))^m1*(tanh(uu^n1*k2))*1;
% %     df=k1*m1*uu^(m-1)*tanh(uu^n1*k2)*uu1+k1*k2*n*uu^(m+n-1)*(1-(tanh(uu^n1*k2))^2)*uu1;
%     df=10*(1-(tanh(uu^1*k2))^2)*uu1;
%     
    s1=s;
    xi=s;

    
    %%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %神经网络输入输出计算
    net1=xi;
    for j=1:3
        for k=1:1
            ujk(j,k)=xi(k)+h_1(j,k)*bjk(j,k);
        end
    end
    
    zjk=djk.^2.*(ujk-mjk).^2;
    h=(1-djk.^2.*ujk.^2).*exp(-zjk);%小波输出
    for j=1:3
        f1(j)=(1-djk(j,1)^2*ujk(j,1)^2)*exp(-zjk(j,1))*1;%乘积
            
    end
    temp=0;
    for i=1:3
        temp=temp+wko(i)*f1(i);%wko=1x7
    end
    output4=temp;temp=0;
    for k=1:1
        temp=temp+a(k)*xi(k);%wko=1x7
    end
    out=output4+temp;  %小波输出
k=1.0;
%%
% o=0.5;
ww=wko(1,1)^2+wko(1,2)^2+wko(1,3)^2;%权值更新
ff=f1(1)^2+f1(2)^2+f1(3)^2;   
aa=a^2;
a_w=.1;a_a=.1;
%控制器设计
df=0.8*m11*(uu)^(m11-1)*(tanh((uu^n11)/0.08))*0*uu1+0.8/0.08*n11*abs(uu)^(n11+m11-1)*(1-(tanh(uu/0.08))^2)*uu1;
%    FF1=dre*uu1/re+re*(-e2*dp1/p1-e1*ddp1/p1+e1*dp1^2/p1^2); 
   F2=c2*B*(abs(uu1))^(B-1);
   F1=1*uu1+c1*r*uu1/(r*abs(uu)+1)+df+F2*(F-re*ddyd);
   ur=tt*tanh(F2*re*s/0.01)/2;
   drt=F2*re*s*0.05+F2*re*s*ww*ff/2/a_w^2+F2*re*s*aa*xi^2/2/a_a^2+0*ur+1*(a1+a2*abs(uu)+a3*abs(uu1)+0.520)*tanh(F2*re*s/0.01);
   tyr=drt*tanh(F2*re*s*drt/0.01);
   
   
   v=-(1+l1)*((F1+F2*re*f+s/0.64^2)/F2/re+u1*tanh(F2*re*s*u1/0.01)+1*tyr*tanh(F2*re*s*tyr/0.02)+0/(F2*re)*(s)^(1-2*0.02)*exp(s^(2*0.02))/0.02/5);
%    if abs(v-ut)<l1*abs(ut)+u2
%         ut=(v-l3*u2)*(1+l1*l2)/1;
%    else
%        ut=v;
%    end
%    u2=0.05;l1=0.5;
   if 1*abs(vT-uT)>=l1*abs(uT)+u2
        jko=jko+1;vT=v;
        tk(jko)=t-Ts(jko-1);
        Ts(jko)=t;
   end
   uT=(v-l3*u2)*(1+l1*l2)/1;
   ut=(v-l3*u2)*(1+l1*l2)/1;

    dk=-F2*re*s*ut;
    kt=kt+dk*Dt*1;
    bb1=0.52*(F2*re*s*tanh(F2*re*s/0.02)-0.1*tt)/2;%0.02;0.1
    tt=tt+bb1*Dt;
   %参数自适应 
    ww=ww+0.02*(F2^2*re^2*s^2*ff/2/a_w^2-0.4*ww)*Dt;
    aa=aa+0.02*(F2^2*re^2*s^2*xi^2/2/a_a^2-0.4*aa)*Dt;
    a1=a1+0.005*((re*(s)*(abs(uu1))^(B-1)+0.52)*tanh(F2*re*s/0.01)-0.1*a1)*Dt;
    a2=a2+0.005*(re*(s)*abs(uu)*(abs(uu1))^(B-1)*tanh(F2*re*s/0.01)-0.1*a2)*Dt;
    a3=a3+0.005*(re*(s)*(abs(uu1))^(B)*tanh(F2*re*s/0.01)-0.1*a3)*Dt;
    
   u_BP1=(((kt^2+2)*exp(kt^2/2)*sin(kt)+1)*ut/g);N=((kt^2+2)*exp(kt^2/2)*sin(kt)+1);%kt^2*cos(kt);
%    u_BP1=((1+kt)*cos(pi/2*(log(1+kt))^0.5)*ut/g);N=(1+kt)*cos(pi/2*(log(1+kt))^0.5);%kt^2*cos(kt);
%    u_BP1=((1+2*kt^2)*exp(kt^2)*sin(2^(n-0.5)*kt)+2^(n-0.5)*kt*exp(kt^2)*cos(2^(n-0.5)*kt))*ut/g;N=(1+2*kt^2)*exp(kt^2)*sin(2^(n-0.5)*kt)+2^(n-0.5)*kt*exp(kt^2)*cos(2^(n-0.5)*kt);%kt^2*cos(kt);
%%
    u=u_BP1;
    fpp=u^2;
    fp=fp+fpp*Dt;
    %%%倒立摆模型
    f=(9.81*sin(x1)-0.492*0.21*0.65*x2^2*sin(2*x1))/(0.65*2/3-0.492*0.21*0.65/2*(cos(x1))^2);
    g=(-0.492*1*cos(x1))/(0.65*2/3-0.492*0.21*0.65/2*(cos(x1))^2);
    Dx1=x2;
    Dx2=g*(u)+1*f+1*d;

%%
%     ds=c2*B*(abs(uu1))^(B-1)*(re*d+f*re);
ds=d+f;out= d_est;
    x1=x1+Dx1*Dt;
    x2=x2+Dx2*Dt;
    E14=E14+e1^2;
    if a3<=abs(e1)
        a3=abs(e1);end
%      p1=(0.015-0.0005)*exp(-11*t)+0.0005;%0.2,0.0008,11
    out_store(n)=out;
    d1_store(:,n)=[ds;out];
    dy1_store(:,n)=[yd;y];
    du1_store(n)=u;ds1_store(n)=s1;dv1_store(n)=v;dut1_store(n)=ut;
    de1_store(n)=e1;N1_store(n)=N;kt1_store(n)=kt;
    ddp11_store(n)=p1;
     w1_store(n)=ww;
    if n>=2000
        iae=iae+abs(e1)*Dt;
        itae=itae+t*abs(e1)*Dt;
        ce=ce+uu^2*Dt;
    end
%          yd11_store(:,n)=[yd1;yd];
%     T3=(de3_store(n)^2+de1_store(n)^2)^0.5;
%     TEm1=TEm1+T3/200000;
%     rr=T3-TEm1;
%     RR1=RR1+(rr)^0.5/200000;
    u_1=u;
    n=n+1; % u2=u2+(-0.001*u2-0.7*abs(vT-uT)+l1*abs(uT))*Dt;
    t=t+Dt;
end
iae
itae
ce
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iae=0;itae=0;ce=0;Ts1=Ts;tk1=tk;jko
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(88)
% scatter(Ts,tk);xlabel('Time/s');ylabel('u/m')
wko=0.5*rand(2,3);%10x4
mjk=1*rand(3,2);%3x10
bjk=1*rand(3,2);
djk=1*rand(3,2);
a=1*rand(1,2);
na=0.1;nb=0.1;nw=0.1;nd=0.1;nm=0.1;h_1=zeros(3,2);

xi=0;
n=1;
t=0;
Dt=0.0001;
tk=[];Ts=[];jko=1;Ts(1)=0;tk(1)=0;
x1=0;
x2=0;x3=0;
e0=0;
b=100;
z=00;a3=0;
u1=0;u_1=0;
e13=0;e1=0;a1=0;a2=0;a3=0;
c=2;%10
m=5/3;fp=0;
k1=0.1;tt=0;E14=0;
k2=0.01;k=0;kt=0;
k3=350;%350
k4=400;%400
dx1=0.0;v1=0;v2=0;ut=0;vT=0.1;uT=0;
delta=0;TEm3=0;RR3=0;u1=1;u2=0.1;l1=0.1;l2=0.1;l3=0.1;
for i=1:200000


    d=0.1*sin(t);    
    f=(9.81*sin(x1)-0.492*0.21*0.65*x2^2*sin(2*x1))/(0.65*2/3-0.492*0.21*0.65/2*(cos(x1))^2);
    g=(-0.492*1*cos(x1))/(0.65*2/3-0.492*0.21*0.65/2*(cos(x1))^2);
    p=b*x2;
    dz=-b*z+b*(-b*x2-f-g*(u_1));
    z=z+dz*Dt;
    d_est=z+p;
    yd=1*cos(6*0.2*t)*cos(0.2*t); %输入指令 期望信号20°
    dyd=-1*6*0.2*sin(6*0.2*t)*cos(0.2*t)-1*0.2*cos(6*0.2*t)*sin(0.2*t); %输入信号的一阶导数
    ddyd=-1*6*0.2*6*0.2*cos(6*0.2*t)*cos(0.2*t)+1*6*0.2*0.2*sin(6*0.2*t)*sin(0.2*t)-1*6*0.2*0.2*sin(6*0.2*t)*sin(0.2*t)+1*0.2*0.2*cos(6*0.2*t)*cos(0.2*t); %输入信号的二阶导数
    y=x1; %实际输出
    dy=x2; %实际输出的一阶导数
    e1=y-yd;
    e2=dy-dyd;
 if t<0.753
    p1=(1.15^(2/13)+1.3-1.30*exp(5*2/13*t))^(13/2)+0.05;%0.002
  else
      p1=0.05;
  end
%     p1=(1.2-0.05)*exp(-11*t)+0.05;%0.002
%     T=9;
%      if t<T
%         p1=(0.02^1-0.0011*t)^(7)+0.0005;
%      else
%          p1=0.0005;
%      end
%      dp1=-1.4*(0.2^7-0.2*t)^(6);
%      dp1=-11*(0.02-0.001)*exp(-11*t);
%      nn=e1/p1;
%      uu=tan(p1*nn/2);
%      re=pi*(1+uu^2)/2/p1;
%      uu1=re*(e2-2*dp1*atan(uu)/pi);
%      ddp1=121*(0.02-0.001)*exp(-11*t);
% %      ddp1=1.68*(0.2^7-0.2*t)^(5);
%      F=pi*(2*uu*uu1*p1-dp1*(1+uu^2))*(e2-2*dp1*atan(uu)/pi)/2/p1^2-re*(2*ddp1*atan(uu)/pi+2*dp1*uu1/(1+uu^2)/pi);

    nn=e1/p1;
    uu=0.5*log((1+nn)/(1-nn));
    re=1/(1-nn^2)/p1;
    dp1=-11*(1.2-0.05)*exp(-11*t);%转换误差一阶导z;0.0021,8
    uu1=re*(e2-e1*dp1/p1);%转换误差一阶导dz
    ddp1=121*(1.2-0.05)*exp(-11*t);

        dre=-2*nn*dp1/(1-nn^2)/p1^2;
    F=dre*(e2-dp1*e1/p1)+re*(-e2*dp1/p1-e1*ddp1/p1+e1*dp1^2/p1^2);

    c1=0.8510;c2=0.051;k1=2.0;k2=30;k3=350;k4=500;%1,200
    r=5;   
    r=10;
    r1=0.5;B=33.000001/33;k=8;E=0.1;BB=2;%ks;p
%     uu=e1;uu1=e2;re=0;
    s=uu+c1*log(r*abs(uu)+1)*sign(uu)+c2*(abs(uu1))^(1*B)*sign(uu1);
    
%     k1=0.8;k2=1/0.08;m1=17/19;n1=19/33;
%         s=uu+c1*(abs(uu))^r*sign(uu)+c2*(abs(uu1))^(1*B)*sign(uu1)+0.8*1*(tanh(uu^1*k2))*1;
% %     s=uu+c1*(abs(uu))^r*sign(uu)+c2*(abs(uu1))^(1*B)*sign(uu1)+k1*(abs(uu))^m1*(tanh(uu^n1*k2))*1;
% %     df=k1*m1*uu^(m-1)*tanh(uu^n1*k2)*uu1+k1*k2*n*uu^(m+n-1)*(1-(tanh(uu^n1*k2))^2)*uu1;
%     df=10*(1-(tanh(uu^1*k2))^2)*uu1;
%     
    s1=s;
    xi=[s,e1]';

    %2;0.5;0.5
    %%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %神经网络输入输出计算
    net1=xi;
    for j=1:3
        for k=1:2
            ujk(j,k)=xi(k)+h_1(j,k)*bjk(j,k);
        end
    end
    
    zjk=djk.^2.*(ujk-mjk).^2;
    h=(1-djk.^2.*ujk.^2).*exp(-zjk);%小波输出
    for j=1:3
        f1(j)=(1-djk(j,1)^2*ujk(j,1)^2)*exp(-zjk(j,1))*(1-djk(j,2)^2*ujk(j,2)^2)*exp(-zjk(j,2));%乘积
            
    end
    temp=0;
    for i=1:3
        temp=temp+wko(i)*f1(i);%wko=1x7
    end
    output4=temp;temp=0;
    for k=1:2
        temp=temp+a(k)*xi(k);%wko=1x7
    end
    out=output4+temp;  %小波输出

%%
%o=0.5;
ww=wko(1,1)^2+wko(1,2)^2+wko(1,3)^2+wko(2,1)^2+wko(2,2)^2+wko(2,3)^2;
ff=f1(1)^2+f1(2)^2+f1(3)^2;   
aa=a(1,1)^2+a(1,2)^2;
a_w=.1;a_a=.1;xx=xi(1,1)^2+xi(2,1)^2;
%    FF1=dre*uu1/re+re*(-e2*dp1/p1-e1*ddp1/p1+e1*dp1^2/p1^2); 
   F2=c2*B*(abs(uu1))^(B-1);
   F1=1*uu1+c1*r*uu1/(1*abs(uu)+1)+F2*(F-re*ddyd);
   ur=tt*tanh(F2*re*s/0.01)/2;
   drt=F2*re*s*0.05+F2*re*s*ww*ff/2/a_w^2+F2*re*s*aa*xx/2/a_a^2+0*ur+1*(a1+a2*abs(uu)+a3*abs(uu1)+0.520)*tanh(F2*re*s/0.01);
   tyr=drt*tanh(F2*re*s*drt/0.01);
   dd=1*(a1+a2*abs(uu)+a3*abs(uu1)+0.520)*tanh(F2*re*s/0.01);
   
   v=-(1+l1)*((F1+F2*re*f+s/0.64^2)/F2/re+u1*tanh(F2*re*s*u1/0.01)+1*tyr*tanh(F2*re*s*tyr/0.02));
  
   if 1*abs(vT-uT)>=l1*abs(ut)+u2
        jko=jko+1;vT=v;
        tk(jko)=t-Ts(jko-1);
        Ts(jko)=t;
   end
   uT=(v-l3*u2)*(1+l1*l2)/1;
   ut=(v-l3*u2)*(1+l1*l2)/1;
    dk=-F2*re*s*ut;
    kt=kt+dk*Dt*1;
    bb1=0.52*(F2*re*s*tanh(F2*re*s/0.02)-0.1*tt)/2;%0.02;0.1
    tt=tt+bb1*Dt;
        g0=F2*re*g*(u);
    da=g0*(s1).*xi; 
    for j=1:3
        dm(j,1)=s1*(g0)*wko(j)*h(j,1)*djk(j,1)^2*(ujk(j,1)-mjk(j,1))*(1-djk(j,2)^2*ujk(j,2)^2)*exp(-zjk(j,2));
    end
    for j=1:3 
            dm(j,2)=s1*(g0)*wko(j)*h(j,2)*djk(j,1)^2*(ujk(j,1)-mjk(j,1))*(1-djk(j,1)^2*ujk(j,1).^2)*exp(-zjk(j,1));
    end
    for j=1:3
            d_d(j,1)=2*s1*(g0)*wko(j)*djk(j,1)*((ujk(j,1)-mjk(j,1))^2*h(j,1)+ujk(j,1)^2*exp(-zjk(j,1)))*(1-djk(j,2)^2*ujk(j,2)^2)*exp(-zjk(j,2));
    end
    for j=1:3    
            d_d(j,2)=2*s1*(g0)*wko(j)*djk(j,2)*((ujk(j,2)-mjk(j,2))^2*h(j,2)+ujk(j,2)^2*exp(-zjk(j,2)))*(1-djk(j,1)^2*ujk(j,1)^2)*exp(-zjk(j,1));
    end
    for j=1:3 
            db(j,1)=2*s1*(g0)*wko(j)*djk(j,1)^2*h(j,1)*((ujk(j,1)-mjk(j,1))*h(j,1)+ujk(j,1)*exp(-zjk(j,1)))*(1-djk(j,2)^2*ujk(j,2)^2)*exp(-zjk(j,2));
    end
    for j=1:3  
            db(j,2)=2*s1*(g0)*wko(j)*djk(j,2)^2*h(j,2)*((ujk(j,2)-mjk(j,2))*h(j,2)+ujk(j,2)*exp(-zjk(j,2)))*(1-djk(j,1)^2*ujk(j,1)^2)*exp(-zjk(j,1));
    end
    bjk=bjk+nb*db*Dt;
    mjk=mjk+nm*dm*Dt;
%     wko=wko+nw*s1*g0*f1*Dt;
%     a=a+na*da*Dt;
    djk=djk+nd*d_d*Dt;
    ff;
       u2=u2+(-0.001*u2-abs(vT-uT)+l1*abs(uT))*Dt;
    ww=ww+0.02*(F2^2*re^2*s^2*ff/2/a_w^2-0.4*ww)*Dt;
    aa=aa+0.02*(F2^2*re^2*s^2*xx/2/a_a^2-0.4*aa)*Dt;
    a1=a1+0.005*(re*(s)*(abs(uu1))^(B-1)*tanh(F2*re*s/0.01)-0.4*a1)*Dt;
    a2=a2+0.005*(re*(s)*abs(uu)*(abs(uu1))^(B-1)*tanh(F2*re*s/0.01)-0.4*a2)*Dt;
    a3=a3+0.005*(re*(s)*(abs(uu1))^(B)*tanh(F2*re*s/0.01)-0.4*a3)*Dt;
    
   u_BP1=(((kt^2+2)*exp(kt^2/2)*sin(kt)+1)*ut/g);N=((kt^2+2)*exp(kt^2/2)*sin(kt)+1);%kt^2*cos(kt);
%%
    u=u_BP1;
    fpp=u^2;
    fp=fp+fpp*Dt;
    %%%倒立摆模型
    f=(9.81*sin(x1)-0.492*0.21*0.65*x2^2*sin(2*x1))/(0.65*2/3-0.492*0.21*0.65/2*(cos(x1))^2);
    g=(-0.492*1*cos(x1))/(0.65*2/3-0.492*0.21*0.65/2*(cos(x1))^2);
    Dx1=x2;
    Dx2=g*(u)+1*f+1*d;

%%
%     ds=c2*B*(abs(uu1))^(B-1)*(re*d+f*re);
ds=d+f;out= d_est;
    x1=x1+Dx1*Dt;
    x2=x2+Dx2*Dt;
    E14=E14+e1^2;
    if a3<=abs(e1)
        a3=abs(e1);end
%      p1=(0.015-0.0005)*exp(-11*t)+0.0005;%0.2,0.0008,11
    out1_store(n)=out;
d_store(n)=d;
    dy11_store(:,n)=[yd;y];
    du11_store(n)=u;ds11_store(n)=s1;dv11_store(n)=v;dut11_store(n)=ut;
    de11_store(n)=e1;N_store(n)=N;kt_store(n)=kt;
    ddp1_store(n)=p1;
     w11_store(n)=ww;a11_store(n)=aa;
    if n>=2000
        iae=iae+abs(e1)*Dt;
        itae=itae+t*abs(e1)*Dt;
        ce=ce+uu^2*Dt;
    end
%          yd11_store(:,n)=[yd1;yd];
%     T3=(de3_store(n)^2+de1_store(n)^2)^0.5;
%     TEm1=TEm1+T3/200000;
%     rr=T3-TEm1;
%     RR1=RR1+(rr)^0.5/200000;
    u_1=u;
    n=n+1;%  u2=u2+(-0.001*u2-0.7*abs(vT-uT)+l1*abs(uT))*Dt;
    t=t+Dt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iae
itae
ce
jko


figure(1)
plot((1:n-1)*Dt,dy1_store(1,:),'k',(1:n-1)*Dt,dy1_store(2,:),':r',(1:n-1)*Dt,yy1_store(2,:),'r',(1:n-1)*Dt,dy11_store(2,:),':b','linewidth',3);set(gcf,'color','w');xlabel('Time/s');ylabel('yd/m')
legend('y_d','Proposed method','ANFTSM[28]','Unmodified NFTSM');box off
figure(2)
plot((1:n-1)*Dt,ddp11_store,'k',(1:n-1)*Dt,-ddp11_store,'k',(1:n-1)*Dt,de1_store,'r',(1:n-1)*Dt,ee1_store,'b',(1:n-1)*Dt,de11_store,':k','linewidth',2);set(gcf,'color','w');xlabel('Time/s');ylabel('e/m')
legend('FTPF','-FTPF','Proposed method','ANFTSM[28]','Unmodified NFTSM');box off
figure(3)
plot((1:n-1)*Dt,du1_store,'r',(1:n-1)*Dt,uu1_store/1,'k',(1:n-1)*Dt,du11_store/1,'b','linewidth',2);set(gcf,'color','w');xlabel('Time/s');ylabel('u/N')
legend('Proposed method','ANFTSM[28]','无改进NFTSM');box off
figure(4)
plot((1:n-1)*Dt,ds1_store/6,'r',(1:n-1)*Dt,s1_store/1,'k',(1:n-1)*Dt,ds11_store/6,'b','linewidth',2);set(gcf,'color','w');xlabel('Time/s');ylabel('s/m')
legend('Proposed method','ANFTSM[28]','无改进NFTSM');box off
figure(5)
plot((1:n-1)*Dt,N1_store,':r',(1:n-1)*Dt,kt1_store,':k','linewidth',2);set(gcf,'color','w');xlabel('Time/s');ylabel('N(\xi), \xi')
legend('N(\xi)','\xi');box off
figure(6)
stem(Ts1,tk1);xlabel('Time/s');ylabel('t_k_+_1-t_k/s')
legend("Event-triggered times");
figure(7)
stem(Ts,tk);xlabel('Time/s');ylabel('t_k_+_1-t_k/s')
legend("Event-triggered times");
% figure
% plot((1:n-1)*Dt,1*w1_store,'r',(1:n-1)*Dt,a11_store,'b','linewidth',2);xlabel('Time/s');ylabel('||w||^2,||a||^2')
% legend('||w||^2','||a||^2');
