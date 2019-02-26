clc
clear all
close all

% Initial Parameters
gsyn=0.2;
gL=0.05;
gpir=0.5;
C=1;
tau=10;
phi=2;
thetasyn=-35;
ksyn=2;
kr=0.005;
Vpir=120;
VL=-60;
Vsyn=-80;
dt=0.1;
t=(0:dt:3000)';
num=length(t);

V1=zeros(num,1);
V2=zeros(num,1);
V11=-60;
V22=-72;
h1=zeros(num,1);
h2=zeros(num,1);
h11=0.2;
h22=0.15;
s12=1;
s21=0;
I1=heaviside(t-300)-heaviside(t-350)+heaviside(t-1100)-heaviside(t-1150);
I2=heaviside(t-300)-heaviside(t-350)-heaviside(t-1100)+heaviside(t-1150);

for j=1:num
   h1(j)=h11;
   h2(j)=h22;
   V1(j)=V11;
   V2(j)=V22;
   
   S12(j)=s12;
   S21(j)=s21;
   
  
   minf1=1/(1+exp(-(V1(j)+65)/7.8));
   hinf1=1/(1+exp((V1(j)+81)/11));
   tauh1=hinf1*exp((V1(j)+162.3)/17.8);
   sgiven12=1/(1+exp(-(V1(j)-thetasyn)/ksyn));
   
   s12=1+heaviside(dt*j-300)*(s12+dt*(sgiven12*(1-s12)-kr*s12)-1);
    h11=h11+dt*(phi*(hinf1-h1(j))/tauh1);
    
   minf2=1/(1+exp(-(V2(j)+65)/7.8));
   hinf2=1/(1+exp((V2(j)+81)/11));
   tauh2=hinf2*exp((V2(j)+162.3)/17.8);
   sgiven21=1/(1+exp(-(V2(j)-thetasyn)/ksyn));
   
   s21=heaviside(dt*j-300)*(s21+dt*(sgiven21*(1-s21)-kr*s21));
   h22=h22+dt*(phi*(hinf2-h2(j))/tauh2);

  V11=V11+dt*((-gpir*minf1^3*h1(j)*(V11-Vpir)-gL*(V11-VL)-gsyn*s21*(V11-Vsyn)+I1(j))/C);
   V22=V22+dt*((-gpir*minf2^3*h2(j)*(V11-Vpir)-gL*(V22-VL)-gsyn*s12*(V22-Vsyn)+I2(j))/C);
   
  
  
   
   
end

figure (1)
   subplot(211);
   plot(t',V1)
   title('Cell 1 Voltage vs. Time')
   xlabel('Time (ms)')
   ylabel('Voltage of the Cell 1 (mV)')
   xlim([-20 3000])
   
   subplot(212);
   plot(t',V2)
   title('Cell 2 Voltage vs. Time')
   xlabel('Time (ms)')
   ylabel('Voltage of the Cell 2 (mV)')
   xlim([-20 3000])
   
   

figure (2)
   subplot(611);
   plot(t',S12)
   title(' Variables of Neuron 1')
   ylabel('S12')
   ylim([-0.1 1.1])
   
   subplot(612);
   plot(t',V1)
   ylabel('V1 (mV)')

   subplot(613);
   plot(t',I1)
   ylabel('I1 (uA)')
   ylim([-1.1 1.1])
   
   subplot(614);
   plot(t',S21)
   title(' Variables of Neuron 2')
   ylabel('S21')
   ylim([-0.1 1.1])
   
   subplot(615);
   plot(t',V2)
   ylabel('V2 (mV)')

   subplot(616);
   plot(t',I2)
   xlabel('Time (ms)')
   ylabel('I2 (uA)')
   ylim([-1.1 1.1])
   
   



