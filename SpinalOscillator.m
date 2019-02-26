clc
clear all
close all

%initialization of variables 
ksyn=2;
gsyn=0.3;
Vpir=120;
VL=-60;
Vsyn=-80;
gL=0.1;
tau=10;
phi=3;
thetasyn=-44;
gpir=0.3;
C=1;
dt=0.1;
t=(0:dt:500)';
num=length(t);

%From the graph 
V1=zeros(num,1);
V2=zeros(num,1);
V11=-45;
V22=-60;

%Corresponding initial values for h
h1=zeros(num,1);
h2=zeros(num,1);
h11=0.2;
h22=0.1;

minf1=0;
minf2=0;

for j=1:num
   
   h1(j)=h11;
   h2(j)=h22;
   V1(j)=V11;
   V2(j)=V22;
   
   %equations for the two cells
   minf1=1/(1+exp(-(V1(j)+65)/7.8));
   hinf1=1/(1+exp((V1(j)+81)/11));
   tauh1=hinf1*exp((V1(j)+162.3)/17.8);
   s12=1/(1+exp(-(V1(j)-thetasyn)/ksyn));
  
   minf2=1/(1+exp(-(V2(j)+65)/7.8));
   hinf2=1/(1+exp((V2(j)+81)/11));
   tauh2=hinf2*exp((V2(j)+162.3)/17.8);
   s21=1/(1+exp(-(V2(j)-thetasyn)/ksyn));
   
   V11=V11+dt*((-gpir*minf1^3*h1(j)*(V11-Vpir)-gL*(V11-VL)-gsyn*s21*(V11-Vsyn))/C);
   h11=h11+dt*(phi*(hinf1-h1(j))/tauh1);
   V22=V22+dt*((-gpir*minf2^3*h2(j)*(V11-Vpir)-gL*(V22-VL)-gsyn*s12*(V22-Vsyn))/C);
   h22=h22+dt*(phi*(hinf2-h2(j))/tauh2);
   
end

figure(1)
   subplot(211);
   plot(t',V1,'m')
   title('Voltage of cell 1 and 2 vs. Time')
   ylabel('V1 (mV)')
   xlabel('Time (ms)')
   
   subplot(212);
   plot(t',V2)
   ylabel('V2 (mV)')
   xlabel('Time (ms)')

   
   figure (2)
   plot([V1],[h1])
   title('Voltage vs. h for neuron 1')
   ylabel('V1 (mV)')
   xlabel('h1 ')
   

   
   
   
   



