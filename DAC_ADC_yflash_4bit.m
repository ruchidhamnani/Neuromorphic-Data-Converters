close all;
clc;
clear;
%%
%defining Y-flash parameters
Vthinit = [1,1,1,1];
Vth1 = 1;
Vth2 = 1;
Vth3 = 1;
Vth4 = 1;

error = 1;
Vth = Vthinit;
Iread = 1e-9;
mVT = 0.144765;
CRprog = 0.48;
CRread = 1;
b5 = 2.1e-4;
dt = 1e-6;
%Writing pulse for 10us
T_write = 10e-6;
%%
%Feedback resistor
Rf = 4.05e6;
%Reference voltage for DAC
Vref = 0.1125;
alpha = 0.1125;
D = 1.8*[0 0 0 0; 0 0 0 1; 0 0 1 0; 0 0 1 1; 0 1 0 0; 0 1 0 1; 0 1 1 0; 0 1 1 1; 
         1 0 0 0; 1 0 0 1; 1 0 1 0; 1 0 1 1; 1 1 0 0; 1 1 0 1; 1 1 1 0; 1 1 1 1];
vin = [0:Vref:1.6875];

Vout = zeros(1,4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_rate_gd = 1;
Vdig_ideal = 0;
n = 0;
A_gd = 0;                           %analog output of the GD algorithm,  in ADC based feedback the analog aoutput is considered
tgd=0;                              %time of gd algorithm
W_prev = [0.4, 0.3, 0.2, 0.2];      %initial width of the memristor device
max_DNL = 0;
max_INL = 0;

for l = 1:200
for x = 1:16
    k = x; 
    n = n+1;

A_gd(end+1) = (Rf*(Iread*(exp(D(x,1)/mVT))*exp(-Vth(1)/mVT)*(D(x,1)/1.8) + Iread*(exp(D(x,2)/mVT))*exp(-Vth(2)/mVT)*(D(x,2)/1.8) + Iread*(exp(D(x,3)/mVT))*exp(-Vth(3)/mVT)*(D(x,3)/1.8) + Iread*(exp(D(x,4)/mVT))*exp(-Vth(4)/mVT)*(D(x,4)/1.8)));
A_gd_static(x) = (Rf*(Iread*(exp(D(x,1)/mVT))*exp(-Vth(1)/mVT)*(D(x,1)/1.8) + Iread*(exp(D(x,2)/mVT))*exp(-Vth(2)/mVT)*(D(x,2)/1.8) + Iread*(exp(D(x,3)/mVT))*exp(-Vth(3)/mVT)*(D(x,3)/1.8) + Iread*(exp(D(x,4)/mVT))*exp(-Vth(4)/mVT)*(D(x,4)/1.8)));


%%
%ideal ADC conversion in the feedback
if(A_gd(end) < Vref)
    Vdig_ideal = [0 0 0 0];
elseif((2*Vref > A_gd(end)) &&  A_gd(end) >= Vref*1)
    Vdig_ideal = [0 0 0 1];
elseif((3*Vref > A_gd(end)) &&  A_gd(end) >= Vref*2)
    Vdig_ideal = [0 0 1 0];
elseif((4*Vref > A_gd(end)) &&  A_gd(end) >= Vref*3)
    Vdig_ideal = [0 0 1 1];
elseif((5*Vref > A_gd(end)) &&  A_gd(end) >= Vref*4)
    Vdig_ideal = [0 1 0 0];
elseif((6*Vref > A_gd(end)) &&  A_gd(end) >= Vref*5)
    Vdig_ideal = [0 1 0 1];
elseif((7*Vref > A_gd(end)) &&  A_gd(end) >= Vref*6)
    Vdig_ideal = [0 1 1 0];
elseif((8*Vref > A_gd(end)) &&  A_gd(end) >= Vref*7)
    Vdig_ideal = [0 1 1 1];
elseif((9*Vref > A_gd(end)) &&  A_gd(end) >= Vref*8)
    Vdig_ideal = [1 0 0 0];
elseif((10*Vref > A_gd(end)) &&  A_gd(end) >= Vref*9)
    Vdig_ideal = [1 0 0 1];
elseif((11*Vref > A_gd(end)) &&  A_gd(end) >= Vref*10)
    Vdig_ideal = [1 0 1 0];
elseif((12*Vref > A_gd(end)) &&  A_gd(end) >= Vref*11)
    Vdig_ideal = [1 0 1 1];
elseif((13*Vref > A_gd(end)) &&  A_gd(end) >= Vref*12)
    Vdig_ideal = [1 1 0 0];
elseif((14*Vref > A_gd(end)) &&  A_gd(end) >= Vref*13)
    Vdig_ideal = [1 1 0 1];
elseif((15*Vref > A_gd(end)) &&  A_gd(end) >= Vref*14)
    Vdig_ideal = [1 1 1 0];
elseif(A_gd(end) >= 15*Vref)
    Vdig_ideal = [1 1 1 1];
end
Vdig_ideal = Vdig_ideal*1.8;                                %converting the digital values to analog voltages

%%
%Vdig_ideal goes from MSB to LSB <-- left to right
%calculating the error function for the ADC with feedback
error(end+1) = error(end) + 0.5*((D(k,1)-Vdig_ideal(1))^2 +(D(k,2)-Vdig_ideal(2))^2 +(D(k,3)-Vdig_ideal(3))^2 +(D(k,4)-Vdig_ideal(4))^2 );
tgd(end+1) = tgd(end) + 1;
E_rate_gd(end+1) = (1/n)*error(end);
%%
end

%updating the memristors
for i = 1:1:4
    W_f = W_prev(i)*(1-W_prev(i));                          %this is a simple window function
        if(D(k,i)-Vdig_ideal(i)<0)

            Vw = 5;
                        dVthdt(i) = -(vin(x) - A_gd(end))*b5*(CRprog*Vw - Vth(i))*D(x,i)/1.8;
                        Vth(i) = Vth(i) + dVthdt(i);
            
        end
        if(D(k,i)-Vdig_ideal(i)>0)

             dVthdt(c) = (vin(x) - A_gd(end))*4.643e-4*(0.9531 -  Vth(c))*exp(-(0.07/0.9531-Vth(c)))*D(x,c)/1.8;
                        Vth(c) = Vth(c) + dVthdt(c);
        end

end
    Vth1(end+1) = Vth(1);
    Vth2(end+1) = Vth(2);
    Vth3(end+1) = Vth(3);
    Vth4(end+1) = Vth(4);   
end
vin = [0:0.1125:1.6875];
A_static=zeros(1,1*length(vin));

for i=1:length(vin)
    for j=1:1
       A_gd_static(i) = Rf*(Iread*(exp(D(i,1)/mVT))*exp(-Vth(1)/mVT)*(D(i,1)/1.8) + Iread*(exp(D(i,2)/mVT))*exp(-Vth(2)/mVT)*(D(i,2)/1.8) + Iread*(exp(D(i,3)/mVT))*exp(-Vth(3)/mVT)*(D(i,3)/1.8) + Iread*(exp(D(i,4)/mVT))*exp(-Vth(4)/mVT)*(D(i,4)/1.8));
           end
end

out = A_gd_static;

for i=1:length(out)
    INL(i) = (out(i))/Vref - (i-1);
end
for i=1:length(out)-1
    DNL(i) = (out(i+1) - out(i))/Vref -1;
end

maxINL=max(abs(INL));
maxDNL=max(abs(DNL));

%Plotting DNL
figure(1)
plot((0:14),DNL,'b-^','linewidth',2);
xticks([0 2 4 6 8 10 12 14])
xticklabels({'0','2','4','6','8','10','12', '14'})
hold on
ylabel('DNL (LSB)');
xlabel('Input code')
box on

%Plotting INL
figure(2)
plot((0:15),INL,'r-^','linewidth',2);
xticks([0 2 4 6 8 10 12 14])
xticklabels({'0','2','4','6','8','10','12', '14'})
hold on
ylabel('INL (LSB)');
xlabel('Input code')
box on
