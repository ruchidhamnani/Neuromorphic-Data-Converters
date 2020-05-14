clc;
clear all;
close all;

%
Vthinit = [1,1,1,1,1,1,1,1,1,1,1,1; 1,1,1,1,1,1,1,1,1,1,1,1;  1,1,1,1,1,1,1,1,1,1,1,1;  1,1,1,1,1,1,1,1,1,1,1,1; 
1,1,1,1,1,1,1,1,1,1,1,1;  1,1,1,1,1,1,1,1,1,1,1,1;  1,1,1,1,1,1,1,1,1,1,1,1;  1,1,1,1,1,1,1,1,1,1,1,1;];
%defining Y-flash parameters
Vth = Vthinit;
Iread = 1e-9;
mVT = 0.144765;
CRprog = 0.48;
CRread = 1;
b5 = 2.1e-4;
dt = 1e-6;


%Writing pulse for 10us
T_write = 10e-6;

%Feedback resistor
Rf = 1e3;


%Read Voltage Y-flash for Sub-threshold mode operation
Vr = 1.8;


 
%output bits 
Y = zeros(8,1);
Neurons = zeros(8,1);




for i = 1:1:8
    for j = 1:1:12
	    
	  I(i,j)= Iread*exp(CRread*(Vr-Vth(i,j))/0.144765);
       Y(i,1)=Y(i,1)+ (I(i,j)*Rf);
        
    end
    
   
end
Y(:,1)
        
	
	    
      
	   
           
                