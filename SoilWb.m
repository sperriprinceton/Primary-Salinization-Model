%This function models soil water and salt mass balance and with previously generated stochastic rain (compound
%poisson proccess)
%Author Saverio Perri
% the model solves the water balance 
%Last update 24/01/2021

function [Leak,M,sT,Theta,E,s,C,CMax,Leaks,Ms,Es,ss,Cs,ETmax] = SoilWb(t,beta,CT,n,Zr,Etmax,Rain,in,sw,omega_e,phi_e)

b=1; %dilution efficiency
s1=0.8; 

%variables preallocation
Leak=zeros(1,t); % leakage WITH SALINITY FEEDBACK [g/day]
Leak(1)=0; 
Leaks=zeros(1,t); % leakage WITHOUT SALINITY FEEDBACK[g/day]
Leaks(1)=0;
ss=zeros(1,t);
s=zeros(1,t);
E=zeros(1,t); % Evapotranspiration WITH SALINITY FEEDBACK [mm/day]
M=zeros(1,t); % Salt mass WITH SALINITY FEEDBACK[g/m^2]
C=zeros(1,t); % Salt concentration WITH SALINITY FEEDBACK [g/l]
Es=zeros(1,t);% Evapotranspiration WITHOUT SALINITY FEEDBACK [mm/day]
Ms=zeros(1,t); % Salt mass WITHOUT SALINITY FEEDBACK [g/m^2]
Cs=zeros(1,t); % Salt concentration WITHOUT SALINITY FEEDBACK [g/l]
Theta=zeros(1,t); % Water balance energy limit WITH SALINITY FEEDBACK [-]
sT=zeros(1,t); % salt stress threshold value[-]
C0=0; % initial concentration value g/l
C(1)=C0;
Cs(1)=C0;
s0=0.8; % Soil moisture at the beginnig of the simulation [-]
s(1)=s0;
E(1)=2.5/(n*Zr); % mm/day
M(1)=C0*s0*n*Zr;
ss(1)=s0;
Es(1)=2.5/(n*Zr); % mm/day
Ms(1)=C0*s0*n*Zr;
Theta(1)=(beta*M(1)/(n*Zr))/(1+beta*CT);
sT(1)=M(1)/(n*Zr*CT);

CMax=(1+beta*CT)/beta*ones(t,1);
  
  
for i=2:t
ETmax(i)=Etmax/(n*Zr)*(1+0.5*sin(omega_e*i+ phi_e)); %time-dependent maximum value of normalized ET [-]

   sT(i)=sT(i-1);
   M(i)=M(i-1)+in;
   Theta(i)=Theta(i-1);
             
  
   Ms(i)=Ms(i-1)+in;


       ss(i)=ss(i-1)+Rain(i);
       s(i)=s(i-1)+Rain(i);

      %%%% NO FEEDBACK MODEL STEPWISE FUNCTION%%%
    if ss(i)<=sw

     Es(i)=0;
%      ss(i)=ss(i)-Es(i)+Rain(i);
     ss(i)=ss(i)-Es(i);
     Ms(i)= Ms(i);
     Cs(i)=Ms(i)/(ss(i)*n*Zr);     
    
    elseif ss(i)>sw && ss(i)< s1

     Es(i)=ETmax(i)*(ss(i)/(s1));
%      ss(i)=ss(i)-Es(i)+Rain(i);
     ss(i)=ss(i)-Es(i);
     Ms(i)= Ms(i);
     Cs(i)=Ms(i)/(ss(i)*n*Zr); 
     Es(i)=ETmax(i)*(ss(i)/(s1));        
           
   elseif ss(i)>=s1 && ss(i)<=1

     Es(i)=ETmax(i);
     Leaks(i)=b*((ss(i)-s1))*Ms(i);
%      ss(i)=s1-Es(i)+Rain(i); 
     ss(i)=s1-Es(i); 
     Ms(i)=Ms(i)-Leaks(i);    
     Cs(i)=Ms(i)/(ss(i)*n*Zr); 
            
    elseif ss(i)> 1 
      Es(i)=ETmax(i);
      Leaks(i)=b*((1-s1))*Ms(i); 
%      ss(i)=s1-Es(i)+Rain(i); 
      ss(i)=s1-Es(i); 
      Ms(i)=Ms(i)-Leaks(i);    
      Cs(i)=Ms(i)/(ss(i)*n*Zr); 

     end
           
    if Ms(i)<=0
        Ms(i)=0;
    end
    if ss(i)>1
       ss(i)=1;       
    end


%%%% FEEDBACK MODEL STEPWISE FUNCTION%%%

      if s(i)<= Theta(i)
        
       M(i)= M(i);
       sT(i)= M(i)/(n*Zr*CT);
       Theta(i)=(beta*(M(i))/(n*Zr))/(1+beta*CT);
       E(i)=0;
%        s(i)=s(i)-E(i)+Rain(i);
       s(i)=s(i)-E(i);
       C(i)=M(i)/(n*Zr*s(i));
       
      elseif s(i)< s1 && s(i)<= sT(i)

       M(i)= M(i);
       sT(i)= M(i)/(n*Zr*CT);
       Theta(i)=((beta*M(i))/(n*Zr))/(1+beta*CT);
       E(i)=ETmax(i)/s1*(1+beta*CT)*(s(i)-Theta(i));
%        s(i)=s(i)-E(i)+Rain(i);
       s(i)=s(i)-E(i);
       C(i)=M(i)/(n*Zr*s(i));
       E(i)=ETmax(i)/s1*(1+beta*CT)*(s(i)-Theta(i));
              

    elseif s(i)< s1 && s(i)> sT(i)
       M(i)= M(i);
       sT(i)= M(i)/(n*Zr*CT);
       Theta(i)=((beta*M(i))/(n*Zr))/(1+beta*CT);
       E(i)=ETmax(i)*s(i)/s1;
%        s(i)=s(i)-E(i)+Rain(i);
       s(i)=s(i)-E(i);
       C(i)=M(i)/(n*Zr*s(i));
       
       
    elseif s(i)>=s1 && s(i)<=1
        
       Leak(i)=b*((s(i)-s1))*M(i);
       M(i)= M(i)-Leak(i);
       sT(i)= M(i)/(n*Zr*CT);
       Theta(i)=((beta*M(i))/(n*Zr))/(1+beta*CT);
       E(i)=ETmax(i);
%        s(i)=s1-E(i)+Rain(i);
       s(i)=s1-E(i);
       C(i)=M(i)/(n*Zr*s(i));
        
     elseif s(i)> 1 
     
      Leak(i)=b*((1-s1))*M(i);  
      M(i)= M(i)-Leak(i);      
      sT(i)= M(i)/(n*Zr*CT);
      Theta(i)=((beta*M(i))/(n*Zr))/(1+beta*CT);
      E(i)=ETmax(i);
%        s(i)=s1-E(i)+Rain(i);
      s(i)=s1-E(i);
      C(i)=M(i)/(n*Zr*s(i));
      
      end
     
    if M(i)<=0
        M(i)=0;
    else
    end
    if s(i)>1
       s(i)=1;
    end
    

    if M(i)<=0
        M(i)=0;
    else
    end
    if s(i)>1
       s(i)=1;
    end
    if E(i)<0
       E(i)=0;
    end
    
    
    end
end
% end