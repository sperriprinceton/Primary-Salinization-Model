%This function models soil water and salt mass balance and generates stochastic rain (compound
%poisson proccess)
%Author Saverio Perri
%Last update 25/04/2019

function [Leak,M,sT,Theta,E,s,C,CMax,Leaks,Ms,Es,ss,Cs,M_mean,C_mean,E_mean,s_mean] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rain,in,sw,omega_e,phi_e)
% sT= Relative soil moisture threshold 
% Theta= Energetic limit in the water balance
% E= Evaporation
% L= Leaching
% s= Realative soil moisture
% M= Salt mass
% C= Salt concentration
b=1;
% %%%% PHYSICAL VARIABLES%%%%%% 
% ETmax=Etmax/(n*Zr)*(1+0.5*sin(omega_e*t+ phi_e)); %time-dependent maximum value of normalized ET [-]
s1=0.8; 
% sw=0;% wilting point
% beta=0.084; % l/g average tolerance Maas et al.
% CT=3.84;% g/l average tolerance Maas et al.
% gamma_P=s1*n*Zr./alpha_P;
%%%%%%variables preallocation and initial values%%%%
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
% %Salt input parameters
% in=60*10^(-3);%salt input  g/day*m^2 in coastal area (200 kg/ha/yr)
% b=0.8;% manovella di samir
  tot=0;
CMax=(1+beta*CT)/beta*ones(t,1);
  
  
for i=2:t
ETmax(i)=Etmax/(n*Zr)*(1+0.5*sin(omega_e*i+ phi_e)); %time-dependent maximum value of normalized ET [-]

   sT(i)=sT(i-1);
   M(i)=M(i-1)+in;
   Theta(i)=Theta(i-1);
             
  
   Ms(i)=Ms(i-1)+in;

%     r=rand;
%      
%     if r<lambda_P
% %        tot=tot+1;
%        Rain(i)=(-log(rand)./gamma_P);
       ss(i)=ss(i-1)+Rain(i);
       s(i)=s(i-1)+Rain(i);

   %%%% NO FEEDBACK MODEL STEPWISE FUNCTION%%%
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
%     if Cs(i)>45
%        Cs(i)=45; % with sw=0 the concentration in the No-feedback system could reach unrealistic values
%     end
    

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
    
      %%%% No FEEDBACK MODEL STEPWISE FUNCTION%%%
%     else
%         
%          Rain(i)=0;
%          ss(i)=ss(i-1)+Rain(i);
%          s(i)=s(i-1)+Rain(i);
%            if ss(i)<=sw
% 
%      Es(i)=0;
% %      ss(i)=ss(i)-Es(i)+Rain(i);
%      ss(i)=ss(i)-Es(i);
%      Ms(i)= Ms(i);
%      Cs(i)=Ms(i)/(ss(i)*n*Zr);     
%     
%     elseif ss(i)>sw && ss(i)< s1
% 
%      Es(i)=ETmax*(ss(i)/(s1));
% %      ss(i)=ss(i)-Es(i)+Rain(i);
%      ss(i)=ss(i)-Es(i);
%      Ms(i)= Ms(i);
%      Cs(i)=Ms(i)/(ss(i)*n*Zr); 
%      Es(i)=ETmax*(ss(i)/(s1));
%              
%            
%    elseif ss(i)>=s1 && ss(i)<=1
% 
%      Es(i)=ETmax;
%      Leaks(i)=b*((ss(i)-s1))*Ms(i);
% %      ss(i)=s1-Es(i)+Rain(i); 
%      ss(i)=s1-Es(i); 
%      Ms(i)=Ms(i)-Leaks(i);    
%      Cs(i)=Ms(i)/(ss(i)*n*Zr); 
%             
%     elseif ss(i)> 1 
%       Es(i)=ETmax;
%       Leaks(i)=b*((1-s1))*Ms(i); 
% %      ss(i)=s1-Es(i)+Rain(i); 
%       ss(i)=s1-Es(i); 
%       Ms(i)=Ms(i)-Leaks(i);    
%       Cs(i)=Ms(i)/(ss(i)*n*Zr); 
% 
%     end
%            
%     if Ms(i)<=0
%         Ms(i)=0;
%     end
%     if ss(i)>1
%        ss(i)=1;       
%     end
%         if Es(i)<1
%        Es(i)=0;       
%     end
% %     if Cs(i)>45
% %        Cs(i)=45; % with sw=0 the concentration in the No-feedback system could reach unrealistic values
% %     end
%     
% 
% %%%% FEEDBACK MODEL STEPWISE FUNCTION%%%
% 
%       if s(i)<= Theta(i)
%         
%        M(i)= M(i);
%        sT(i)= M(i)/(n*Zr*CT);
%        Theta(i)=(beta*(M(i))/(n*Zr))/(1+beta*CT);
%        E(i)=0;
% %        s(i)=s(i)-E(i)+Rain(i);
%        s(i)=s(i)-E(i);
%        C(i)=M(i)/(n*Zr*s(i));
%        
%       elseif s(i)< s1 && s(i)<= sT(i)
% 
%        M(i)= M(i);
%        sT(i)= M(i)/(n*Zr*CT);
%        Theta(i)=((beta*M(i))/(n*Zr))/(1+beta*CT);
%        E(i)=ETmax/s1*(1+beta*CT)*(s(i)-Theta(i));
% %        s(i)=s(i)-E(i)+Rain(i);
%        s(i)=s(i)-E(i);
%        C(i)=M(i)/(n*Zr*s(i));
%        E(i)=ETmax/s1*(1+beta*CT)*(s(i)-Theta(i));
%        
%     elseif s(i)< s1 && s(i)> sT(i)
%        M(i)= M(i);
%        sT(i)= M(i)/(n*Zr*CT);
%        Theta(i)=((beta*M(i))/(n*Zr))/(1+beta*CT);
%        E(i)=ETmax*s(i)/s1;
% %        s(i)=s(i)-E(i)+Rain(i);
%        s(i)=s(i)-E(i);
%        C(i)=M(i)/(n*Zr*s(i));
%        E(i)=ETmax*s(i)/s1;
%        
%        
%     elseif s(i)>=s1 && s(i)<=1
%         
%        Leak(i)=b*((s(i)-s1))*M(i);
%        M(i)= M(i)-Leak(i);
%        sT(i)= M(i)/(n*Zr*CT);
%        Theta(i)=((beta*M(i))/(n*Zr))/(1+beta*CT);
%        E(i)=ETmax;
% %        s(i)=s1-E(i)+Rain(i);
%        s(i)=s1-E(i);
%        C(i)=M(i)/(n*Zr*s(i));
%        
%         
%      elseif s(i)> 1 
%      
%       Leak(i)=b*((1-s1))*M(i);  
%       M(i)= M(i)-Leak(i);      
%       sT(i)= M(i)/(n*Zr*CT);
%       Theta(i)=((beta*M(i))/(n*Zr))/(1+beta*CT);
%       E(i)=ETmax;
% %        s(i)=s1-E(i)+Rain(i);
%       s(i)=s1-E(i);
%       C(i)=M(i)/(n*Zr*s(i));
%       
%       end
     
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
    

  M_mean=nanmean(M);  
C_mean=nanmean(C);
E_mean=nanmean(E);
s_mean=nanmean(s);

end
% end