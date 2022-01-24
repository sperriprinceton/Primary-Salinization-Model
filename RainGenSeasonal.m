%This function generates stochastic rain (compound poisson proccess)with
%seasonality
%Author Saverio Perri
%Created 18/04/2020

function [Rain,lambda_Pt,alpha_Pt] = RainGenSeasonal(t,s1,n,Zr,lambda_P,omega_p,phi_p,alpha_p,A)
%preallocate the variables
alpha_Pt=zeros(1,t);
gamma_P=zeros(1,t); 
lambda_Pt=zeros(1,t); 
  
for i=2:t
    
alpha_Pt(i)=(1+A*sin(omega_p*i+ phi_p))*alpha_p; % time-dependent average precipitation per event[mm]
gamma_P(i)=s1*n*Zr./alpha_Pt(i);  
lambda_Pt(i)=lambda_P*(1+A*sin(omega_p*i+ phi_p));

    r=rand;
     
    if r<lambda_Pt(i)
%   
      Rain(i)=(-log(rand)./gamma_P(i));

    else
        
      Rain(i)=0;

 end
end