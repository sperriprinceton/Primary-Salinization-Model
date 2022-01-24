%This function generates stochastic rain (compound poisson proccess)with
%seasonality
%Author Saverio Perri
%Created 18/04/2020

function [Rain] = RainGenNS(t,s1,n,Zr,lambda_P,alpha_P)
gamma_P=s1*n*Zr./alpha_P;  

for i=2:t

    r=rand;
     
    if r<lambda_P

       Rain(i)=(-log(rand)./gamma_P);

    else
        
         Rain(i)=0;

 end
end