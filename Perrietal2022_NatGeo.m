% This code includes a parsimonious model of primary (i.e., natural) 
% soil salinization. It solves the balance equations for the vertically-averaged 
% relative soil moisture and salt mass dynamics under stochastic rainfall inputs 
% and accounts for seasonality and vegetation-salinity feedbacks.
% Author: Saverio Perri (sperri@princeton.edu)

% Citations: 
% [1] Perri, S., Suweis, S., Entekhabi, D. & Molini, A. Vegetation controls on dryland salinity. Geophys. Res. Lett. 45, 669-682 (2018).
% [2] Perri, S., et al. River basin salinization as a form of aridity. Proc. Natl. Acad. Sci. USA. 117, 17635-17642 (2020)
% [3] Perri, S., Molini A., Hedin O. L. & Porporato. Contrasting effects of
%aridity and seasonality on global salinisation, Nat. Geo. (2022)


clear all
close all
tic
Fsize=14;

%% Inputs
Y=1000;% Simulation time [years]
t=365*Y; % Simulation time [days]
dt=1; % time interval [day]

%%% PHYSICAL VARIABLES%%%%%% 
n=0.45; % soil porosity [-]
Zr=300; % Depth of the rooting zone [mm]
Etmax=4.5; %maximum value of ET [mm/day]
s1=0.8; %soil moisture at field capacity
sw=0;% wilting point
beta=0.084; %[l/g]  slope of decrease in ET (or yield) with incresing salinity. This value refer to plants with moderate tolerance. Maas et al.
CT=3.84;% [g/l] salt stress threshold. This value refer to plants with moderate tolerance. Maas et al.
betaLT=0.12; % l/g low tolerance Maas et al.
CTLT=1.92;% g/l low tolerance Maas et al.
betaHT=0.071; % l/g high tolerance Maas et al.
CTHT=6.4;% g/l high tolerance Maas et al.

lambda_P=0.1; % 0.1 frequency of precipitation [1/day]
alpha_p=9; % average precipitation per event[mm] - Semi arid climate
alpha_p_Arid=5.5;% average precipitation per event[mm] - Arid climate

omega_p=2*pi/365; % Frequency of the sinusoidal function for precipitation; Laio 2002. 
phi_p=11/12*2*pi; % Phase shift(maximum of the precipitation in January); Laio 2002.

% lambda_Pt=0.1; % time-dependent frequency of precipitation [1/day]

omega_e=2*pi/365; % Frequency of the sinusoidal function for potential evapotranspiration; Laio 2002. 
phi_e=5/12*2*pi;% Laio 2002 phase shift(maximum of the potential evapotranspiration in July).

%Salt input parameters
in=60*10^(-3);%salt input  g/day*m^2 in coastal area (200 kg/ha/yr)
inCont=6*10^(-3);%salt input  g/day*m^2 in continental area (20 kg/ha/yr)

%% Generate the Rainfal input without seasonality
[RainNS] = RainGenNS(t,s1,n,Zr,lambda_P,1.025*alpha_p); %call the function that generate stochastic precipitation

%% Generate the Rainfal input with seasonality
[Rain,lambda_Pt,alpha_Pt] = RainGenSeasonal(t,s1,n,Zr,lambda_P,omega_p,phi_p,alpha_p,0.25);

%precipitation and PET in phase
[RainIP,lambda_PtIP,alpha_PtIP] = RainGenSeasonal(t,s1,n,Zr,lambda_P,omega_p,phi_e,alpha_p,0.25);


RainCumul=(sum(Rain))/Y*n*Zr*s1; % cumulative rain in the case with seasonal fluctuations
RainNSCumul=sum(RainNS)/Y*n*Zr*s1; % cumulative rain in the case with seasonal fluctuations
%RainCumul shoul be equal to RainNSCumul if all the parameters are
%maintained constant.
 
%% water and salt balance with rainfall input without seasonality
[LeakNS,MNS,sTNS,ThetaNS,ENS,sNS,CNS,CMaxNS,~,~,~,~,~,ETmaxNS] = SoilWb(t,beta,CT,n,Zr,Etmax,RainNS,in,sw,omega_e,phi_e);
[LeakNSLT,MNSLT,sTNSLT,ThetaNSLT,ENSLT,sNSLT,CNSLT,CMaxNSLT,~,~,~,~,~,ETmaxNSLT] = SoilWb(t,betaLT,CTLT,n,Zr,Etmax,RainNS,in,sw,omega_e,phi_e);
[LeakNSHT,MNSHT,sTNSHT,ThetaNSHT,ENSHT,sNSHT,CNSHT,CMaxNSHT,~,~,~,~,~,ETmaxNSHT] = SoilWb(t,betaHT,CTHT,n,Zr,Etmax,RainNS,in,sw,omega_e,phi_e);

%% water and salt balance with rainfall input with seasonality
%precipitation and PET out of phase
[Leak,M,sT,Theta,E,s,C,CMax,~,~,~,ss,~,ETmax] = SoilWb(t,beta,CT,n,Zr,Etmax,Rain,in,sw,omega_e,phi_e);
[LeakLT,MLT,sTLT,ThetaLT,ELT,sLT,CLT,CMaxLT,~,~,~,~,~,ETmaxLT] = SoilWb(t,betaLT,CTLT,n,Zr,Etmax,Rain,in,sw,omega_e,phi_e);
[LeakHT,MHT,sTHT,ThetaHT,EHT,sHT,CHT,CMaxHT,~,~,~,~,~,ETmaxHT] = SoilWb(t,betaHT,CTHT,n,Zr,Etmax,Rain,in,sw,omega_e,phi_e);

%precipitation and PET in phase
[LeakIP,MIP,sTIP,ThetaIP,EIP,sIP,CIP,CMaxIP,~,~,~,ssIP,~,ETmaxIP] = SoilWb(t,beta,CT,n,Zr,Etmax,RainIP,in,sw,omega_e,phi_e);
[LeakIPLT,MIPLT,sTIPLT,ThetaIPLT,EIPLT,sIPLT,CIPLT,CMaxIPLT,~,~,~,~,~,ETmaxIPLT] = SoilWb(t,betaLT,CTLT,n,Zr,Etmax,RainIP,in,sw,omega_e,phi_e);
[LeakIPHT,MIPHT,sTIPHT,ThetaIPHT,EIPHT,sIPHT,CIPHT,CMaxIPHT,~,~,~,~,~,ETmaxIPHT] = SoilWb(t,betaHT,CTHT,n,Zr,Etmax,RainIP,in,sw,omega_e,phi_e);

%% monthly means
[M_YearMonths,~] = MonthlyMeans(t,M);
[s_YearMonths,s_YearMonthsSD] = MonthlyMeans(t,s);
[ss_YearMonths,ss_YearMonthsSD] = MonthlyMeans(t,ss);
[C_YearMonths,~] = MonthlyMeans(t,C);
[E_YearMonths,~] = MonthlyMeans(t,E);
[Rain_YearMonths,~] = MonthlyMeans(t,Rain);
[lambda_Pt_YearMonths,~] = MonthlyMeans(t,lambda_Pt);
[alpha_Pt_YearMonths,~] = MonthlyMeans(t,alpha_Pt);
[PET_YearMonths,~] = MonthlyMeans(t,ETmax);

[sIP_YearMonths,sIP_YearMonthsSD] = MonthlyMeans(t,sIP);
[ssIP_YearMonths,~] = MonthlyMeans(t,ssIP);
[lambda_PtIP_YearMonths,~] = MonthlyMeans(t,lambda_PtIP);
[alpha_PtIP_YearMonths,~] = MonthlyMeans(t,alpha_PtIP);

[M_YearMonthsNS,~] = MonthlyMeans(t,MNS);
[s_YearMonthsNS,s_YearMonthsNSSD] = MonthlyMeans(t,sNS);
[C_YearMonthsNS,~] = MonthlyMeans(t,CNS);
[E_YearMonthsNS,~] = MonthlyMeans(t,ENS);
[Rain_YearMonthsNS,~] = MonthlyMeans(t,RainNS);
[PET_YearMonthsNS,~] = MonthlyMeans(t,ETmaxNS);

%% plots
figure (1)
subplot (1,2,1)
x=[1:1:12];
scatter(x,PET_YearMonths*n*Zr*s1/10,'r')
hold on
scatter(x,alpha_Pt_YearMonths.*lambda_Pt_YearMonths*30/(n*Zr*s1),'b')
hold on
scatter(x,s_YearMonths,'k')
hold on
plot(x,PET_YearMonths*n*Zr*s1/10,'--r','linewidth',1) %PET in cm/day
hold on
plot(x,alpha_Pt_YearMonths.*lambda_Pt_YearMonths*30/(n*Zr*s1),'--b','linewidth',1)
hold on
plot(x,s_YearMonths,'--k','linewidth',1)
hold off
xlim([1 12])
axis square
box on
ylim([0.15 0.5])
xticks([1:1:12])
xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
legend('Potential Evapotranspiration [cm/day]','Effective Rainfall [-]','Relative Soil Moisture [-]')
legend boxoff     
title('PET and Rain out of phase')

subplot (1,2,2)
plot(x,PET_YearMonths*n*Zr*s1/10,'--r','linewidth',1) %PET in cm/day
hold on
scatter(x,PET_YearMonths*n*Zr*s1/10,'r')
hold on
plot(x,alpha_PtIP_YearMonths.*lambda_PtIP_YearMonths*30/(n*Zr*s1),'--b','linewidth',1)
hold on
scatter(x,alpha_PtIP_YearMonths.*lambda_PtIP_YearMonths*30/(n*Zr*s1),'b')
hold on
plot(x,sIP_YearMonths,'--k','linewidth',1)
hold on
scatter(x,sIP_YearMonths,'k')
hold off
xlim([1 12])
ylim([0.15 0.5])
axis square
box on
xticks([1:1:12])
xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
title('PET and Rain in phase')

figure (2)
[f,xi] = ksdensity(s); 
plot(xi,f,'--b','linewidth',2)
hold on
[f,xi] = ksdensity(sNS); 
plot(xi,f,'--k','linewidth',2)
% hold on
% [f,xi] = ksdensity(M_I(21,:)); 
% plot(xi,f,'--r','linewidth',2)
xlim([0 1])
% ylim([0 0.007])
xlabel('Realtive soil moisture [-]') 
ylabel('p(s)') 
legend('With seasonality','Without seasonality')
title('p(s) with and without seasonality')
hold off



figure (3)
subplot (1,2,1)
[f,xi] = ksdensity(MLT); 
plot(xi,f,'r','linewidth',2)
hold on
[f,xi] = ksdensity(MNSLT); 
plot(xi,f,'--r','linewidth',2)
hold on
[f,xi] = ksdensity(MHT); 
plot(xi,f,'g','linewidth',2)
hold on
[f,xi] = ksdensity(MNSHT); 
plot(xi,f,'--g','linewidth',2)
xlim([0 850])
ylim([0 0.008])
xlabel('Salt mass [g/m^2]') 
ylabel('p(m)') 
axis square
box on
legend('With seasonality','Without seasonality')
title('p(m) with and without seasonality')
hold off

subplot (1,2,2)
[f,xi] = ksdensity(CLT); 
plot(xi,f,'k','linewidth',2)
hold on
[f,xi] = ksdensity(CIPLT); 
plot(xi,f,'--k','linewidth',2)
hold on
% xline(CMaxLT(1))
% hold on
[f,xi] = ksdensity(CHT); 
plot(xi,f,'b','linewidth',2)
hold on
[f,xi] = ksdensity(CIPHT); 
plot(xi,f,'--b','linewidth',2)
% hold on
% [f,xi] = ksdensity(M_I(21,:)); 
% plot(xi,f,'--r','linewidth',2)
 xlim([0 21])
% ylim([0 0.007])
xlabel('C [g/L]') 
ylabel('p(C)') 
axis square
box on
legend('\Delta\Phi=\pi Salt-sensitive vegetation','\Delta\Phi=0 Salt-sensitive vegetation',...
    '\Delta\Phi=\pi Salt-tolerant vegetation','\Delta\Phi=0 Salt-tolerant vegetation')
legend boxoff 
title('p(C) with and without seasonality')
hold off

% figure (4)
% [f,xi] = ksdensity(CLT); 
% plot(xi,f,'r','linewidth',2)
% hold on
% [f,xi] = ksdensity(CNSLT); 
% plot(xi,f,'--r','linewidth',2)
% hold on
% [f,xi] = ksdensity(CHT); 
% plot(xi,f,'g','linewidth',2)
% hold on
% [f,xi] = ksdensity(CNSHT); 
% plot(xi,f,'--g','linewidth',2)
% % hold on
% % [f,xi] = ksdensity(M_I(21,:)); 
% % plot(xi,f,'--r','linewidth',2)
%  xlim([0 21])
% % ylim([0 0.007])
% xlabel('Salt concentration [g/L]') 
% ylabel('p(C)') 
% legend('With seasonality','Without seasonality')
% title('p(C) with and without seasonality')
% hold off

% figure (333)

% figure (5)
% [f,xi] = ksdensity(C); 
% plot(xi,f,'--b','linewidth',2)
% hold on
% [f,xi] = ksdensity(CNS); 
% plot(xi,f,'--k','linewidth',2)
% % hold on
% % [f,xi] = ksdensity(M_I(21,:)); 
% % plot(xi,f,'--r','linewidth',2)
%  xlim([0 17])
% % ylim([0 0.007])
% xlabel('Salt concentration [g/L]') 
% ylabel('p(C)') 
% legend('With seasonality','Without seasonality')
% title('p(C) with and without seasonality')
% hold off

toc
return

tic
%% Generate seasonal Rainfall inputs with varying phase shift
phi_e=5/12*2*pi;%
Delta=[1:0.1:12]/12;
phi_pV=Delta*2*pi;
Delta_phi=phi_pV-phi_e;


for j=1:length(phi_pV)
[Rains(:,j),~,~] = RainGenSeasonal(t,s1,n,Zr,lambda_P,omega_p,phi_pV(j),alpha_p,0.25);
[Rains_Arid(:,j),~,~] = RainGenSeasonal(t,s1,n,Zr,lambda_P,omega_p,phi_pV(j),alpha_p_Arid,0.25);
end
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean,C_mean,~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains,in,sw,omega_e,phi_e);

[M_meanDisc,phi_pVDisc,Delta_phi_pVDisc,C_meanDisc] = Mean_DeltaPhi(t,beta,CT,n,Zr,Etmax,Rains,in,sw,omega_e,phi_e,phi_pV);
[M_meanDisc_Arid,phi_pVDisc_Arid,Delta_phi_pVDisc_Arid,C_meanDisc_Arid] = Mean_DeltaPhi(t,beta,CT,n,Zr,Etmax,Rains_Arid,in,sw,omega_e,phi_e,phi_pV);
[M_meanDiscLT,phi_pVDiscLT,Delta_phi_pVDiscLT,C_meanDiscLT] = Mean_DeltaPhi(t,betaLT,CTLT,n,Zr,Etmax,Rains,in,sw,omega_e,phi_e,phi_pV);
[M_meanDiscLT_Arid,phi_pVDiscLT_Arid,Delta_phi_pVDiscLT_Arid,C_meanDiscLT_Arid] = Mean_DeltaPhi(t,betaLT,CTLT,n,Zr,Etmax,Rains_Arid,in,sw,omega_e,phi_e,phi_pV);
[M_meanDiscHT,phi_pVDiscHT,Delta_phi_pVDiscHT,C_meanDiscHT] = Mean_DeltaPhi(t,betaHT,CTHT,n,Zr,Etmax,Rains,in,sw,omega_e,phi_e,phi_pV);
[M_meanDiscHT_Arid,phi_pVDiscHT_Arid,Delta_phi_pVDiscHT_Arid,C_meanDiscHT_Arid] = Mean_DeltaPhi(t,betaHT,CTHT,n,Zr,Etmax,Rains_Arid,in,sw,omega_e,phi_e,phi_pV);


Delta_phi_pVDisc_=[-Delta_phi_pVDisc(end);Delta_phi_pVDisc];
M_meanDiscLT_=[M_meanDiscLT(end);M_meanDiscLT];
M_meanDiscHT_=[M_meanDiscHT(end);M_meanDiscHT];
C_meanDiscLT_=[C_meanDiscLT(end);C_meanDiscLT];
C_meanDiscHT_=[C_meanDiscHT(end);C_meanDiscHT];
Delta_phi_pVDisc_Arid_=[-Delta_phi_pVDisc_Arid(end);Delta_phi_pVDisc_Arid];
M_meanDiscLT_Arid_=[M_meanDiscLT_Arid(end);M_meanDiscLT_Arid];
M_meanDiscHT_Arid_=[M_meanDiscHT_Arid(end);M_meanDiscHT_Arid];
C_meanDiscLT_Arid_=[C_meanDiscLT_Arid(end);C_meanDiscLT_Arid];
C_meanDiscHT_Arid_=[C_meanDiscHT_Arid(end);C_meanDiscHT_Arid];
% 
% figure (6)
% subplot (2,2,1)
% % scatter(Delta_phi_pVDisc,M_meanDisc,'k')
% % hold on
% scatter(Delta_phi_pVDisc_,M_meanDiscLT_,'r')
% hold on
% scatter(Delta_phi_pVDisc_,M_meanDiscHT_,'g')
% hold on
% scatter(Delta_phi_pVDisc_,M_meanDiscLT_,'r')
% hold on
% % plot(Delta_phi_pVDisc,M_meanDisc,'--k','linewidth',1) 
% % hold on
% plot(Delta_phi_pVDisc_,M_meanDiscLT_,'--r','linewidth',1) 
% hold on
% plot(Delta_phi_pVDisc_,M_meanDiscHT_,'--g','linewidth',1) 
% hold on
% % xlim([-2.1 3.7])
% % xlim([-3.14 3.14])
% axis square
% ylim([0 600])
%  set(gca,'XTick',-pi:pi/2:pi) 
%  set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
% xlabel('\Delta\Phi = \Phi_p - \Phi_e') 
% ylabel('Average salt mass [g/m^2]') 
% legend('Salt-tolerant vegetation','Salt-sensitive vegetation')
% legend boxoff     
% title('Semi-arid')
% 
% subplot (2,2,3)
% % scatter(Delta_phi_pVDisc_Arid,M_meanDisc_Arid,'k')
% % hold on
% scatter(Delta_phi_pVDisc_Arid_,M_meanDiscLT_Arid_,'r')
% hold on
% scatter(Delta_phi_pVDisc_Arid_,M_meanDiscHT_Arid_,'g')
% hold on
% % plot(Delta_phi_pVDisc_Arid,M_meanDisc_Arid,'--k','linewidth',1) 
% % hold on
% plot(Delta_phi_pVDisc_Arid_,M_meanDiscLT_Arid_,'--r','linewidth',1) 
% hold on
% plot(Delta_phi_pVDisc_Arid_,M_meanDiscHT_Arid_,'--g','linewidth',1) 
% hold on
% axis square
% ylim([150 1200])
% % xlim([-3.14 3.14])
%  set(gca,'XTick',-pi:pi/2:pi) 
%  set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
% xlabel('\Delta\Phi = \Phi_p - \Phi_e') 
% ylabel('Average salt mass [g/m^2]') 
% title('Arid')
% 
% 
% 
figure (7)
subplot (2,2,1)
% scatter(Delta_phi_pVDisc,M_meanDisc,'k')
% hold on
scatter(Delta_phi_pVDisc_,C_meanDiscLT_,'r')
hold on
scatter(Delta_phi_pVDisc_,C_meanDiscHT_,'g')
hold on
scatter(Delta_phi_pVDisc_,C_meanDiscLT_,'r')
hold on
% plot(Delta_phi_pVDisc,M_meanDisc,'--k','linewidth',1) 
% hold on
plot(Delta_phi_pVDisc_,C_meanDiscLT_,'--r','linewidth',1) 
hold on
plot(Delta_phi_pVDisc_,C_meanDiscHT_,'--g','linewidth',1) 
hold on
% xlim([-2.1 3.7])
% xlim([-3.14 3.14])
axis square
ylim([3 15])
 set(gca,'XTick',-pi:pi/2:pi) 
 set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
xlabel('\Delta\Phi = \Phi_p - \Phi_e') 
ylabel('C [g/L]') 
legend('Salt-tolerant vegetation','Salt-sensitive vegetation')
legend boxoff     
title('Semi-arid')

subplot (2,2,3)
% scatter(Delta_phi_pVDisc_Arid,M_meanDisc_Arid,'k')
% hold on
scatter(Delta_phi_pVDisc_Arid_,C_meanDiscLT_Arid_,'r')
hold on
scatter(Delta_phi_pVDisc_Arid_,C_meanDiscHT_Arid_,'g')
hold on
% plot(Delta_phi_pVDisc_Arid,M_meanDisc_Arid,'--k','linewidth',1) 
% hold on
plot(Delta_phi_pVDisc_Arid_,C_meanDiscLT_Arid_,'--r','linewidth',1) 
hold on
plot(Delta_phi_pVDisc_Arid_,C_meanDiscHT_Arid_,'--g','linewidth',1) 
hold on
axis square
ylim([5 20])
% xlim([-3.14 3.14])
 set(gca,'XTick',-pi:pi/2:pi) 
 set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
xlabel('\Delta\Phi = \Phi_p - \Phi_e') 
ylabel('C [g/L]') 
title('Arid')

%% Generate seasonal Rainfall inputs with varying Seasonal amplitute

A=[0:0.05:0.55];


for j=1:length(A)
[Rains(:,j),~,~] = RainGenSeasonal(t,s1,n,Zr,lambda_P,omega_p,phi_p,alpha_p,A(j));
[Rains_Arid(:,j),~,~] = RainGenSeasonal(t,s1,n,Zr,lambda_P,omega_p,phi_p,alpha_p_Arid,A(j));
end

[M_meanDisc_A,C_meanDisc_A] = Mean_DeltaPhi_A(t,beta,CT,n,Zr,Etmax,Rains,in,sw,omega_e,phi_e);
[M_meanDisc_Arid_A,C_meanDisc_Arid_A] = Mean_DeltaPhi_A(t,beta,CT,n,Zr,Etmax,Rains_Arid,in,sw,omega_e,phi_e);
[M_meanDiscLT_A,C_meanDiscLT_A] = Mean_DeltaPhi_A(t,betaLT,CTLT,n,Zr,Etmax,Rains,in,sw,omega_e,phi_e);
[M_meanDiscLT_Arid_A,C_meanDiscLT_Arid_A] = Mean_DeltaPhi_A(t,betaLT,CTLT,n,Zr,Etmax,Rains_Arid,in,sw,omega_e,phi_e);
[M_meanDiscHT_A,C_meanDiscHT_A] = Mean_DeltaPhi_A(t,betaHT,CTHT,n,Zr,Etmax,Rains,in,sw,omega_e,phi_e);
[M_meanDiscHT_Arid_A,C_meanDiscHT_Arid_A] = Mean_DeltaPhi_A(t,betaHT,CTHT,n,Zr,Etmax,Rains_Arid,in,sw,omega_e,phi_e);


% figure (6)
% subplot (2,2,2)
% % scatter(Delta_phi_pVDisc,M_meanDisc,'k')
% % hold on
% scatter(A,M_meanDiscLT_A,'r')
% hold on
% scatter(A,M_meanDiscHT_A,'g')
% hold on
% % plot(Delta_phi_pVDisc,M_meanDisc,'--k','linewidth',1) 
% % hold on
% plot(A,M_meanDiscLT_A,'--r','linewidth',1) 
% hold on
% plot(A,M_meanDiscHT_A,'--g','linewidth',1) 
% hold on
% % xlim([-2.1 3.7])
% xlim([0 0.5])
% axis square
%  ylim([0 600])
% %  set(gca,'XTick',-pi:pi/2:pi) 
% %  set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
% xlabel('Seasonality amplitude') 
% ylabel('Average salt mass [g/m^2]') 
% legend('Salt-tolerant vegetation','Salt-sensitive vegetation')
% legend boxoff     
% title('Semi-arid')
% 
% subplot (2,2,4)
% % scatter(Delta_phi_pVDisc_Arid,M_meanDisc_Arid,'k')
% % hold on
% scatter(A,M_meanDiscLT_Arid_A,'r')
% hold on
% scatter(A,M_meanDiscHT_Arid_A,'g')
% hold on
% % plot(Delta_phi_pVDisc_Arid,M_meanDisc_Arid,'--k','linewidth',1) 
% % hold on
% plot(A,M_meanDiscLT_Arid_A,'--r','linewidth',1) 
% hold on
% plot(A,M_meanDiscHT_Arid_A,'--g','linewidth',1) 
% hold on
% axis square
% ylim([150 1200])
% xlim([0 0.5])
% %  set(gca,'XTick',-pi:pi/2:pi) 
% %  set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
% xlabel('Seasonality amplitude') 
% ylabel('Average salt mass [g/m^2]') 
% title('Arid')


figure (7)
subplot (2,2,2)
% scatter(Delta_phi_pVDisc,M_meanDisc,'k')
% hold on
scatter(A,C_meanDiscLT_A,'r')
hold on
scatter(A,C_meanDiscHT_A,'g')
hold on
% plot(Delta_phi_pVDisc,M_meanDisc,'--k','linewidth',1) 
% hold on
plot(A,C_meanDiscLT_A,'--r','linewidth',1) 
hold on
plot(A,C_meanDiscHT_A,'--g','linewidth',1) 
hold on
% xlim([-2.1 3.7])
xlim([0 0.5])
axis square
 ylim([3 15])
%  set(gca,'XTick',-pi:pi/2:pi) 
%  set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
xlabel('\delta_p \equiv \delta_{et}') 
ylabel('C [g/L]') 
legend('Salt-tolerant vegetation','Salt-sensitive vegetation')
legend boxoff     
title('Semi-arid')

subplot (2,2,4)
% scatter(Delta_phi_pVDisc_Arid,M_meanDisc_Arid,'k')
% hold on
scatter(A,C_meanDiscLT_Arid_A,'r')
hold on
scatter(A,C_meanDiscHT_Arid_A,'g')
hold on
% plot(Delta_phi_pVDisc_Arid,M_meanDisc_Arid,'--k','linewidth',1) 
% hold on
plot(A,C_meanDiscLT_Arid_A,'--r','linewidth',1) 
hold on
plot(A,C_meanDiscHT_Arid_A,'--g','linewidth',1) 
hold on
axis square
ylim([5 20])
xlim([0 0.5])
%  set(gca,'XTick',-pi:pi/2:pi) 
%  set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
xlabel('\delta_p \equiv \delta_{et}') 
ylabel('C [g/L]') 
title('Arid')

%% Generate seasonal Rainfall inputs with varying Dryness
tic

Y=4000;% Simulation time [years]
t=365*Y; % Simulation time [days]
% A=[0:0.05:0.55];

in=6*10^(-3);%salt input  g/day*m^2 in coastal area (200 kg/ha/yr)



alpha_p_V=[5:0.1:25];
lambda_P_V=[0.05:0.0015:.35];
DI=Etmax./(lambda_P_V.*alpha_p_V);
Rains=zeros(t,1);
Rains_A=zeros(t,1);
Rains_V=zeros(t,1);
in_V=[0:20^(-2):600*10^(-3)];


beta_V=[0.05:0.001:0.2];
CT_V=[1:0.05:10];
lambda_P6=Etmax./(6.*alpha_p);
lambda_P1=Etmax./(1.*alpha_p);
lambda_P3=Etmax./(3.*alpha_p);

[Rain6,~,~] = RainGenSeasonal(t,s1,n,Zr,lambda_P6,omega_p,phi_p,alpha_p,0.25);
[Rain1,~,~] = RainGenSeasonal(t,s1,n,Zr,lambda_P1,omega_p,phi_p,alpha_p,0.25);
[Rain3,~,~] = RainGenSeasonal(t,s1,n,Zr,lambda_P3,omega_p,phi_p,alpha_p,0.25);

[Rain,~,~] = RainGenSeasonal(t,s1,n,Zr,lambda_P,omega_p,phi_p,alpha_p,0.25);

for j=1:length(beta_V)
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean_beta_V1(j),C_mean_beta_V1(j),~,~] = SoilWb_Mean(t,beta_V(j),CT,n,Zr,Etmax,Rain1,in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean_beta_V3(j),C_mean_beta_V3(j),~,~] = SoilWb_Mean(t,beta_V(j),CT,n,Zr,Etmax,Rain3,in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean_beta_V6(j),C_mean_beta_V6(j),~,~] = SoilWb_Mean(t,beta_V(j),CT,n,Zr,Etmax,Rain6,in,sw,omega_e,phi_e);
end

for j=1:length(CT_V)
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean_CT_V1(j),C_mean_CT_V1(j),~,~] = SoilWb_Mean(t,beta,CT_V(j),n,Zr,Etmax,Rain1,in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean_CT_V3(j),C_mean_CT_V3(j),~,~] = SoilWb_Mean(t,beta,CT_V(j),n,Zr,Etmax,Rain3,in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean_CT_V6(j),C_mean_CT_V6(j),~,~] = SoilWb_Mean(t,beta,CT_V(j),n,Zr,Etmax,Rain6,in,sw,omega_e,phi_e);
end

figure (10)
subplot (1,2,1)
k=1;
plot(movmean(beta_V,k),movmean(C_mean_beta_V6,k),'k','linewidth',2) 
hold on
plot(movmean(beta_V,k),movmean(C_mean_beta_V3,k),'b','linewidth',2) 
hold on
plot(movmean(beta_V,k),movmean(C_mean_beta_V1,k),'g','linewidth',2) 
hold off
ylim([0 13])
axis square
xlabel('\beta [L/g]') 
ylabel('C [g/L]') 
legend('DI=6','DI=3','DI=1')
legend boxoff    

subplot (1,2,2)
plot(movmean(CT_V,k),movmean(C_mean_CT_V6,k),'k','linewidth',2) 
hold on
plot(movmean(CT_V,k),movmean(C_mean_CT_V3,k),'b','linewidth',2) 
hold on
plot(movmean(CT_V,k),movmean(C_mean_CT_V1,k),'g','linewidth',2) 
hold off
axis square
xlabel('C_T [g/L]') 
ylabel('C [g/L]') 
xlim([min(CT_V) max(CT_V)])
ylim([0 13])
% legend('DI=6','DI=3','DI=1')
% legend boxoff 





for j=1:length(in_V)
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean_in6(j),C_mean_in6(j),~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rain6,in_V(j),sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean_in3(j),C_mean_in3(j),~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rain3,in_V(j),sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean_in1(j),C_mean_in1(j),~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rain1,in_V(j),sw,omega_e,phi_e);

% [Rains_Arid(:,j),~,~] = RainGenSeasonal(t,s1,n,Zr,lambda_P_V,omega_p,phi_p,alpha_p_Arid,.25);
end
% 
figure (11)
% subplot (1,3,1)
k=1;
plot(movmean(in_V,k),movmean(C_mean_in6,k),'k','linewidth',2) 
hold on
plot(movmean(in_V,k),movmean(C_mean_in3,k),'b','linewidth',2) 
hold on
plot(movmean(in_V,k),movmean(C_mean_in1,k),'g','linewidth',2) 
hold off
% xlim([-2.1 3.7])
% xlim([0.5 10])
axis square
% ylim([0 20])
%  set(gca,'XTick',-pi:pi/2:pi) 
%  set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
xlabel('\Upsilon [g/(m^2 day)]') 
ylabel('C [g/L]') 
legend('DI=6','DI=3','DI=1')
legend boxoff     
% title('Semi-arid')



for j=1:length(DI)
[Rains(:,j),~,~] = RainGenSeasonal(t,s1,n,Zr,lambda_P_V(j),omega_p,phi_p,alpha_p_V(j),.25);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_meanHT(j),C_meanHT(j),~,~] = SoilWb_Mean(t,betaHT,CTHT,n,Zr,Etmax,Rains(:,j),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean(j),C_mean(j),~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains(:,j),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_meanLT(j),C_meanLT(j),~,~] = SoilWb_Mean(t,betaLT,CTLT,n,Zr,Etmax,Rains(:,j),in,sw,omega_e,phi_e);
end




A=[0:0.00275:0.55];


for j=1:length(A)
[Rains_A(:,j),~,~] = RainGenSeasonal(t,s1,n,Zr,lambda_P,omega_p,phi_p,alpha_p,A(j));
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean_A(j),C_mean_A(j),~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains_A(:,j),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_meanLT_A(j),C_meanLT_A(j),~,~] = SoilWb_Mean(t,betaLT,CTLT,n,Zr,Etmax,Rains_A(:,j),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_meanHT_A(j),C_meanHT_A(j),~,~] = SoilWb_Mean(t,betaHT,CTHT,n,Zr,Etmax,Rains_A(:,j),in,sw,omega_e,phi_e);

end


phi_e=6/12*2*pi;%
Delta=[0:0.0600:12]/12;
phi_pV=Delta*2*pi;
Delta_phi=phi_pV-phi_e;


for j=1:length(phi_pV)
[Rains_V(:,j),~,~] = RainGenSeasonal(t,s1,n,Zr,lambda_P,omega_p,phi_pV(j),alpha_p,0.25);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean_V(j),C_mean_V(j),~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains_V(:,j),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_meanLT_V(j),C_meanLT_V(j),~,~] = SoilWb_Mean(t,betaLT,CTLT,n,Zr,Etmax,Rains_V(:,j),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_meanHT_V(j),C_meanHT_V(j),~,~] = SoilWb_Mean(t,betaHT,CTHT,n,Zr,Etmax,Rains_V(:,j),in,sw,omega_e,phi_e);

end


figure (12)
subplot (1,3,1)
k=20;
plot(movmean(DI,k),movmean(C_meanHT,k),'k','linewidth',1) 
hold on
plot(movmean(DI,k),movmean(C_meanLT,k),'b','linewidth',1) 
hold on
plot(movmean(DI,k),movmean(C_mean,k),'g','linewidth',1) 
hold off
% xlim([-2.1 3.7])
xlim([0.5 10])
axis square
ylim([0 20])
xlabel('DI') 
ylabel('C [g/L]') 
legend('Salt-tolerant','Salt-sensitive')
legend boxoff     


subplot (1,3,2)
plot(movmean(A,k),movmean(C_meanHT_A,k),'k','linewidth',1) 
hold on
plot(movmean(A,k),movmean(C_meanLT_A,k),'b','linewidth',1) 
hold on
plot(movmean(A,k),movmean(C_mean_A,k),'g','linewidth',1) 
hold off
% xlim([-2.1 3.7])
xlim([0 .5])
axis square
ylim([0 20])

xlabel('\delta') 
ylabel('C [g/L]') 
  

subplot (1,3,3)
plot(movmean(Delta_phi,k),movmean(C_meanHT_V,k),'k','linewidth',1) 
hold on
plot(movmean(Delta_phi,k),movmean(C_meanLT_V,k),'b','linewidth',1) 
hold on
plot(movmean(Delta_phi,k),movmean(C_mean_V,k),'g','linewidth',1) 
hold off
% xlim([-2.1 3.7])
xlim([-pi pi])
axis square
ylim([0 20])
 set(gca,'XTick',-pi:pi/2:pi) 
 set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
xlabel('\Delta \phi') 
ylabel('C [g/L]') 


toc


