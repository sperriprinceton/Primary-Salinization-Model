
function[M_meanDisc,phi_pVDisc,Delta_phi_pVDisc,C_meanDisc] = Mean_DeltaPhi(t,beta,CT,n,Zr,Etmax,Rains,in,sw,omega_e,phi_e,phi_pV)



[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean1,C_mean1,~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains(:,1),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean11,C_mean11,~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains(:,11),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean21,C_mean21,~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains(:,21),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean31,C_mean31,~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains(:,31),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean41,C_mean41,~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains(:,41),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean51,C_mean51,~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains(:,51),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean61,C_mean61,~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains(:,61),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean71,C_mean71,~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains(:,71),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean81,C_mean81,~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains(:,81),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean91,C_mean91,~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains(:,91),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean101,C_mean101,~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains(:,101),in,sw,omega_e,phi_e);
[~,~,~,~,~,~,~,~,~,~,~,~,~,M_mean111,C_mean111,~,~] = SoilWb_Mean(t,beta,CT,n,Zr,Etmax,Rains(:,111),in,sw,omega_e,phi_e);

M_meanDisc=[M_mean1;M_mean11;M_mean21;M_mean31;M_mean41;M_mean51;M_mean61;M_mean71;M_mean81;M_mean91;M_mean101;M_mean111];
C_meanDisc=[C_mean1;C_mean11;C_mean21;C_mean31;C_mean41;C_mean51;C_mean61;C_mean71;C_mean81;C_mean91;C_mean101;C_mean111];

phi_pVDisc=[phi_pV(1);phi_pV(11);phi_pV(21);phi_pV(31);phi_pV(41);phi_pV(51);phi_pV(61);phi_pV(71);phi_pV(81);phi_pV(91);phi_pV(101);phi_pV(111)];
Delta_phi_pVDisc=phi_pVDisc-phi_e;

end