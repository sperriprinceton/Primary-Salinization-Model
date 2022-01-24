%This function calculates the monthly means of a time series

function [YearMonths,YearMonthsSD] = MonthlyMeans(t,Data)
year=2010; % for example
start=datetime(year,1,1); 
m=month(start+caldays(1:t)-1); %

Jan=nanmean(Data(m==1)); 
Feb=nanmean(Data(m==2)); % february is the second month, ie m==2
Mar=nanmean(Data(m==3)); 
Apr=nanmean(Data(m==4)); % february is the second month, ie m==2
May=nanmean(Data(m==5)); 
Jun=nanmean(Data(m==6)); % february is the second month, ie m==2
Jul=nanmean(Data(m==7)); 
Aug=nanmean(Data(m==8)); % february is the second month, ie m==2
Sep=nanmean(Data(m==9)); 
Oct=nanmean(Data(m==10)); 
Nov=nanmean(Data(m==11)); 
Dec=nanmean(Data(m==12)); 

JanSD=nanstd(Data(m==1)); 
FebSD=nanstd(Data(m==2)); % february is the second month, ie m==2
MarSD=nanstd(Data(m==3)); 
AprSD=nanstd(Data(m==4)); % february is the second month, ie m==2
MaySD=nanstd(Data(m==5)); 
JunSD=nanstd(Data(m==6)); % february is the second month, ie m==2
JulSD=nanstd(Data(m==7)); 
AugSD=nanstd(Data(m==8)); % february is the second month, ie m==2
SepSD=nanstd(Data(m==9)); 
OctSD=nanstd(Data(m==10)); 
NovSD=nanstd(Data(m==11)); 
DecSD=nanstd(Data(m==12)); 
YearMonths=[Jan;Feb;Mar;Apr;May;Jun;Jul;Aug;Sep;Oct;Nov;Dec];
YearMonthsSD=[JanSD;FebSD;MarSD;AprSD;MaySD;JunSD;JulSD;AugSD;SepSD;OctSD;NovSD;DecSD];
end
% end