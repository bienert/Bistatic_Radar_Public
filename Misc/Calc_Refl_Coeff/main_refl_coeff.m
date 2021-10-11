%Nicole Bienert
%4-15-2021

%Calculate the reflection coeffecient as a function of angle for common
%basal materials


clc
clear all
close all

addpath('..\..\Processing_Scripts\functions');

addpath('../../Processing_Scripts/functions');

    
%ice properties 
er_ice = 3.18; %real part of permittivity
tan_ice=0.0062; %loss tangent
%basal properties
er_bed=[2.7,2.8,3.43,6.6,18,80]; %real part of permittivity
tan_bed=[0.022,0.035,0.05,0.41,0.82,0.002]; %loss tangent
mytitle={'Frozen Bedrock';'Frozen Till';'Marine Ice';'Unfrozen Bedrock';'Unfrozen Till';'Fresh Water'}
thetaI=linspace(0,pi/2,1001); %incidence angle

figure()
for i=1:length(er_bed)
    refParE=calcReflCoeff_spol(er_ice,tan_ice,er_bed(i),tan_bed(i),thetaI);
    refParH=calcReflCoeff_ppol(er_ice,tan_ice,er_bed(i),tan_bed(i),thetaI);
    
    %remove cases where thetaT >90deg
    thetaT=sqrt(1-er_ice./er_bed(i).*(sin(thetaI)).^2);
    refParE(thetaT>=pi/2)=0;
    refParH(thetaT>=pi/2)=0;
    
    subplot(2,3,i)
    plot(thetaI*180/pi,abs(refParE),'LineWidth',1.5);
    hold on
    plot(thetaI*180/pi,abs(refParH),'--','LineWidth',1.5);
    xlim([0 90])
    hTitle=title(string(mytitle(i)))
    hLegend=legend('$\vec{\bf{E}}$ Parallel to Interface (S-pol)','$\vec{\bf{H}}$ Parallel to Interface (P-pol)')
    hYlabel=ylabel('Reflection Coeffecient Magnitude');
    hXlabel=xlabel('Angle (Degrees)')
    Aesthetics_Script
end
