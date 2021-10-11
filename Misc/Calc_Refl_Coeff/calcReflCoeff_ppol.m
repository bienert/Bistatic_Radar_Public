%%Nicole Bienert
%4-15-2021
%Calculate the reflection coeffecient for a p-polarized wave. We assume
%low loss, non-magnetic medium

%Note that the loss tangent is defined as:
%          conductivity + omega * epsilon''
% tan(d)= --------------------------------
%                 omega * epsilon'


%inputs: 
% er_ice : The real part of the permittivity of the ice. Typically 3.18
% tan_ice : The loss tangent of ice, typically around 0.006
% er_bed : The real part of the permittivity of the bed
% tan_bed : The loss tangent of the bed

function refParH=calcReflCoeff_ppol(er_ice,tan_ice,er_bed,tan_bed,thetaI)

eta_bed=1./sqrt(er_bed*(1-j*tan_bed));
eta_ice=1./sqrt(er_ice*(1-j*tan_ice));
cosThetaT=sqrt(1-(er_ice*(j*tan_ice-1))./(er_bed*(j*tan_bed-1)).*(sin(thetaI)).^2);
refParH = (eta_bed.*cosThetaT-eta_ice.*cos(thetaI))./(eta_bed.*cosThetaT+eta_ice.*cos(thetaI)); %E Parallel to Interface


