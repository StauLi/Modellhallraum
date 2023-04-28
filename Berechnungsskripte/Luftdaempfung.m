function [alpha] = Luftdaempfung(Temp, H, pa, f)
%% Zusammenfassung
% Berechnet die Dämpfung von Schall bei der Ausbreitung im Freien nach ISO 9613-1.
% Autor: Linus Staubach; Version 1; Erstellung 2022
% 
%% Nähere Beschreibung
% Syntax: 
% [alpha] = Luftdaempfung(Temp, H, pa, f)
% 
% Input:
% Temp          Lufttemperatur [°C];
% H             relative Luftfeuchte [%];
% pa            Luftdruck [kPa];
% f             Betrachtungsfrequenz [Hz]
% 
% Output:
% alpha         Luftdämpfung alpha [dB/m]
%% Berechnung
T = Temp+273.15;    % [K]
pr = 101.325;       % [kPa]
T0 = 293.15;        % [K]

h = (H.*10.^(-6.8346.*((273.16./T).^1.261)+4.6151))./(pa./pr);

fro = (pa./pr).*(24+(4.04.*h.*(10^4).*((0.02+h)./(0.391+h))));
frn = (pa./pr).*((T./T0).^(-0.5)).*(9+280.*h.*exp(-4.17.*(((T./T0).^(-1/3))-1)));

alpha = 8.686.*(f.^2).*((1.84.*(10.^(-11)).*((pa./pr).^-1).*((T./T0).^0.5))+((T./T0).^(-5/2)).*(0.01275.*(exp(-2239.1./T)).*((fro+((f.^2)./fro)).^-1))+(0.1068.*(exp(-3352./T)).*((frn+((f.^2)./frn)).^(-1))));

end