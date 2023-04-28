function [T,Reg,Fehler] = RT_Regression(Impulsantwort,fs,L_5dB,L_vardB)
%% Zusammenfassung
% Nachhallzeitbestimmung mittels linearer Regression
% Autor: Linus Staubach; Version 1; Erstellung 2022
% 
%% Nähere Beschreibung
% Syntax: 
% [T,Reg,Fehler] = RT_Regression(Impulsantwort,fs,L_5dB,L_vardB)
% 
% Input:
% Impulsantwort     zu Analysierende Impulsantwort;
% fs                Samplingfrequenz [Hz];
% L_5dB             Anfang Analysebereich [Sample];
% L_vardB           Ende Analysebereich [Sample]
% 
% Output:
% T                 Nachhallzeit [s];
% Reg               Regressionsgerade;
% Fehler            Standardabweichung während der Nachhallzeitbestimmung
% [Sample]
%% Berechnung
EDC2 = 20*log10(sqrt(flipud(cumsum(flipud(Impulsantwort.^2))))/(2*10^(-5)));
EDC3 = rot90(EDC2);

L_zeit = 0:(1/fs):(length(EDC2)-1)/fs;
[p,S] = polyfit( L_zeit(:,L_5dB:L_vardB),EDC3(:,L_5dB:L_vardB),1);
linreg = p(2)+L_zeit'.*p(1); 
T = (60/abs(p(1)));
[Reg,Fehler] = polyval(p,L_zeit,S);

end

