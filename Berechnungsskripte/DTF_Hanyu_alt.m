function [DTF] = DTF_Hanyu_alt(Impulsantwort,L_5dB,L_vardB,Sw)
%% Zusammenfassung
% Berechnung des DTF (degree of time fluctuation) nach Hanyu 2014
% Autor: Linus Staubach; Version 1; Erstellung 2022
% 
%% Nähere Beschreibung
% Syntax: 
% [DTF] = DTF_Hanyu_alt(Impulsantwort,L_5dB,L_vardB,Sw)
% 
% Input:
% Impulsantwort     zu Analysierende Impulsantwort;
% L_5dB             Anfang Analysebereich [Sample];
% L_vardB           Ende Analysebereich [Sample];
% Sw                Schrittweite der Schwellenwertberechnung (verändert die Berechnungspräzision)
% 
% Output:
% DTF               Degree of Time Series Fluctuation [dB/m]
%% Berechnung
g_t = Impulsantwort./sqrt(flipud(cumsum(flipud(Impulsantwort.^2))));
g_t_ber = g_t([L_5dB:L_vardB],:);

h_t = (1./sqrt(mean(g_t_ber.^2))).*g_t_ber;
R_total = sum(h_t.^2);

z_k = zeros(Sw,1);
for tr = 1:Sw
    z_k(tr,:) = sum(h_t(h_t.^2 > ((max(h_t.^2)*tr)/Sw)).^2)./R_total;
end
DTF = (length(z_k(z_k >= 0.01))/Sw)*max(h_t.^2);

end

