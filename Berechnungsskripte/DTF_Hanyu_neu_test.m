function [DTF, A_t] = DTF_Hanyu_neu_test(Impulsantwort,L_5dB,L_vardB,Sw)
%% Zusammenfassung
% Berechnung des DTF (degree of time fluctuation) nach Hanyu 2018
% Autor: Linus Staubach; Version 2; Erstellung 2022
% 
%% Nähere Beschreibung
% Syntax: 
% [DTF, A_t] = DTF_Hanyu_neu_test(Impulsantwort,L_5dB,L_vardB,Sw)
% 
% Input:
% Impulsantwort     zu Analysierende Impulsantwort;
% L_5dB             Anfang Analysebereich [Sample];
% L_vardB           Ende Analysebereich [Sample];
% Sw                Schrittweite der Schwellenwertberechnung (verändert die Berechnungspräzision)
% 
% Output:
% DTF               Degree of Time Series Fluctuation [dB/m];
% A_t               Momentane Abklingrate
%% Berechnung
                g_t = Impulsantwort./sqrt(flipud(cumsum(flipud(Impulsantwort.^2))));
                g_t_ber = g_t([L_5dB:L_vardB],:);

                EDC2 = 20*log10(sqrt(flipud(cumsum(flipud(Impulsantwort.^2))))/(2*10^(-5)));
                A_t = zeros(L_vardB-L_5dB+1,1);
                
%                     for t2 = L_5dB:L_vardB
%                         A_t(t2-L_5dB+1) = (log(10))/(find(EDC2 < (EDC2(t2)-5),1)-find(EDC2 < (EDC2(t2)+5),1)); 
%                     end
                    %---- Optimierte Version ----
                    for t2 = 1:(L_vardB-L_5dB+1)
                        %A_t(t2) = log(10)/(find(EDC2 <= (EDC2(t2)-10),1)-t2);
                        A_t(t2) = log(10)/(length(EDC2(EDC2 >= (EDC2(t2)-10),1))-t2);
                    end
                    %----
                h_tv2 = (g_t_ber./sqrt(A_t));
                h_t2v2 = h_tv2.^2;
                R_total2 = sum(h_t2v2);

                z_k2 = zeros(Sw,1);
                    for tr2 = 1:Sw
                        z_k2(tr2,:) = sum(h_t2v2(h_t2v2 > (max((h_t2v2)*(tr2/Sw)))))./R_total2;
                    end
                DTF = (length(z_k2(z_k2 >= 0.01))/Sw)*max(h_t2v2);

end