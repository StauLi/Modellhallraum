function [T] = RT_Hanyu_DCIR(Impulsantwort,L_5dB,L_vardB,fs)
%% Zusammenfassung
% Nachhallzeitbestimmung mittels der hallkorrigierten Impulsantwort nach
% Hanyu 2014
% Autor: Linus Staubach; Version 1; Erstellung 2022
% 
%% Nähere Beschreibung
% Syntax: 
% [T] = RT_Hanyu_DCIR(Impulsantwort,L_5dB,L_vardB,fs)
% 
% Input:
% Impulsantwort     zu Analysierende Impulsantwort;
% L_5dB             Anfang Analysebereich [Sample];
% L_vardB           Ende Analysebereich [Sample];
% fs                Samplingfrequenz [Hz]
% 
% Output:
% T                 Nachhallzeit [s]
%% Berechnung
g_t2_T = (Impulsantwort.^2)./(flipud(cumsum(flipud(Impulsantwort.^2))));
g_t2_T = g_t2_T([L_5dB:L_vardB],:);
T = 13.82/(mean(g_t2_T)*fs);
end

