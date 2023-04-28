function [Volumen] = VolModellhallraum(Zustand)
%% Zusammenfassung
% Berechnet das Volumen des Modellhallraumes pro Diffusorbelegungszustand
% Autor: Linus Staubach; Version 1; Erstellung 2022
% 
%% Nähere Beschreibung
% Syntax: 
% [Volumen] = VolModellhallraum(Zustand)
% 
% Input:
% Zustand           Belegungszustand mit Diffusorelementen;
% 
% Output:
% Volumen           Volumen des Modellhallraumes [m³];
%% Berechnung
VolAnfang = 0.3175875;
VolKl = 0.00019563;
VolGr = 0.00083148;
ZustandsVerlauf = [2,2,2,2,1,1,1,1,2,2,1,2,1,2,2,1,1,1];


if Zustand ==0
    Volumen = VolAnfang;
else
    VKorr = 0;
    for i = 1:Zustand
        switch ZustandsVerlauf(i)
            case 1
            VKorr = VKorr + VolKl;
            case 2 
            VKorr = VKorr + VolGr;
        end
    end
    Volumen = VolAnfang-VKorr;
end

end

