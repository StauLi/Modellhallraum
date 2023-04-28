%% Messskript zur Qualifizierung von Diffusitätsmessparametern (Gesamtspektrum)(2022/23)
% Autor: Linus Staubach
% Erstellungsdatum: Dezember 2022
% Version: V6.1, Analyse des gesamten Frequenzbereiches

%% Erläuterung
% Dieses Messskript zur Qualifizierung von Diffusitätsparametern wurde im
% Rahmen einer Bachelorarbeit erstellt. Es dient dazu, drei Parameter auf
% ihre Aussagefähigkeit im Zusammenhang mit Diffusitätsänderungen zu
% qualifizieren. Dabei werden der Degree of Time Series Fluctuation (DTF),
% die Standardabweichung der Nachhallzeit und die Standardabweichung der
% Schallenergie betrachtet. Ebenfalls wurden mehrere Varianten der
% Parameter berücksichtigt. In diesem Skript findet die Auswertung im
% gesamten Frequenzbereich der aufgenommenen Impulsantworten statt.

%% Variablendeklaration
clear all
close all

Z = 18;         %Zustand (Diffusoranzahl)
C = 1; Ca = 1;  %C: Raumzustand; Ca: Anfangsraumzustand (Absorberbelegung)
S = 2;          %Lautsprecherpositionen
R = 6;          %Mikrofonpositionen

Ab = 20;        %Analysebereich in dB
Sw = 25;        %Schrittweite des DTF-Thresholds

DTFvX = 0;      %Berechnungsvariante DTF (0=aus; 1=alt(schneller, evtl. ungenauer); 2=neu)
TvX = 0;        %Berechnungsvariante T (0=aus; 1=Regression; 2=aus dcIR; 3=aus A(t)(nur mit DTFvX=2 berechenbar))
KlKorrvX = 0;   %Klimakorrektur T (0=aus; 1=an)(nicht implementiert)
AvX = 0;        %Aeq (0=aus; 1=an)

load('Klimadaten.mat')

V = 0.318;      %Volumen des diffusorunbesetzten Raumes

fs = 96000;     %Samplingfrequenz

%Grenzfreqzenzen der Filterung der Impulsantwort
OGF = 31500;    %Obere Grenzfrequenz [Hz]
UGF = 315;      %Untere Grenzfrequenz [Hz]

%% Ergebnisvariableninitialisierung
DTF_ges = zeros((R*S*(C-Ca+1))+2,Z);
T_ges = zeros((R*S*(C-Ca+1))+3,Z);
p_ges = zeros((R*S*(C-Ca+1))+3,Z);
A_ges = zeros(3,Z);

[b,a] = butter(1,[UGF OGF]/(fs/2),'bandpass');

%% Berechnung
for x1 = 0:Z                    %Anfang der Zustandsschleife der Diffusoren
    for x2 = Ca:C               %Anfang der Raumzustandsschleife
        for x3 = 1:S            %Anfang der Lautsprecherpositionsschleife
            for x4 = 1:R        %Anfang der Mikrofonpositionsschleife
                
%[Impulsantwort,fs] = audioread('c1s1r3_ir_1.wav');
[lisi, Impulsantwort, fs] = MBBMwavread(append("C:\...\Impulsantworten\.Z",num2str(x1),"_roh\c",num2str(x2),"s",num2str(x3),"r",num2str(x4),"_ir.wav"));
Impulsantwort = rot90(Impulsantwort,3)*lisi.peakAmplitude;

%Filterung 
Impulsantwort = filter(b,a,Impulsantwort);

i=1;
                                %Berechnung des Analysebereichs 
    EDC = flipud(cumsum(flipud(Impulsantwort.^2)));
    EDC2 = 20*log10(sqrt(EDC)/(2*10^(-5)));
    L_5dB = find(EDC2 < max(EDC2)-5,1);
    L_vardB = find(EDC2 < max(EDC2)-5-Ab,1);

    switch DTFvX                %DTF-Variantenbestimmung
        case 0
            DTF = 0;
        case 1
            % Berechnung des DTF mittels hallkorrigierter Impulsantwort (Hanyu 2014)
            DTF = DTF_Hanyu_alt(Impulsantwort,L_5dB,L_vardB,Sw);

        case 2
            % Alternative Berechnung des DTF (Hanyu 2018)
            [DTF,A_t] = DTF_Hanyu_neu_test(Impulsantwort,L_5dB,L_vardB,Sw);
            
    end

    switch TvX                  %Nachhallzeitvariantenbestimmung
        case 0
            T = 0;
        case 1
            % Nachhallzeitberechnung mittels Regression
            [T,Reg,Fehler] = RT_Regression(Impulsantwort,fs,L_5dB,L_vardB);

        case 2
            % Nachhallzeitberechnung aus der hallkorrigierten Impulsantwort (nach Hanyu 2014)
            T = RT_Hanyu_DCIR(Impulsantwort,L_5dB,L_vardB,fs);
            
        case 3
            % Alternative Nachhallzeitbestimmung (mittels Hanyu 2018)
            T = 13.82/(mean(A_t)*fs);
    end
    
        %Platzhalter für Tests zur Klimakorrektur der Nachhallzeit
            switch KlKorrvX
                case 1
                    c = 331+0.6*((Klimadaten(x1+1,1)+Klimadaten(x1+1,4))/2);
                    alpha = Luftdaempfung((Klimadaten(x1+1,1)+Klimadaten(x1+1,4))/2,(Klimadaten(x1+1,2)+Klimadaten(x1+1,5))/2,(Klimadaten(x1+1,3)+Klimadaten(x1+1,6))/2,1000);
                    m = alpha/(10*log10(exp(1)));
                    %T = ...
                case 0
                    T = T;
            end
        %TEST ende
    

    % Schallenergieberechnung aus der Impulsantwort
    p = sqrt(sum(Impulsantwort.^2));

    % Ergebnisstrukturierung
    DTF_ges(x4+((x3-1)*(C-Ca+1)*R)+((x2-Ca)*R),x1+1,i) = DTF;
    T_ges(x4+((x3-1)*(C-Ca+1)*R)+((x2-Ca)*R),x1+1,i) = T;
    p_ges(x4+((x3-1)*(C-Ca+1)*R)+((x2-Ca)*R),x1+1,i) = p;

            end                 %Ende der Mikrofonpositionsschleife
        end                     %Ende der Lautsprecherpositionsschleife
    end                         %Ende der Raumzustandsschleife
    
    i2 = 1;
                                %Ergebnisstrukturierung
        DTF_ges(((C-Ca+1)*S*R)+2,x1+1,i2) = mean(DTF_ges(1:(C-Ca+1)*S*R,x1+1,i2));
        T_ges(((C-Ca+1)*S*R)+2,x1+1,i2) = mean(T_ges(1:(C-Ca+1)*S*R,x1+1,i2));
        T_ges(((C-Ca+1)*S*R)+3,x1+1,i2) = std(T_ges(1:(C-Ca+1)*S*R,x1+1,i2))/T_ges(((C-Ca+1)*S*R)+2,x1+1,i2);
        p_ges(((C-Ca+1)*S*R)+2,x1+1,i2) = mean(p_ges(1:(C-Ca+1)*S*R,x1+1,i2));
        p_ges(((C-Ca+1)*S*R)+3,x1+1,i2) = std(p_ges(1:(C-Ca+1)*S*R,x1+1,i2))/p_ges(((C-Ca+1)*S*R)+2,x1+1,i2);
        
        switch AvX              %Testweise Berechnung von Aeq
            case 0
                A = 0;
            case 1 && (TvX ~= 0)
                c = 331+0.6*((Klimadaten(x1+1,1)+Klimadaten(x1+1,4))/2);
                alpha = Luftdaempfung((Klimadaten(x1+1,1)+Klimadaten(x1+1,4))/2,(Klimadaten(x1+1,2)+Klimadaten(x1+1,5))/2,(Klimadaten(x1+1,3)+Klimadaten(x1+1,6))/2,1000);
                m = alpha/(10*log10(exp(1)));
                A = 55.3*(VolModellhallraum(x1)/(c*T_ges(((C-Ca+1)*S*R)+2,x1+1,i2)))-4*VolModellhallraum(x1)*m;
                A_ges(i2,x1+1) = A;

        end
        
end                             %Ende der Zustandsschleife der Diffusoren
                    
                                %Plotten der vorläufigen Ergebnisse
figure(1)
plot(0:18,DTF_ges(14,:));xlabel('Z (Anzahl Diffusoren)'); ylabel('DTF');title('Degree of Time Fluctuation (DTF) mit zunehmender Diffusoroberfläche');axis([0,20,0,12.5]);grid on; grid minor
figure(2)
plot(0:18,T_ges(14,:));xlabel('Z (Anzahl Diffusoren)'); ylabel('T [s]');title('Nachhallzeit mit zunehmender Diffusoroberfläche');axis([0,20,0,0.7])
figure(3)
plot(0:18,T_ges(15,:));xlabel('Z (Anzahl Diffusoren)'); ylabel('\sigma');title('rel. Standardabweichung von T20 mit zunehmender Diffusoroberfläche');axis([0,20,0.0,0.7])
figure(4)
plot(0:18,p_ges(14,:));xlabel('Z (Anzahl Diffusoren)'); ylabel('E [Pa s²]');title('Schallenergie mit zunehmender Diffusoroberfläche');
figure(5)
plot(0:18,p_ges(15,:));xlabel('Z (Anzahl Diffusoren)'); ylabel('\sigma');title('rel. Standardabweichung der Schallenergie mit zunehmender Diffusoroberfläche');axis([0,20,0.0,3])
figure(6)
fschroeder = 2000*sqrt(T_ges(14,:)/VolModellhallraum(1:18));
plot(0:18,fschroeder);xlabel('Z (Anzahl Diffusoren)'); ylabel('f_S_c_h_r_o_e_d_e_r [Hz]');title('Schröderfrequez mit zunehmender Diffusoroberfläche');

%Weitere Plots wurden mit einem seperaten Skript erstellt.