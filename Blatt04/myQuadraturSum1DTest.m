%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Numerische Mathematik fuer Physik und Ingenieurwissenschaften 2018     %%%
%%   Programmierabgaben (Praktischer Teil des Uebungungsplattes)            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Student 1: Jonah, Blank
%%  Unimail-adress: jonah.blank@tu-dortmund.de
%%
%%  Student 2: Lars, Kolk
%%  Unimail-adress: lars.kolk@tu-dortmund.de
%%
%%  Student 3: David, Rolf
%%  Unimail-adress: david.rolf@tu-dortmund.de
%%
%% Uebungszettel-Nr: Blatt 04
%% Aufgabennummer:   4.1
%% Program name:     myQuadratur1DTest.m 
%%
%% Program(version): Octave-4.22
%% OS:               Windows 10 64bit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description of the program
%                  
% Intput: none
%
% Output: Die genäherten Werte für die Integrale der Definierten Funktionen für
% N=2 und die Graphen für die Fehler in Abhängigkeit von N
%
clear all;
close all;
format long;

% ---------------------------- function handles -------------------------------
f = @(x) e.^x.+1; 
g = @(x) 1./(1.+25.*x.^2); % Runge-Funktion

% ------------------------------ constants ------------------------------------
Intg = 0.5493603067780063443445088; % Wert des Integrals über g von -1 bis 1
Intf = 2.7182818284590452353602875; % Wert des Integrals über f von 0 bis 1
% ------------------------------ functions ------------------------------------
function [c]  = myQuadratur1Dsub(f,w,x)
  c = sum(w.*f(x));
endfunction

function [v]  = myQuadratur1D(f,w,p,a,b,N)
  v = 0;
  h = (b-a)/N; % Länge der Teilintervalle
  for i=1:N
    Ia = (a.+(i.-1).*h); % untere Teilintervallgrenze
    Ib = (a.+i.*h); % obere Teilintervallgrenze
    x = (Ib-Ia)*p.+Ia; % Berechnung der Stützstellen im Teilintervall nach übergebenem Schema
    z = (Ib-Ia)*w; % Berechnung der Gewiichte im Teilinterwall
    v += myQuadratur1Dsub(f,z,x); % aufsummieren der Teilintevalle
  end
endfunction

% ------------------------------- main ----------------------------------------

% ---Plotte Fehler für f---
for i = 1:1000
  FehlerTrapezf(i) = abs((myQuadratur1D(f,[0.5 0.5],[0 1],0,1,i)-Intf)/Intf*100);
  FehlerSimpsonf(i) = abs((myQuadratur1D(f,[1/6 4/6 1/6],[0 1/2 1],0,1,i)-Intf)/Intf*100);
  if FehlerSimpsonf(i)==0 
    FehlerSimpsonf(i) = 1*10^(-13); % ln(0) gäbe sonst einen Fehler
  end
end
N = 1:1000;
cla();
clf();
semilogy(N, FehlerTrapezf, ";Trapez;");
hold on;
semilogy(N, FehlerSimpsonf, ";Simpson;");
xlabel("N");
ylabel("Fehler(N)");
xlim([0 1000]);
legend("show");
title("Fehler zu f mit der Trapez- und Simpson-Regel");
print ("myQuadratsum1DPlotf.fig");
print ("myQuadratsum1DPlotf.jpg");

% ---Plotte Fehler für g---
for i = 1:1000
  FehlerTrapezg(i) = abs((myQuadratur1D(g,[0.5 0.5],[0 1],-1,1,i)-Intg)/Intg*100);
  FehlerSimpsong(i) = abs((myQuadratur1D(g,[1/6 4/6 1/6],[0 1/2 1],-1,1,i)-Intg)/Intg*100);
end
N = 1:1000;
figure;
cla();
clf();
semilogy(N, FehlerTrapezg, ";Trapez;");
hold on;
semilogy(N, FehlerSimpsong, ";Simpson;");
xlabel("N");
ylabel("Fehler(N)");
xlim([0 1000]);
legend("show");
title("Fehler zu g mit der Trapez- und Simpson-Regel");
print ("myQuadratsum1DPlotg.fig");
print ("myQuadratsum1DPlotg.jpg");

% berechnung über Aufruf der Funktion nach Aufgabenstellung für N=2
printf("Naeherung der Integrale fuer N=2: \n")
IntfTrapez = myQuadratur1D(f,[0.5 0.5],[0 1],0,1,2)
IntfSimpson = myQuadratur1D(f,[1/6 4/6 1/6],[0 1/2 1],0,1,2)
IntgTrapez = myQuadratur1D(g,[0.5 0.5],[0 1],-1,1,2)
IntgSimpson = myQuadratur1D(g,[1/6 4/6 1/6],[0 1/2 1],-1,1,2)

% ---------------------------------Ausgabe-------------------------------------
%{
Naeherung der Integrale fuer N=2:
IntfTrapez =  2.75393109246483
IntfSimpson =  2.71831884192175
IntgTrapez =  1.03846153846154
IntgSimpson =  0.530061892130858
%}

% ---------------------------------Anmerkungen---------------------------------
%{ 
Die Simpson-Regel ist besser.
Wenn die Abweichung zu gering wird entsteht ein "Untergrundrauschen".
Bei der Rungefunktion hat die Simpson-Regel einen Ausreisser bei ca. N=20.
%}