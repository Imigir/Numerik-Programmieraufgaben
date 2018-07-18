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
%% Uebungszettel-Nr: Blatt 02
%% Aufgabennummer:   1.1, 1.2
%% Program name:     Blatt02.m 
%%
%% Program(version): Octave-4.22
%% OS:               Windows 10 64bit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description of the program
%                  
% Intput: none
%
% Output: Graphen der nähernden Polynome, sowie deren Abweichungen zur Originalfunktion
%

clear all;
format long;

% -----------------------------Funktionen--------------------------------------
function y = DividierteDifferenzen(xi, fi) % Rekursive Implementation der Dividierten Differenzen
  if length(xi)==1 % Abbruchbedingung
    y = fi(1);
  else
    % Rekursionsformel aus Vorlesung
    y = (DividierteDifferenzen(xi(2:length(xi)), fi(2:length(fi)))-DividierteDifferenzen(xi(1:length(xi)-1), fi(1:length(fi)-1)))/(xi(length(xi))-xi(1));
  end
endfunction

function c = myNewtonInterpol(xi, fi) % Berechnung der Newton-Koeffizienten
  for i=1:length(xi)
    c(i) = DividierteDifferenzen(xi(1:i), fi(1:i));
  end
endfunction

function N = NewtonPolynom(x, i, xi) % Rekursive Implementierung des Newton-Polynoms
  if i==1 % Abbruchbedingung
    N = 1;
  else
    % Rekursionsformel aus Vorlesung
    N = (x.-xi(i-1)).*NewtonPolynom(x, i-1, xi);
  end
endfunction

function p = Polynom(x, xi, ci) % Berechnung des Interpolationspolynoms
  p = 0;
  for i=1:length(xi)
    p .+= ci(i).*NewtonPolynom(x, i, xi);
  end
endfunction

function x1 = aequidistanteKnoten(n) % Berechnung der aequisistanten Stützstellen
  for i=0:n
    x1(i+1) = -1+2*i/n;
  end
endfunction

function x2 = TschebyscheffKnoten(n) % Berechnung der Tschebyscheff´schen Stützstellen
  for i=0:n
    x2(i+1) = cos((2*i+1)*pi/(2*n+2));
  end
endfunction

% --------------------------Function_handles-----------------------------------
f = @(x) 1./(1.+25.*x.^2); % Die Runge-Funktion

% -----------------------Ausführen des Programms-------------------------------
for i=7:5:17
  h=figure;
  x = -1:5*10^(-5):1;
  % Berechnen der Stützstellen
  x1 = aequidistanteKnoten(i);
  x2 = TschebyscheffKnoten(i);
  % Berechnen der Koeffizienten
  c1 = myNewtonInterpol(x1, f(x1));
  c2 = myNewtonInterpol(x2, f(x2));
  % Berechnen der Interpolationspolynome
  p1 = @(x) Polynom(x, x1, c1);
  p2 = @(x) Polynom(x, x2, c2);
  % Plooten der Graphen
  cla();
  clf();
  hold on;
  plot(x, f(x), 'g');
  scatter(x1,f(x1), 'r', '*');
  plot(x, p1(x), 'r');
  scatter(x2,f(x2), 'b', '*');
  plot(x, p2(x), 'b');
  hold off;
  xlabel("x");
  ylabel("f(x)");
  xlim([-1 1]);
  title("Newton Naeherung bei verschiedenen Stuetzpunkten");
  saveas(h,sprintf('PA2-1-N%d.png',i)); 
  % Plotten mit y beschränkt auf [0 1] zur besseren Anschauung (unnötig, kann aber bei Bedarf einkommentiert werden)
  %cla();
  %clf();
  %hold on;
  %plot(x, f(x), 'g');
  %scatter(x1,f(x1), 'r', '*');
  %plot(x, p1(x), 'r');
  %scatter(x2,f(x2), 'b', '*');
  %plot(x, p2(x), 'b');
  %hold off;
  %xlabel("x");
  %ylabel("f(x)");
  %xlim([-1 1]);
  %ylim([0 1]);
  %title("Newton Näherung");
  %saveas(h,sprintf('PA2-1-N%d-2.png',i)); 
  
  x = -1:2*10^(-2):1;
  Fp1 = max(abs(f(x)-p1(x)));
  Fp2 = max(abs(f(x)-p2(x)));
  printf("n=%2.0f:\n",i);
  printf("Fehler p1(rot)=%2.10f\n",Fp1);
  printf("Fehler p2(blau)=%2.10f\n",Fp2);
end

% ---------------------------------Ausgabe-------------------------------------
%n= 7:
%Fehler p1(rot)=0.2473586066
%Fehler p2(blau)=0.3917402846
%n=12:
%Fehler p1(rot)=3.6052744506
%Fehler p2(blau)=0.0690810331
%n=17:
%Fehler p1(rot)=4.0635507433
%Fehler p2(blau)=0.0559074779

% -------------------------------Beobachtungen---------------------------------
% Nah um 0 liefert der rote Graph (p1) die besseren Ergebnisse, doch je weiter 
% sich x von 0 entfernt, desto schlechter wird die Näherung. 
% Bei dem blauen Graphen (p2) ist dies nicht der Fall, da hier umso mehr
% Stützstellen genommen werden, desto weiter sich x von 0 entfernt.
% Bei größerem n wird die Näherung des blauen Graphen besser, und auch der rote 
% Graph wird in der Nähe von 0 besser. Weiter außen wird der rote Graph dafür
% drastisch schlechter, was den maximalen Fehler im Endeffekt erhöht.

% ---------------------------------Anmerkungen---------------------------------
% Das Kompilieren kann etwas länger dauern