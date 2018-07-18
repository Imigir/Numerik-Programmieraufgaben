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
%% Uebungszettel-Nr: Blatt 06
%% Aufgabennummer:   6.1
%% Program name:     myBisectTest.m 
%%
%% Program(version): Octave-4.22
%% OS:               Windows 10 64bit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description of the program
%                  
% Intput: none
%
% Output: siehe Ausgabe
%
clear all;
close all;

% ---------------------------- function handles -------------------------------
g = @(x) x; 
f = @(x) cos(2.*x).^2-x.^2;

% ------------------------------ constants ------------------------------------

% ------------------------------ functions ------------------------------------
function [x,e,v] = mybisect(f,x00,x0)
  x = (x00+x0)/2;
  v = f(x);
  if f(x) == 0
    e = 0;
    return;
  else
    e = x0-x00;
  end
  if e < 10^(-12)
    return;
  end
  if (f(x00)<0 && f(x)<0) || (f(x00)>0 && f(x)>0) 
    [x1,e1,v1] = mybisect(f,x,x0);
  else
    [x1,e1,v1] = mybisect(f,x00,x);
  end
  x = [x,x1];
  e = [e,e1];
  v = [v,v1];
endfunction

% ------------------------------- main ----------------------------------------
[x,e,v] = mybisect(f,0,0.75);
xEnd = x(end)
eEnd = e(end)
vEnd = v(end)

cla();
clf();
k = 1:length(e);
semilogy(k, e(k));
hold on;
xlabel("k");
ylabel("Fehler(k)");
xlim([1 k(end)]);
title("Fehler e der Nullstelle in Abhängigkeit des Iterationsschrittes k");
print ("PA7.1.fig");
print ("PA7.1.jpg");

% ---------------------------------Ausgabe-------------------------------------
%{
xEnd =  0.51493
eEnd =   6.8212e-013
vEnd =   4.9133e-013
%}

% ---------------------------------Anmerkungen---------------------------------
%{ 
Exponentielle Abnahme des Fehlers (e^(-k)).
Die numerische Konvergenzrate stimmt mit der aus der Theorie zu erwartenden 
Konvergenzgeschwindigkeit überein.
%}