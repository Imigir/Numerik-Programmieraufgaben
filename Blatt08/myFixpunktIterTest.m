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
%% Uebungszettel-Nr: Blatt 08
%% Aufgabennummer:   8.1
%% Program name:     myFixpunktIter.m 
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
format long;

% ---------------------------- function handles -------------------------------
f = @(x) 1./6.*x.^3.+1./3.*x.^2.+1./6; %Die Funktion des Fixpunktproblems

% ------------------------------ constants ------------------------------------
x_lit = 0.17819408002075; %Der Wert der Nullstelle von Wolfram Alpha

% ------------------------------ functions ------------------------------------
function [x,e] = myFixpunktIter(f,x0,N) %Stelle x, Schrittweite e (Maß für den Fehler)
  x = f(x0);
  e = abs(x-x0);
  for i=2:N
    x = [x,f(x(i-1))];
    e = [e,abs(x(i)-x(i-1))];
  end
endfunction

% ------------------------------- main ----------------------------------------
[x,e] = myFixpunktIter(f,0,8); 
v = abs(x.-x_lit); %wirklicher Fehler

%output
for i=1:length(e)
  printf("x%1.0f = %1.15f\n" , i, x(i));
  printf("e%1.0f = %1.15f\n" , i, e(i));
  printf("v%1.0f = %1.15f\n" , i, v(i));
  printf("\n");
end

%plot von e und v in Abhängigkeit von k
cla();
clf();
k = 1:length(e);
semilogy(k, e(k));
hold on;
semilogy(k, v(k));
xlabel("k");
ylabel("e(k) und v(k)");
xlim([1 k(end)]);
legend({"e(k)", "v(k)"});
title("e und v in Abhaengigkeit des Iterationsschrittes k");
print ("PA7.1.jpg");

% ---------------------------------Ausgabe-------------------------------------
%{
x1 = 0.166666666666667
e1 = 0.166666666666667
v1 = 0.011527413354083

x2 = 0.176697530864198
e2 = 0.010030864197531
v2 = 0.001496549156552

x3 = 0.177993481368762
e3 = 0.001295950504565
v3 = 0.000200598651988

x4 = 0.178167081872216
e4 = 0.000173600503454
v4 = 0.000026998148534

x5 = 0.178190444417306
e5 = 0.000023362545090
v5 = 0.000003635603444

x6 = 0.178193590410104
e6 = 0.000003145992798
v6 = 0.000000489610646

x7 = 0.178194014083709
e7 = 0.000000423673605
v7 = 0.000000065937041

x8 = 0.178194071140843
e8 = 0.000000057057134
v8 = 0.000000008879907
%}

% ---------------------------------Anmerkungen---------------------------------
%{ 
Es werden tatsächlich nur drei Iterationen für einen Fehler kleiner als 10^(-3)
benötigt (anstatt 8).
Der Fehler konvergiert linear.
%}