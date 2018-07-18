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
%% Uebungszettel-Nr: Blatt 07
%% Aufgabennummer:   7.1
%% Program name:     myNewtonTest.m 
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
g = @(x) x; 
f = @(x) cos(2.*x).^2.-x.^2;
df = @(x) -4.*sin(2.*x).*cos(2.*x).-2.*x;

% ------------------------------ constants ------------------------------------

% ------------------------------ functions ------------------------------------
function [x,e,v] = myNewton(f,df,x0)
  x = x0-f(x0)/df(x0);
  v = f(x);
  e = abs(x-x0);
  if e < 10^(-12) || length(x)>=50 
    return;
  end
  [x1,e1,v1] = myNewton(f,df,x);
  x = [x,x1];
  e = [e,e1];
  v = [v,v1];
endfunction

function [x,e,v] = mybisect(f,x00,x0)
  x = (x00+x0)/2;
  v = f(x);
  e = x0-x00;
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
[xN,eN,vN] = myNewton(f,df,0.75);
xEndNewton = xN(end)
eEndNewton = eN(end)
vEndNewton = vN(end)
iterationsNewton = length(xN)
[xb,eb,vb] = mybisect(f,0,0.75);
xEndbisect = xb(end)
eEndbisect = eb(end)
vEndbisect = vb(end)
iterationsbisect = length(xb)

cla();
clf();
k = 1:length(eN);
semilogy(k, eN(k));
hold on;
semilogy(k, eb(k));
xlabel("k");
ylabel("Fehler(k)");
xlim([1 k(end)]);
legend({"Newton", "Bisect"},"location", "southwest");
title("Fehler e der Nullstelle in Abhängigkeit des Iterationsschrittes k");
print ("PA7.1.fig");
print ("PA7.1.jpg");

% ---------------------------------Ausgabe-------------------------------------
%{
xEndNewton =  0.514933264661129
eEndNewton =   1.11022302462516e-016
vEndNewton = 0
iterationsNewton =  5
xEndbisect =  0.514933264660954
eEndbisect =   6.82121026329696e-013
vEndbisect =   4.91329199547863e-013
iterationsbisect =  41
%}

% ---------------------------------Anmerkungen---------------------------------
%{ 
Das Newtonverfahren konvergiert deutlich schneller (quadratisch).
Die numerische Konvergenzrate stimmt mit der aus der Theorie zu erwartenden 
Konvergenzgeschwindigkeit überein.
%}