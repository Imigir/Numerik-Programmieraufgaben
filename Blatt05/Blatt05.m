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
%% Uebungszettel-Nr: Blatt 05
%% Aufgabennummer:   5.1
%% Program name:     Blatt05.m 
%%
%% Program(version): Octave-4.22
%% OS:               Windows 10 64bit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description of the program
%                  
% Intput: none
%
% Output: Dimension n, Zeit t, Residuen Res
%
clear all;
close all;

% ---------------------------- function handles -------------------------------
f = @(x) e.^x.+1; 
g = @(x) 1./(1.+25.*x.^2); % Runge-Funktion

% ------------------------------ constants ------------------------------------
Intg = 0.5493603067780063443445088; % Wert des Integrals über g von -1 bis 1
Intf = 2.7182818284590452353602875; % Wert des Integrals über f von 0 bis 1
% ------------------------------ functions ------------------------------------
function LR  = thomas_decompose(A)
  LR = A;
  for i=2:length(A)
    a = LR(i,i-1)/LR(i-1,i-1);
    LR(i,i) = LR(i-1,i)*(-a)+LR(i,i);
    LR(i,i-1) = a;
  end
endfunction

function y  = forward_solve(LR,b)
  y = b;
  for i=2:length(b)
    y(i) = b(i)-LR(i,i-1)*y(i-1);
  end
endfunction

function x  = backward_solve(LR,y)
  x = y;
  n = length(y)
  x(n) = y(n)/LR(n,n); 
  for i=n-1:-1:1
    x(i) = (y(i)-LR(i,i+1)*x(i+1))/LR(i,i);
  end
endfunction

function t = LR_test(n)
  A = gallery('tridiag', n, -2, 8, -2);
  tic;
  LR = thomas_decompose(A);
  t = toc;
  b = 1:n;
  for i=1:floor(n/2)
    b(i) = 1;
  end
  for i=floor(n/2)+1:n
    b(i) = 2;
  end
  y = forward_solve(LR,b);
  x = backward_solve(LR,y);
  Res = sqrt(sum((A*x'-b').^2))
endfunction


% ------------------------------- main ----------------------------------------
for n=[10000,100000,1000000,10000000]
  t = LR_test(n)
end


% ---------------------------------Ausgabe-------------------------------------
%{
n =  10000
Res =   3.1310e-014
t =  0.33280
>> Blatt05
n =  10000
Res =   3.1310e-014
t =  0.33280
n =  100000
Res =   9.9272e-014
t =  3.3161
>> Blatt05
n =  10000
Res =   3.1310e-014
t =  0.33481
n =  100000
Res =   9.9272e-014
t =  3.3441
n =  1000000
Res =   3.1401e-013
t =  33.277
n =  10000000
Res =   9.9301e-013
t =  336.50
%}

% ---------------------------------Anmerkungen---------------------------------
%{ 
Es ist zu erkennen, dass die Laufzeit der erwarteten Laufzeit von O(n) entspricht.
=> alpha ist 1.
Die LR-Zerlegung.zip ist deutlih langsamer und braucht bereits bei n=1000 
eine längere Zeit, als die in der Aufgabe implementierte Funktion für n=100000 
Die Laufzeit scheint hier exponentiell (zumindest sehr stark) zu wachsen.
%}