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
%% Uebungszettel-Nr: Blatt 03
%% Aufgabennummer:   3.1
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
% Output: Die genäherten Werte für das Integral aus Aufgabe 3.1a)
%

clear all;

% -----------------------------Funktionen--------------------------------------
function [c]  = myQuadratur1D(f,w,x)
  c = sum(w.*f(x));
endfunction

% --------------------------Function_handles-----------------------------------
f = @(x) e.^x.+1;

% --------------------------------Main-----------------------------------------
I1 = myQuadratur1D(f,[0.5 0.5],[0 1])
I2 = myQuadratur1D(f,[1/6 4/6 1/6],[0 1/2 1])

% ---------------------------------Ausgabe-------------------------------------
%I1 =  2.8591
%I2 =  2.7189

% ---------------------------------Anmerkungen---------------------------------
% Der Test verlief erfolgreich