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
%% Uebungszettel-Nr: Blatt 1
%% Aufgabennummer:   1.1, 1.2
%% Program name:     Blatt1 
%%
%% Program(version): Octave-4.22
%% OS:               Windows 10 64bit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description of the program
%                  
% Intput: none
%
% Output: Werte der Naeherungen und deren Abweichungen
%         Plots des Polynoms P_7(x) mit unterschiedlicher Darstellung
%

clear all;
format long;

% ------Die Funktion zu Naeherung der Exponentialfunktion über die Summe-------
function y = expsum(x, n)
  k = 0:n;
  y = sum(x.^k./factorial(k));
endfunction

% -----------------------------Konstanten--------------------------------------
el = 0.00408677143846406; % Der Wert für e^(-5.5) von Wolfram Alpha

% --------------------------Function_handles-----------------------------------
P = @(x) (x.-1).^7; % Das Polynom P_7
PH = @(x) ((((((x.-7).*x+21).*x.-35).*x.+35).*x.-21).*x.+7).*x.-1; % Das Polynom P_7 in Horner-Form

% -------Ausführen der expsum Funktion für die verschiedenen Naeherungen-------
% ---a)---
printf("Naeherung a): \n"); 
for i=3:3:30
  ew = expsum(-5.5,i);
  printf("n=%2.0f: e^(-5.5) = %21.17f, Abweichung: %25.17f%% \n",i,ew, (ew-el)/el*100);
end
% ---b)---
printf("Naeherung b): \n"); 
for i=3:3:30
  ew = 1/expsum(5.5,i);
  printf("n=%2.0f: e^(-5.5) = %21.17f, Abweichung: %25.17f%% \n",i,ew, (ew-el)/el*100);
end
% ---c)---
printf("Naeherung c): \n"); 
for i=3:3:30
  ew = expsum(-0.5,i)^11;
  printf("n=%2.0f: e^(-5.5) = %21.17f, Abweichung: %25.17f%% \n",i,ew, (ew-el)/el*100);
end

% -----------Plotten von P(x) und PH(x) im Bereich von 0.8 bis 1.2-------------
x = 0.8:5*10^(-5):1.2;
x = single(x);
% ---Plotte P(x)---
cla();
clf();
plot(x, P(x));
xlabel("x");
ylabel("P(x)");
xlim([0.8 1.2]);
title("P(x) genaehert");
print ("Graph1.jpg");
hold on;
figure;
% ---Plotte PH(x)---
cla();
clf();
plot(x, PH(x));
xlabel("x");
ylabel("P(x)");
xlim([0.8 1.2]);
title("P(x) genaehert in Horner-Form");
print ("Graph2.jpg"); 

% ---------------------------------Ausgabe-------------------------------------
%Naeherung a):
%n= 3: e^(-5.5) = -17.10416666666666785, Abweichung: -418625.15914359438465908%
%n= 6: e^(-5.5) =  17.52840711805555074, Abweichung:  428805.98072309355484322%
%n= 9: e^(-5.5) =  -4.60245118586895430, Abweichung: -112718.26738220290280879%
%n=12: e^(-5.5) =   0.48723105192884963, Abweichung:   11822.15075555991461442%
%n=15: e^(-5.5) =  -0.02114457410215519, Abweichung:    -617.39066939603549145%
%n=18: e^(-5.5) =   0.00483716591494927, Abweichung:      18.36154744115653870%
%n=21: e^(-5.5) =   0.00407285902525824, Abweichung:      -0.34042552697894757%
%n=24: e^(-5.5) =   0.00408694306429718, Abweichung:       0.00419954567329898%
%n=27: e^(-5.5) =   0.00408676995860855, Abweichung:      -0.00003621087031123%
%n=30: e^(-5.5) =   0.00408677144773050, Abweichung:       0.00000022674233323%
%Naeherung b):
%n= 3: e^(-5.5) =   0.02026171380329253, Abweichung:     395.78779015122836427%
%n= 6: e^(-5.5) =   0.00595708032220462, Abweichung:      45.76494946934145247%
%n= 9: e^(-5.5) =   0.00431903839403954, Abweichung:       5.68338501608929736%
%n=12: e^(-5.5) =   0.00410504250047283, Abweichung:       0.44707814674442431%
%n=15: e^(-5.5) =   0.00408758964263671, Abweichung:       0.02002079599927505%
%n=18: e^(-5.5) =   0.00408679339021944, Abweichung:       0.00053714174408482%
%n=21: e^(-5.5) =   0.00408677181601334, Abweichung:       0.00000923832622198%
%n=24: e^(-5.5) =   0.00408677144286244, Abweichung:       0.00000010762487616%
%n=27: e^(-5.5) =   0.00408677143850032, Abweichung:       0.00000000088727558%
%n=30: e^(-5.5) =   0.00408677143846429, Abweichung:       0.00000000000556059%
%Naeherung c):
%n= 3: e^(-5.5) =   0.00391493328663423, Abweichung:      -4.20474094079527205%
%n= 6: e^(-5.5) =   0.00408687952837399, Abweichung:       0.00264487289189706%
%n= 9: e^(-5.5) =   0.00408677141938826, Abweichung:      -0.00000046676939290%
%n=12: e^(-5.5) =   0.00408677143846546, Abweichung:       0.00000000003438230%
%n=15: e^(-5.5) =   0.00408677143846406, Abweichung:      -0.00000000000002122%
%n=18: e^(-5.5) =   0.00408677143846406, Abweichung:      -0.00000000000002122%
%n=21: e^(-5.5) =   0.00408677143846406, Abweichung:      -0.00000000000002122%
%n=24: e^(-5.5) =   0.00408677143846406, Abweichung:      -0.00000000000002122%
%n=27: e^(-5.5) =   0.00408677143846406, Abweichung:      -0.00000000000002122%
%n=30: e^(-5.5) =   0.00408677143846406, Abweichung:      -0.00000000000002122%

% -------------------------------Erklärungen-----------------------------------
% In a) wechselt das Vorzeichen, was wegen der grossen Werte zu grossen Abweichungen führt.
% In b) wechselt das Vorzeichen nicht, weshalb die Naeherung genauer ist.
% In c) wechselt das Vorzeichen wieder, doch x ist kleiner als 1, weshalb hier gute Naeherungen erziehlt werden.
% Graph2 ist schlechter genaehert, da hier durch die vielen Rechenoperationen groessere Abweichungen entstehen.
