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
%% Uebungszettel-Nr: Blatt 09
%% Aufgabennummer:   9.1
%% Program name:     PA10.m 
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

% ------------------------------ constants ------------------------------------
eps = 1*10^(-6);  % error tolerance
maxit = 1*10^8;   % maximal number of iterations
start = 10;       % start number for n
ende = 20;        % end number for n

% ----------------------------- variables -------------------------------------
vectorOnly = true;      % changing this to false is not recommended

% ---------------------------- function handles -------------------------------
f = @(x) 2./(1.+sin(pi./x)); %function to calculate omega(n)

% ------------------------------ functions ------------------------------------
function [x,numit] = my_jacobi(A,b,x0,eps,maxit) 
  % initialisierung
  numit = 0;
  xalt = x0;
  x = x0;
  relEinDef = eps*norm(A*x0-b)/norm(x0); % relativer Eingangsdefekt mal eps
  while norm(A*x-b) > relEinDef
    % Fehlermeldung
    if numit > maxit
      errordlg('maximale Anzahl Iterationen ueberschritten!')
      return;
    end
    % Zeilenweise implementation des Jacobi Algorithmus
    x = b;
    for i=1:length(b)
      for j=1:length(b)
        if j != i
          x(i) = x(i)-A(i,j)*xalt(j);
        end
      end
      x(i) = x(i)/A(i,i);
    end
    xalt = x;
    numit++;
  end
endfunction

function [x,numit] = my_jacobi_vector(A,b,x0,eps,maxit) 
  % initialisierung
  numit = 0;
  D = diag(diag(A));  % Diagonalmatrix
  invD = inv(D);
  LU = A-D;           % Matrix ohne Diagonale
  J = -invD*LU;       % Jacobi-Iterationsmatrix
  xalt = x0;
  x = x0;
  relEinDef = eps*norm(A*x0-b)/norm(x0); % relativer Eingangsdefekt mal eps
  while norm(A*x-b) > relEinDef
    % Fehlermeldung
    if numit > maxit
      errordlg('maximale Anzahl Iterationen ueberschritten!')
      return;
    end
    % Vektorielle implementation des Jacobi Algorithmus
    x = invD*b+J*xalt;
    xalt = x;
    numit++;
  end
endfunction

function [x,numit] = my_gauss_seidel(A,b,x0,eps,maxit) 
  % initialisierung
  numit = 0;
  xalt = x0;
  x = x0;
  relEinDef = eps*norm(A*x0-b)/norm(x0); % relativer Eingangsdefekt mal eps
  while norm(A*x-b) > relEinDef
    % Fehlermeldung
    if numit > maxit
      errordlg('maximale Anzahl Iterationen ueberschritten!')
      return;
    end
    % Zeilenweise implementation des Gauss Seidel Algorithmus
    x = b;
    for i=1:length(b)
      for j=1:i-1
        x(i) = x(i)-A(i,j)*x(j);
      end
      for j=i+1:length(b)
        x(i) = x(i)-A(i,j)*xalt(j);
      end
      x(i) = x(i)/A(i,i);
    end
    xalt = x;
    numit++;
  end
endfunction

function [x,numit] = my_gauss_seidel_vector(A,b,x0,eps,maxit) 
  % initialisierung
  numit = 0;
  DL = tril(A);       % untere Dreiecksmatrix mit Diagonale
  invDL = inv(DL);
  U = A-DL;      % obere Dreiecksmatrix ohne Diagonale
  H = -invDL*U;       % Gausß-Seidel-Iterationsmatrix
  xalt = x0;
  x = x0;
  relEinDef = eps*norm(A*x0-b)/norm(x0); % relativer Eingangsdefekt mal eps
  while norm(A*x-b) > relEinDef
    % Fehlermeldung
    if numit > maxit
      errordlg('maximale Anzahl Iterationen ueberschritten!')
      return;
    end
    % Vektorielle implementation des Gauss Seidel Algorithmus
    x = invDL*b+H*xalt;
    xalt = x;
    numit++;
  end
endfunction

function [x,numit] = my_sor(A,b,x0,eps,maxit,omega) 
  % initialisierung
  numit = 0;
  xalt = x0;
  x = x0;
  relEinDef = eps*norm(A*x0-b)/norm(x0); % relativer Eingangsdefekt mal eps
  while norm(A*x-b) > relEinDef
    % Fehlermeldung
    if numit > maxit
      errordlg('maximale Anzahl Iterationen ueberschritten!')
      return;
    end
    % Zeilenweise implementation des Sor Algorithmus
    x = b;
    for i=1:length(b)
      for j=1:i-1
        x(i) = x(i)-A(i,j)*x(j);
      end
      for j=i+1:length(b)
        x(i) = x(i)-A(i,j)*xalt(j);
      end
      x(i) = x(i)*omega/A(i,i)+(1-omega)*xalt(i);
    end
    xalt = x;
    numit++;
  end
endfunction

function [x,numit] = my_sor_vector(A,b,x0,eps,maxit,omega) 
  % initialisierung
  numit = 0;
  D = diag(diag(A));  % Diagonalmatrix
  L = tril(A)-D;      % untere Dreiecksmatrix ohne Diagonale
  U = triu(A)-D;      % obere Dreiecksmatrix ohne Diagonale
  invDwL = inv(D+omega*L);
  H = invDwL*((1-omega)*D-omega*U);
  xalt = x0;
  x = x0;
  relEinDef = eps*norm(A*x0-b)/norm(x0); % relativer Eingangsdefekt mal eps
  while norm(A*x-b) > relEinDef
    % Fehlermeldung
    if numit > maxit
      errordlg('maximale Anzahl Iterationen ueberschritten!')
      return;
    end
    % Vektorielle implementation des Sor Algorithmus
    x = H*xalt+omega*invDwL*b;
    xalt = x;
    numit++;
  end
endfunction

function [A,b] = my_test_system(n)
  b = ones(n*n,1);
  T = gallery('tridiag', n, -1, 4, -1);
  Id = eye(n);    % unit matrix
  for i=0:n-1
    % add diagonal n*n submatrices T
    A(i*n+1:n*(i+1),i*n+1:n*(i+1)) = T;       
  end
  for i=0:n-2
    % add lower and upper n*n diagonal submatrices -Id
    A(i*n+1:n*(i+1),(i+1)*n+1:n*(i+2)) = -Id; 
    A((i+1)*n+1:n*(i+2),i*n+1:n*(i+1)) = -Id;
  end
endfunction

% ------------------------------- main ----------------------------------------

for n=start:ende
  % initialisierung
  [A,b] = my_test_system(n);
  x0 = b;
  omega = f(n);
  
  % Aufruf der Zeilenweisen Implementierung falls 'fast' == 'false' (way to slow)
  if !vectorOnly 
    tic;
    [x_jacobi,numit_jacobi(n-start+1)] = my_jacobi(A,b,x0,eps,maxit);
    t_jacobi(n-start+1) = toc;
    
    tic;
    [x_gauss_seidel,numit_gauss_seidel(n-start+1)] = my_gauss_seidel(A,b,x0,eps,maxit);
    t_gauss_seidel(n-start+1) = toc;
    
    tic;
    [x_sor,numit_sor(n-start+1)] = my_sor(A,b,x0,eps,maxit,omega);
    t_sor(n-start+1) = toc;
  end
  
  % Aufruf der Vektoriellen Implementierung
  tic;
  [x_jacobi_vector,numit_jacobi_vector(n-start+1)] = my_jacobi_vector(A,b,x0,eps,maxit);
  t_jacobi_vector(n-start+1) = toc;
  
  tic;
  [x_gauss_seidel_vector,numit_gauss_seidel_vector(n-start+1)] = my_gauss_seidel_vector(A,b,x0,eps,maxit);
  t_gauss_seidel_vector(n-start+1) = toc;
  
  tic;
  [x_sor_vector,numit_sor_vector(n-start+1)] = my_sor_vector(A,b,x0,eps,maxit,omega);
  t_sor_vector(n-start+1) = toc;
end

%output
printf("n = %1.0f:\n" , ende);
if !vectorOnly
  printf("Benoetigte Zeit Jacobi: %f\n" , t_jacobi(end));
  printf("Benoetigte Zeit Gauss Seidel: %f\n" , t_gauss_seidel(end));
  printf("Benoetigte Zeit Sor: %f\n" , t_sor(end));
end
printf("Benoetigte Zeit Jacobi Vektor: %f\n" , t_jacobi_vector(end));
printf("Benöoetigte Zeit Gauss Seidel Vektor: %f\n" , t_gauss_seidel_vector(end));
printf("Benoetigte Zeit Sor Vektor: %f\n" , t_sor_vector(end));

%plot der Anzahl der Iterationen in Abhängigkeit von n
cla();
clf();
figure('visible','off');
n = start:ende;
hold on;
if !vectorOnly
  plot(n, numit_jacobi, ";Jacobi;");
  plot(n, numit_gauss_seidel, ";Gauss Seidel;");
  plot(n, numit_sor, ";Sor;");
end
plot(n, numit_jacobi_vector, ";Jacobi Vektor;");
plot(n, numit_gauss_seidel_vector, ";Gauss Seidel Vektor;");
plot(n, numit_sor_vector, ";Sor Vektor;");
%hold off;
xlabel("n");
ylabel("Anzahl Iterationen");
xlim([n(1) n(end)]);
legend("show");
title("Anzahl der Iterationen in Abhaengigkeit von n");
print ("Fig1.jpg");

%plot der Zeit t in Abhängigkeit von n
%figure;
cla();
clf();
n = start:ende;
hold on;
if !vectorOnly
  plot(n, t_jacobi, ";Jacobi;");
  plot(n, t_gauss_seidel, ";Gauss Seidel;");
  plot(n, t_sor, ";Sor;");
end
plot(n, t_jacobi_vector, ";Jacobi Vektor;");
plot(n, t_gauss_seidel_vector, ";Gauss Seidel Vektor;");
plot(n, t_sor_vector, ";Sor Vektor;");
%hold off;
xlabel("n");
ylabel("t");
xlim([n(1) n(end)]);
legend("show");
title("t in Abhaengigkeit von n");
print ("Fig2.jpg");

% ---------------------------------Ausgabe-------------------------------------
%{
n = 60:
Benoetigte Zeit Jacobi Vektor: 2.557560
Benoetigte Zeit Gauss Seidel Vektor: 140.835529
Benoetigte Zeit Sor Vektor: 15.178807
%}

% ---------------------------------Anmerkungen---------------------------------
%{ 
  Sehe nicht wo der \ Operator benötigt wird.
  Alle Verfahren wurden auf zwei Arten implementiert.
  Zeilenweises Vorgehen ist unendlich langsam.
  Vektorielles Vorgehen hat etwas mehr Iterationen ist aber sehr viel schneller.
  Ausführen dauert bei derzeitigen Einstellungen ca. 20 min (intel core i5).
  Jacobi Vector ist verdammt schnell, obwohl es die meißten Iterationen benötigt.
  Sor ist eine Verbesserung von Gauss Seidel.
  Unter 'constants' können der Start- und Endwert für n gesetzt werden.
  Unter 'variables' kann die Variable 'vectorOnly' gesetzt werden. 
  'vectorOnly' bestimmt, ob die langsamen Implementationen mit genutzt werden.
  Wird 'vectorOnly' auf 'false' gesetzt, sollte 'ende' nicht viel größer als 10 
  gewählt werden, um tagelange Laufzeiten zu vermeiden.
  'start' muss größer gleich 2 gewählt werden und sollte kleiner als 'ende' sein.
%}