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

% ----------------------------- variables -------------------------------------
vectorOnly = true;    % changing this to false is not recommended
useParallel = true;   % only working on linux systems 
isParallelInstalled = false;  % if the parallel package is installed on the system
isStructInstalled = false;    % if the struct package is installed on the system

% ------------------------------ packages -------------------------------------

if useParallel
  %only working on linux systems
  if exist('OCTAVE_VERSION') ~= 0 && isunix ()
    instaled_packages = pkg ("list");
    for i=1:length(instaled_packages)
      if strcmp(instaled_packages{i}.name, "parallel")
        isParallelInstalled = true;
      end
      if strcmp(instaled_packages{i}.name, "struct")
        isStructInstalled = true;
      end
    end
    if !isParallelInstalled
      %{
      you'll need to run this once, to install the package
      if installation below is not working try following in bash console: 
        sudo apt-get install octave-<struct>
        sudo apt-get install octave-<parallel>
      %}
      if !isStructInstalled
        pkg install -forge struct;
      end
      pkg install -forge parallel;
      printf("parallel installed\n");
    else
      pkg load parallel; 
      printf("parallel loaded\n");
    end
  else
    useParallel = false;
  end
end 

% ---------------------------- function handles -------------------------------
f = @(x) 2./(1.+sin(pi./x)); %function to calculate omega(n)

% ------------------------------ constants ------------------------------------
eps = 1*10^(-6);  % error tolerance
maxit = 1*10^8;   % maximal number of iterations
start = 10;       % start number for n
ende = 60;        % end number for n
kerne = nproc-1;  % amount of used cores if 'useParallel' is 'true'

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

function [x, numit, t] = functionHandler(n,X,eps,maxit)
  [A,b] = my_test_system(n);
  x0 = b;
  omega = 2./(1.+sin(pi./n));
  tic;
  % switch on input
  switch (X)
    case "jacobi_vector"
      [x, numit] = my_jacobi_vector(A,b,x0,eps,maxit);
    case "gauss_seidel_vector"
      [x, numit] = my_gauss_seidel_vector(A,b,x0,eps,maxit);
    case "sor_vector"
      [x, numit] = my_sor_vector(A,b,x0,eps,maxit,omega);
    case "jacobi"
      [x, numit] = my_jacobi(A,b,x0,eps,maxit);
    case "gauss_seidel"
      [x, numit] = my_gauss_seidel(A,b,x0,eps,maxit);
    case "sor"
      [x, numit] = my_sor(A,b,x0,eps,maxit,omega);
    otherwise
      errordlg('wrong second argument!')
      x = -1;
      numit = -1;
      t = -1;
      return;
  endswitch   
  t = toc;
endfunction

% ------------------------------- main ----------------------------------------

if !useParallel
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
else
  % initialisierung
  n=start:ende;
  
  % Aufruf der Zeilenweisen Implementierung falls 'fast' == 'false' (way to slow)
  if !vectorOnly 
    [x_jacobi,numit_jacobi,t_jacobi] = pararrayfun(kerne, @(z) functionHandler(z,"jacobi",eps,maxit), n, "UniformOutput", false);
    t_jacobi = cell2mat(t_jacobi);
    numit_jacobi = cell2mat(numit_jacobi);
  
    [x_gauss_seidel,numit_gauss_seidel,t_gauss_seidel] = pararrayfun(kerne, @(z) functionHandler(z,"gauss_seidel",eps,maxit), n, "UniformOutput", false);
    t_gauss_seidel = cell2mat(t_gauss_seidel);
    numit_gauss_seidel = cell2mat(numit_gauss_seidel);
  
    [x_sor,numit_sor,t_sor] = pararrayfun(kerne, @(z) functionHandler(z,"sor",eps,maxit), n, "UniformOutput", false);
    t_sor = cell2mat(t_sor);
    numit_sor = cell2mat(numit_sor);
  end
  
  % Aufruf der Vektoriellen Implementierung
  [x_jacobi_vector,numit_jacobi_vector,t_jacobi_vector] = pararrayfun(kerne, @(z) functionHandler(z,"jacobi_vector",eps,maxit), n, "UniformOutput", false);
  t_jacobi_vector = cell2mat(t_jacobi_vector);
  numit_jacobi_vector = cell2mat(numit_jacobi_vector);
  
  [x_gauss_seidel_vector,numit_gauss_seidel_vector,t_gauss_seidel_vector] = pararrayfun(kerne, @(z) functionHandler(z,"gauss_seidel_vector",eps,maxit), n, "UniformOutput", false);
  t_gauss_seidel_vector = cell2mat(t_gauss_seidel_vector);
  numit_gauss_seidel_vector = cell2mat(numit_gauss_seidel_vector);
  
  [x_sor_vector,numit_sor_vector,t_sor_vector] = pararrayfun(kerne, @(z) functionHandler(z,"sor_vector",eps,maxit), n, "UniformOutput", false);
  t_sor_vector = cell2mat(t_sor_vector);
  numit_sor_vector = cell2mat(numit_sor_vector);
end

%output
printf("n = %1.0f:\n" , ende);
if !vectorOnly
  printf("Benoetigte Zeit Jacobi: %f\n" , t_jacobi(end));
  printf("Benoetigte Zeit Gauss Seidel: %f\n" , t_gauss_seidel(end));
  printf("Benoetigte Zeit Sor: %f\n" , t_sor(end));
end
printf("Benoetigte Zeit Jacobi Vektor: %f\n" , t_jacobi_vector(end));
printf("Benoetigte Zeit Gauss Seidel Vektor: %f\n" , t_gauss_seidel_vector(end));
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
hold off;
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
hold off;
xlabel("n");
ylabel("t");
xlim([n(1) n(end)]);
legend("show");
title("t in Abhaengigkeit von n");
print ("Fig2.jpg");

% ---------------------------------Ausgabe-------------------------------------
%{
parallel loaded
parcellfun: 51/51 jobs done
parcellfun: 51/51 jobs done
parcellfun: 51/51 jobs done
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
  Unter 'constants' können der Start- und Endwert für n gesetzt werden und der
  Parameter 'kerne'.
  'kerne' gibt an, wie viele Kerne zur Berechnung der Werte genutzt werden.
  Unter 'variables' können die Variable 'vectorOnly' und 'useParallel' gesetzt 
  werden. 
  'vectorOnly' bestimmt, ob die langsamen Implementationen mit genutzt werden,
  'useParallel' bestimmt, ob das 'parallel'-Package zur Nutzung mehrerer Kerne 
  benutzt werden soll. 
  'useParallel' funktioniert nur unter Linux Systemen oder Linux Subsystem für
  Windows 10.
  'useParallel' veringert die Gesamtlaufzeit um einen Faktor k, wobei k die Anzahl der
  benutzten Kerne ist (default ist k_max-1).
  Wird 'vectorOnly' auf 'false' gesetzt, sollte 'ende' nicht viel größer als 10 
  gewählt werden, um tagelange Laufzeiten zu vermeiden.
  'start' muss größer gleich 2 gewählt werden und sollte kleiner als 'ende' sein.
%}