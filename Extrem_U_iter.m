clear;

% Script para representar la energía del sistema frente al número de
% iteraciones.

% ------------------------------------
% Primer sistema.

load("1.5241U_iter.mat");
U_0_1 = [0, -1.5241e-17, -95.14, 0];
U_Sist1 = [U_0_1; U_iter];

plot(U_Sist1(:,1), U_Sist1(:,4), '-b', "LineWidth", 2.5,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('U - U_0 (eV)')
title ('U_0 = -95.14eV')
grid on
figure;

% ------------------------------------
% Segundo sistema.

load("1.5431U_iter.mat");
U_0_2 = [0, -1.5431e-17, -96.32, 0];
U_Sist2 = [U_0_2; U_iter];

plot(U_Sist2(:,1), U_Sist2(:,4), '-', "LineWidth", 1.5,"Marker",".","MarkerSize",20,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('U - U_0 (eV)')
title ('U_0 = -96.32eV')
grid on
figure;

% ------------------------------------
% Tercer sistema.

load("1.5544U_iter.mat");
U_0_3 = [0, -1.5544e-17, -97.03, 0];
U_Sist3 = [U_0_3; U_iter];

plot(U_Sist3(:,1), U_Sist3(:,4), '-', "LineWidth", 2.5,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('U - U_0 (eV)')
title ('U_0 = -97.03eV')
grid on
figure;

% ------------------------------------
% Cuarto sistema.

load("1.5706U_iter.mat");
U_0_4 = [0, -1.5706e-17, -98.04, 0];
U_Sist4 = [U_0_4; U_iter];

plot(U_Sist4(:,1), U_Sist4(:,4), '-', "LineWidth", 1.5,"Marker",".","MarkerSize",20,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('U - U_0 (eV)')
title ('U_0 = -98.04eV')
grid on
figure;

% ------------------------------------
% Quinto sistema.

load("1.5859U_iter.mat");
U_0_5 = [0, -1.5859e-17, -99.00, 0];
U_Sist5 = [U_0_5; U_iter];

plot(U_Sist5(:,1), U_Sist5(:,4), '-', "LineWidth", 1.5,"Marker",".","MarkerSize",20,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('U - U_0 (eV)')
title ('U_0 = -99.00eV')
grid on
figure;

% ------------------------------------

load("1.6037U_iter.mat");
U_0_6 = [0, -1.6037e-17, -100.11, 0];
U_Sist6 = [U_0_6; U_iter];
% load("1.4917U_iter.mat");
% U_0_6 = [0, -1.4917e-17, -93.11, 0];
% U_Sist6 = [U_0_6; U_iter];
plot(U_Sist1(:,1),U_Sist1(:,4),'b',"LineWidth", 2.5);
hold on
plot(U_Sist2(:,1),U_Sist2(:,4),'r',"LineWidth", 2.5);
hold on
plot(U_Sist3(:,1),U_Sist3(:,4),'g',"LineWidth", 2.5);
hold on
plot(U_Sist4(:,1),U_Sist4(:,4),'m',"LineWidth", 2.5);
hold on
plot(U_Sist5(:,1),U_Sist5(:,4),'c',"LineWidth", 2.5);
hold on
plot(U_Sist6(:,1),U_Sist6(:,4),'k',"LineWidth", 2.5);
hold off
set(gca, 'FontSize', 15);
xlabel ('Iteración')
ylabel ('U - U_0 (eV)')
legend ('-95.14 eV','-96.32 eV','-97.03 eV','-98.04 eV','-99.00 eV','Extremo')
lgd = legend;
lgd.FontSize = 20;
grid on;
figure;

% ------------------------------------

load("1.6037U_iter.mat");
U_0_6 = [0, -1.6037e-17, -100.11, 0];
U_Sist6 = [U_0_6; U_iter];
% load("1.4917U_iter.mat");
% U_0_6 = [0, -1.4917e-17, -93.11, 0];
% U_Sist6 = [U_0_6; U_iter];
plot(U_Sist1(:,1),U_Sist1(:,3),'b',"LineWidth", 2.5);
hold on
plot(U_Sist2(:,1),U_Sist2(:,3),'r',"LineWidth", 2.5);
hold on
plot(U_Sist3(:,1),U_Sist3(:,3),'g',"LineWidth", 2.5);
hold on
plot(U_Sist4(:,1),U_Sist4(:,3),'m',"LineWidth", 2.5);
hold on
plot(U_Sist5(:,1),U_Sist5(:,3),'c',"LineWidth", 2.5);
hold on
plot(U_Sist6(:,1),U_Sist6(:,3),'k',"LineWidth", 2.5);
hold off
set(gca, 'FontSize', 15);
xlabel ('Iteración')
ylabel ('U (eV)')
legend ('-95.14 eV','-96.32 eV','-97.03 eV','-98.04 eV','-99.00 eV','Extremo')
lgd = legend;
lgd.FontSize = 20;
grid on;
