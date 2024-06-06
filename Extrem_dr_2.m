clear;

% Script para representar el desplazamiento cuadrático medio del sistema frente al número de
% iteraciones.

% ------------------------------------
% Primer sistema.

load("1.5241dr_2.mat");
dr2_0_1 = [0, 0];
dr2_Sist1 = [dr2_0_1; dr_2];

plot(dr2_Sist1(:,1), dr2_Sist1(:,2) * sqrt(3)/2, '-b', "LineWidth", 2.5,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('\deltar^2 (m^2/{\ita}^2)')
title ('U_0 = -95.14eV')
grid on
figure;

% ------------------------------------
% Segundo sistema.

load("1.5431dr_2.mat");
dr2_0_2 = [0, 0];
dr2_Sist2 = [dr2_0_2; dr_2];

plot(dr2_Sist2(:,1), dr2_Sist2(:,2) * sqrt(3)/2, '-', "LineWidth", 1.5,"Marker",".","MarkerSize",20,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('\deltar^2 (m^2/{\ita}^2)')
title ('U_0 = -96.32eV')
grid on
figure;

% ------------------------------------
% Tercer sistema.

load("1.5544dr_2.mat");
dr2_0_3 = [0, 0];
dr2_Sist3 = [dr2_0_3; dr_2];

plot(dr2_Sist3(:,1), dr2_Sist3(:,2) * sqrt(3)/2, '-', "LineWidth", 1.5,"Marker",".","MarkerSize",20,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('\deltar^2 (m^2/{\ita}^2)')
title ('U_0 = -97.03eV')
grid on
figure;

% ------------------------------------
% Cuarto sistema.

load("1.5706dr_2.mat");
dr2_0_4 = [0, 0];
dr2_Sist4 = [dr2_0_4; dr_2];

plot(dr2_Sist4(:,1), dr2_Sist4(:,2) * sqrt(3)/2, '-', "LineWidth", 1.5,"Marker",".","MarkerSize",20,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('\deltar^2 (m^2/{\ita}^2)')
title ('U_0 = -98.04eV')
grid on
figure;

% ------------------------------------
% Quinto sistema.

load("1.5859dr_2.mat");
dr2_0_5 = [0, 0];
dr2_Sist5 = [dr2_0_5; dr_2];

plot(dr2_Sist5(:,1), dr2_Sist5(:,2) * sqrt(3)/2, '-', "LineWidth", 1.5,"Marker",".","MarkerSize",20,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('\deltar^2 (m^2/{\ita}^2)')
title ('U_0 = -99.00eV')
grid on
figure;

% ------------------------------------
% load("1.6037dr_2.mat");
load("1.4917dr_2.mat");
dr2_0_6 = [0, 0];
dr2_Sist6 = [dr2_0_6; dr_2];
plot(dr2_Sist1(:,1),dr2_Sist1(:,2) * sqrt(3)/2,'b',"LineWidth", 2.5);
hold on
plot(dr2_Sist2(:,1),dr2_Sist2(:,2) * sqrt(3)/2,'r',"LineWidth", 2.5);
hold on
plot(dr2_Sist3(:,1),dr2_Sist3(:,2) * sqrt(3)/2,'g',"LineWidth", 2.5);
hold on
plot(dr2_Sist4(:,1),dr2_Sist4(:,2) * sqrt(3)/2,'m',"LineWidth", 2.5);
hold on
plot(dr2_Sist5(:,1),dr2_Sist5(:,2) * sqrt(3)/2,'c',"LineWidth", 2.5);
hold on
plot(dr2_Sist6(:,1),dr2_Sist6(:,2) * sqrt(3)/2,'k',"LineWidth", 2.5);
hold off
set(gca, 'FontSize', 15);
xlabel ('Iteración')
ylabel ('\deltar^2 (m^2/{\ita}^2)')
legend ('-95.14 eV','-96.32 eV','-97.03 eV','-98.04 eV','-99.00 eV','Extremo')
lgd = legend;
lgd.FontSize = 20;
grid on;
