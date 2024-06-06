clear;

% Script para representar el desplazamiento cuadrático medio del sistema frente al número de
% iteraciones.

% ------------------------------------
% Primer sistema.

load("1.5241Rho.mat");
Rho_Sist1 = rho_sat;

plot(Rho_Sist1(:,1), Rho_Sist1(:,2), ':k', "LineWidth", 2.5)%,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
hold on
plot(Rho_Sist1(:,1), Rho_Sist1(:,3), '-b', "LineWidth", 2,"Marker",".","MarkerSize",20,"MarkerEdgeColor",'r','MarkerFaceColor','r')
hold off
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('\rho (m)')
title ('U_0 = -95.14eV')
legend ('Límite superior de \rho', 'Valores medios de \rho')
lgd = legend;
lgd.FontSize = 20;
grid on
figure;

% ------------------------------------
% Segundo sistema.

load("1.5431Rho.mat");
Rho_Sist2 = rho_sat;

plot(Rho_Sist2(:,1), Rho_Sist2(:,2), ':k', "LineWidth", 2.5)%,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
hold on
plot(Rho_Sist2(:,1), Rho_Sist2(:,3), '-b', "LineWidth", 2,"Marker",".","MarkerSize",20,"MarkerEdgeColor",'r','MarkerFaceColor','r')
hold off
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('\rho (m)')
title ('U_0 = -96.32eV')
legend ('Límite superior de \rho', 'Valores medios de \rho')
lgd = legend;
lgd.FontSize = 20;
grid on
figure;

% ------------------------------------
% Tercer sistema.

load("1.5544Rho.mat");
Rho_Sist3 = rho_sat;

plot(Rho_Sist3(:,1), Rho_Sist3(:,2), ':k', "LineWidth", 2.5)%,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
hold on
plot(Rho_Sist3(:,1), Rho_Sist3(:,3), '-b', "LineWidth", 2,"Marker",".","MarkerSize",20,"MarkerEdgeColor",'r','MarkerFaceColor','r')
hold off
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('\rho (m)')
title ('U_0 = -97.03eV')
legend ('Límite superior de \rho', 'Valores medios de \rho')
lgd = legend;
lgd.FontSize = 20;
grid on
figure;

% ------------------------------------
% Cuarto sistema.

load("1.5706Rho.mat");
Rho_Sist4 = rho_sat;

plot(Rho_Sist4(:,1), Rho_Sist4(:,2), ':k', "LineWidth", 2.5)%,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
hold on
plot(Rho_Sist4(:,1), Rho_Sist4(:,3), '-b', "LineWidth", 2,"Marker",".","MarkerSize",20,"MarkerEdgeColor",'r','MarkerFaceColor','r')
hold off
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('\rho (m)')
title ('U_0 = -98.04eV')
legend ('Límite superior de \rho', 'Valores medios de \rho')
lgd = legend;
lgd.FontSize = 20;
grid on
figure;

% ------------------------------------
% Quinto sistema.

load("1.5859Rho.mat");
Rho_Sist5 = rho_sat;

plot(Rho_Sist5(:,1), Rho_Sist5(:,2), ':k', "LineWidth", 2.5)%,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
hold on
plot(Rho_Sist5(:,1), Rho_Sist5(:,3), '-b', "LineWidth", 2,"Marker",".","MarkerSize",20,"MarkerEdgeColor",'r','MarkerFaceColor','r')
hold off
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Iteración')
ylabel ('\rho (m)')
title ('U_0 = -99.00eV')
legend ('Límite superior de \rho', 'Valores medios de \rho')
lgd = legend;
lgd.FontSize = 20;
grid on
figure;

% ------------------------------------
plot(Rho_Sist1(:,1),Rho_Sist1(:,3),'b',"LineWidth", 2);
hold on
plot(Rho_Sist2(:,1),Rho_Sist2(:,3),'r',"LineWidth", 2);
hold on
plot(Rho_Sist3(:,1),Rho_Sist3(:,3),'g',"LineWidth", 2);
hold on
plot(Rho_Sist4(:,1),Rho_Sist4(:,3),'m',"LineWidth", 2);
hold on
plot(Rho_Sist5(:,1),Rho_Sist5(:,3),'c',"LineWidth", 2);
hold on
plot(Rho_Sist5(:,1),Rho_Sist5(:,2),':k',"LineWidth", 2.5);
hold off
set(gca, 'FontSize', 15);
xlabel ('Iteración')
ylabel ('\rho (m)')
legend ('-95.14 eV','-96.32 eV','-97.03 eV','-98.04 eV','-99.00 eV','Límite superior')
lgd = legend;
lgd.FontSize = 20;
grid on;
