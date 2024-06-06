
% Representación del pico de difracción.

load("1.5859Relj_30rep.mat");
load("1.5859SnRelj_30rep.mat");

plot(Matriz_0(:,2) * 180 / pi, Matriz_0(:,4), '-b', "LineWidth", 1.5)%,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('\Theta (\circ)')
ylabel ('|F|^2')
title ('Sin relajar, U_0 = -99.00eV')
grid on

figure;

plot(Matriz_relj(:,2) * 180 / pi, Matriz_relj(:,4), '-b', "LineWidth", 1.5)%,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('\Theta (\circ)')
ylabel ('|F|^2')
title ('Relajada, U_0 = -99.00eV')
grid on

figure;

plot(Matriz_0(:,2) * 180 / pi, Matriz_0(:,4), '-b', "LineWidth", 1.5)
hold on
plot(Matriz_relj(:,2) * 180 / pi, Matriz_relj(:,4), '-r', "LineWidth", 1.5)
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('\Theta (\circ)')
ylabel ('|F|^2')
title ('U_0 = -99.00eV')
legend (' Sin relajar',' Relajada')
lgd = legend;
lgd.FontSize = 20;
grid on