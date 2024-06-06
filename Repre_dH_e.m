clear;

% Representación de la deformación random 0 : 0.07 : 0.49.

load("dH_frente_e.mat")

plot(dH_frente_e(:,2),dH_frente_e(:,4)./dH_frente_e(:,3), '-b', "LineWidth", 2.5,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('\deltar^2')
ylabel ('H/H_0')
grid on

figure;

plot(dH_frente_e(:,2),dH_frente_e(:,5), '-b', "LineWidth", 2.5,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('\deltar^2')
ylabel ('\DeltaH/H_0')
grid on