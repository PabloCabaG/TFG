clear;

load("FWHM_Relj.mat");
N_estruct = FWHM_Relj(:,1) + 1;

plot(N_estruct, FWHM_Relj(:,2)*2.3548/sqrt(2)*180/pi, '-b', "LineWidth", 2.5)%,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Repeticiones')
ylabel ('FWHM (\circ)')
grid on

figure;

N_estruct_invers = N_estruct.^(-1);
t = (N_estruct)*3.08*4*1e-10;
t_invers = (t).^(-1);
beta_medio = FWHM_Relj(:,2)*2.3548/sqrt(2)*180/pi;
plot(t_invers, beta_medio, '-b', "LineWidth", 2.5,"Marker",".","MarkerSize",25,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15);
xlabel('1/t (m^{-1})')
ylabel('FWHM=\beta/2 (\circ)')
grid on

