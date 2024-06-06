clear;

% Representación de la anchura a altura media del pico de difracción frente al número de veces que
% hemos repetido la estructura del sólido en x e y.

load("FWHM_Rep_SnRelj.mat")

plot(FWHM_frente_rep(:,1), FWHM_frente_rep(:,2)/sqrt(2)*180/pi, '-b', "LineWidth", 2.5)%,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Repeticiones')
ylabel ('FWHM (\circ)')
% title ('Sin relajar, U_0 = -99.00eV')
grid on