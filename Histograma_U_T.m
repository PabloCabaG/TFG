clear;

% Código para generar un histograma con los valores del scrip BCC_sin_Iter.
% Primero hay que correr el código anterior para crear la Matriz_U_T.

load("1000U_0.mat")

histogram(Matriz_U_T(:,2)/1.602e-19,'NumBins', 10,'FaceColor', 'r', 'Normalization', 'probability');
set(gca, 'FontSize', 15);
xlabel('U_0 (eV)', 'FontSize', 20)
ylabel('Proporción','FontSize', 20)
title ('Distribución de la energía inicial para distintas configuraciones aleatorias', 'FontSize', 20)