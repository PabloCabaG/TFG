clear;

% Representación del número de publicaciones de artículos de HEAs por año.

N_public = [32, 35, 30, 56, 71, 92, 90, 111, 104, 202, 252, 373, 476, 651, 980, 1266, 1717, 2245, 2970, 3373];
Anio = [2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023];

bar(Anio, N_public, 0.7, 'r')
hold on;
plot(Anio, N_public, '-.ok', 'LineWidth', 2,'MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','b')
set(gca, 'FontSize', 20);
xlabel('Año')
ylabel('Número de publicaciones')
set(gca, 'XTick', Anio);

% for i = 1:length(N_public)
%     text(Anio(i), N_public(i), num2str(N_public(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 15)
% end