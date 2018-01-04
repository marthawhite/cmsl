function plot_compare_cca(X1, Y1, phi, Ccca, Cours, graph_title);

% Normalize Cours to same length as Ccca to facilitate visualization
Cours = Cours * (norm(Ccca)/norm(Cours));
%norm(Cours)
%Cours = Cours / norm(Cours);
%Ccca = Ccca / norm(Ccca);
Ccca
Cours

figure

plot3(X1(:,1),X1(:,2),Y1,'b.');

hold on;

plot_direction(Ccca, phi, '--r*');
plot_direction(Cours, phi, '-gd');

title(graph_title, 'FontWeight','bold','FontSize',24);
legend({'Data', 'CCA', 'Opponent'});
hold off;

end

