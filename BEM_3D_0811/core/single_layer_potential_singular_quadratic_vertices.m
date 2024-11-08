function [values] = single_layer_potential_singular_quadratic_vertices(area, ...
    triangle, index_singularity)

% SINGLE_LAYER_POTENTIAL_SINGULAR Computes the values on the diagonal of G
%
% INPUTS
% area: area of the triangle
% triangle: the coordinates of the points of the triangle
% index_singluarity: the local index of the point in which lies the
%                    singularity
%
% OUTPUTS
% values: the values of the integrals w.r.t. the three basis functions 

% temporary local relabeling of the vertexes

%%forse selezionerei nel triangolo solo i midpoints e fare relabelling su
%%di loro. (Fatto cosi in toeria converge a una soluzione non banale comunque (magari
%%non giusta).
% Nel caso comunque relabelling va fatto cosi (come era prima)
vertex_indexes = mod(index_singularity - 1 + [0 1 2],3) + 1;
vertexes = triangle(vertex_indexes,:);

midpoint_indexes = 3 + mod(index_singularity - 1 + [0 1 2],3) + 1;

% node_indexes = mod(index_singularity - 1 + [0 1 2 3 4 5],6) + 1;
% nodes = triangle(node_indexes,:);
% found = find(ismember(node_indexes,[1,2,3]));
% vertex_indexes = node_indexes(found);
% vertexes = triangle(vertex_indexes,:);

% integrand_1 = @(csi) (1-csi(1)-csi(2)).*(1-2*csi(1)-2*csi(2));
% integrand_2 = @(csi) csi(1).*(2*csi(1)-1);
% integrand_3 = @(csi) csi(2).*(2*csi(2)-1);
% integrand_4 = @(csi) 4*csi(1).*(1-csi(1)-csi(2));
% integrand_5 = @(csi) 4*csi(1).*csi(2);
% integrand_6 = @(csi) 4*csi(2).*(1-csi(1)-csi(2));

% coeff = area/(4*pi);
% upper_bound = @(csi) 1 - csi(1);

% defining integrands and coefficients
B = (vertexes(3,:) - vertexes(1,:)) * (vertexes(2,:) - vertexes(1,:))' ./ ...
    (norm((vertexes(2,:) - vertexes(1,:))));
C = (norm((vertexes(3,:) - vertexes(1,:)))) ./ ...
    (norm((vertexes(2,:) - vertexes(1,:))));
q_1= @(chi) -0.5*(1./(sin(chi)+cos(chi))) + ...
    2/3*(1./(sin(chi)+cos(chi))).^3.*(1+2*sin(chi).*cos(chi));
q_2= @(chi) 2/3*(1./(sin(chi)+cos(chi))).^3.*(cos(chi)).^2 - ...
    0.5*(1./(sin(chi)+cos(chi))).^2.*cos(chi);
q_3= @(chi) 2/3*(1./(sin(chi)+cos(chi))).^3.*(sin(chi)).^2 - ...
    0.5*(1./(sin(chi)+cos(chi))).^2.*sin(chi);
q_4= @(chi) 2*(1./(sin(chi)+cos(chi))).^2.*cos(chi) - ...
    4/3*(1./(sin(chi)+cos(chi))).^3.*((cos(chi)).^2+sin(chi).*cos(chi));
q_5= @(chi) 4/3*(1./(sin(chi)+cos(chi))).^3.*sin(chi).*cos(chi);
q_6= @(chi) 2*(1./(sin(chi)+cos(chi))).^2.*sin(chi) - ...
    4/3*(1./(sin(chi)+cos(chi))).^3.*((sin(chi)).^2+sin(chi).*cos(chi));
denominator = @(chi) sqrt(cos(chi).^2 + B*sin(2*chi) + C*sin(chi).^2);
integrand_1 = @(chi) q_1(chi) ./ denominator(chi);
integrand_2 = @(chi) q_2(chi) ./ denominator(chi);
integrand_3 = @(chi) q_3(chi) ./ denominator(chi);
integrand_4 = @(chi) q_4(chi) ./ denominator(chi);
integrand_5 = @(chi) q_5(chi) ./ denominator(chi);
integrand_6 = @(chi) q_6(chi) ./ denominator(chi);
coeff = 1/(4*pi) * 2 * area / norm((vertexes(2,:) - vertexes(1,:)));
%% quel 2*area non mi torna per niente perche' in teoria lui qui sta lavorando sul triangolo fisico
%% dovrebbe essere gia' giusto. (Non sta usando le coordinate sul triangolo di quadratura)



% computing integrals
values = zeros(1,6);
%%mancano i coefficienti di influenza sui midpoints ma per farlo serve
%%esprimere un punto del triangolo rispetto a 6 punti e non tre. (a meno di
%%spezzare l'integrale.

% values(node_indexes(1)) = coeff * integral(integrand_1,0,pi/2); % singolarità
% values(node_indexes(2)) = coeff * integral(integrand_2,0,pi/2); % secondo vertice
% values(node_indexes(3)) = coeff * integral(integrand_3,0,pi/2); % terzo vertice
% values(node_indexes(4)) = coeff * integral(integrand_4,0,pi/2); % primo mid
% values(node_indexes(5)) = coeff * integral(integrand_5,0,pi/2); % secondo mid
% values(node_indexes(6)) = coeff * integral(integrand_6,0,pi/2); % terzo mid

% relabeling
indexes = [vertex_indexes,midpoint_indexes];

values(indexes(1)) = coeff * integral(integrand_1,0,pi/2); % singolarità
values(indexes(2)) = coeff * integral(integrand_2,0,pi/2); % secondo vertice
values(indexes(3)) = coeff * integral(integrand_3,0,pi/2); % terzo vertice
values(indexes(4)) = coeff * integral(integrand_4,0,pi/2); % primo mid
values(indexes(5)) = coeff * integral(integrand_5,0,pi/2); % secondo mid
values(indexes(6)) = coeff * integral(integrand_6,0,pi/2); % terzo mid


%% qui mancano ovviamente tree assegnazioni, perche' il single layer potential va cambiato.

end