function [values] = single_layer_potential_singular_quadratic(area, ...
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
% indexes = mod(index_singularity - 1 + [0 1 2],3) + 1;

indexes = mod(index_singularity - 1 + [0 1 2 3 4 5],6) + 1;

nodes = triangle(indexes,:);

integrand_1 = @(csi) 1-csi(1)-csi(2)).*(1-2*csi(1)-2*csi(2));
integrand_2 = @(csi) csi(1).*(2*csi(1)-1);
integrand_3 = @(csi) csi(2).*(2*csi(2)-1);
integrand_4 = @(csi) 4*csi(1).*(1-csi(1)-csi(2));
integrand_5 = @(csi) 4*csi(1).*csi(2);
integrand_6 = @(csi) 4*csi(2).*(1-csi(1)-csi(2));

coeff = area/(4*pi);

upper_bound = @(csi) 1 - csi(1);

% defining integrands and coefficients
B = (nodes(3,:) - nodes(1,:)) * (nodes(2,:) - nodes(1,:))' ./ ...
    (norm((nodes(2,:) - nodes(1,:))));
C = (norm((nodes(3,:) - nodes(1,:)))) ./ ...
    (norm((nodes(2,:) - nodes(1,:))));
q_1= @(chi) 1 ./ (cos(chi) + sin(chi))  - ...
    (sin(chi) + cos(chi)) ./ (2 * (cos(chi) + sin(chi)).^2);
q_2= @(chi) cos(chi) ./ (2 * (cos(chi) + sin(chi)).^2);
q_3= @(chi) sin(chi) ./ (2 * (cos(chi) + sin(chi)).^2);
denominator = @(chi) sqrt(cos(chi).^2 + B*sin(2*chi) + C*sin(chi).^2);
integrand_1 = @(chi) q_1(chi) ./ denominator(chi);
integrand_2 = @(chi) q_2(chi) ./ denominator(chi);
integrand_3 = @(chi) q_3(chi) ./ denominator(chi);
coeff = 1/(4*pi) * 2 * area / norm((nodes(2,:) - nodes(1,:)));
%% quel 2*area non mi torna per niente perche' in teoria lui qui sta lavorando sul triangolo fisico
%% dovrebbe essere gia' giusto. (Non sta usando le coordinate sul triangolo di quadratura)



% computing integrals
values = zeros(1,6);
%%mancano i coefficienti di influenza sui midpoints ma per farlo serve
%%esprimere un punto del triangolo rispetto a 6 punti e non tre. (a meno di
%%spezzare l'integrale.


values(indexes(1)) = coeff * integral(integrand_1,0,pi/2);
values(indexes(2)) = coeff * integral(integrand_2,0,pi/2);
values(indexes(3)) = coeff * integral(integrand_3,0,pi/2);
%% qui mancano ovviamente tree assegnazioni, perche' il single layer potential va cambiato.

end

