function [values] = single_layer_potential_singular_quadratic_midpoints(area, ...
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

midpoint_indexes = 3 + mod(index_singularity - 1 + [0 1 2],3) + 1;
midpoints = triangle(midpoint_indexes,:);

% integrand_1 = @(csi) (1-csi(1)-csi(2)).*(1-2*csi(1)-2*csi(2));
% integrand_2 = @(csi) csi(1).*(2*csi(1)-1);
% integrand_3 = @(csi) csi(2).*(2*csi(2)-1);
% integrand_4 = @(csi) 4*csi(1).*(1-csi(1)-csi(2));
% integrand_5 = @(csi) 4*csi(1).*csi(2);
% integrand_6 = @(csi) 4*csi(2).*(1-csi(1)-csi(2));

% coeff = area/(4*pi);
% upper_bound = @(csi) 1 - csi(1);

% defining integrands and coefficients
B = (midpoints(3,:) - midpoints(1,:)) * (midpoints(2,:) - midpoints(1,:))' ./ ...
    (norm((midpoints(2,:) - midpoints(1,:))));
C = (norm((midpoints(3,:) - midpoints(1,:)))) ./ ...
    (norm((midpoints(2,:) - midpoints(1,:))));
q_11 = @(chi) -sec(chi)/12;
q_12 = @(chi) 1/12*cot(chi).*(-3 + 2*cot(chi)).*csc(chi);
q_21 = @(chi) 1/12*sec(chi).*tan(chi).*(-3 + 2*tan(chi));
q_22 = @(chi) -csc(chi)/12;
q_31 = @(chi) 1/12*sec(chi).*(-1 + tan(chi) + 2*tan(chi).^2);
q_32 = @(chi) 1/24*csc(chi).^3.*(1 + 3*cos(2*chi) + sin(2*chi));
q_41 = @(chi) -1/6*sec(chi).*(-3 + tan(chi));
q_42 = @(chi) -1/6*(-3 + cot(chi)).*csc(chi);
q_51 = @(chi) 1/12*sec(chi).^3.*(1 + 5*cos(2*chi) + sin(2*chi));
q_52 = @(chi) 1/6*(1 + cot(chi)).*csc(chi);
q_61 = @(chi) 1/6*sec(chi).*(1 + tan(chi));
q_62 = @(chi) 1/6*(3 + cot(chi) - 2*cot(chi).^2).*csc(chi);


denominator = @(chi) 4*sqrt(cos(chi).^2 + B*sin(2*chi) + C*sin(chi).^2);
integrand_11 = @(chi) q_11(chi) ./ denominator(chi);
integrand_12 = @(chi) q_12(chi) ./ denominator(chi);
integrand_21 = @(chi) q_21(chi) ./ denominator(chi);
integrand_22 = @(chi) q_22(chi) ./ denominator(chi);
integrand_31 = @(chi) q_31(chi) ./ denominator(chi);
integrand_32 = @(chi) q_32(chi) ./ denominator(chi);
integrand_41 = @(chi) q_41(chi) ./ denominator(chi);
integrand_42 = @(chi) q_42(chi) ./ denominator(chi);
integrand_51 = @(chi) q_51(chi) ./ denominator(chi);
integrand_52 = @(chi) q_52(chi) ./ denominator(chi);
integrand_61 = @(chi) q_61(chi) ./ denominator(chi);
integrand_62 = @(chi) q_62(chi) ./ denominator(chi);
coeff = 1/(4*pi) * 2 * area / norm((midpoints(2,:) - midpoints(1,:)));

%% quel 2*area non mi torna per niente perche' in teoria lui qui sta lavorando sul triangolo fisico
%% dovrebbe essere gia' giusto. (Non sta usando le coordinate sul triangolo di quadratura)



% computing integrals
values = zeros(1,6);
%%mancano i coefficienti di influenza sui midpoints ma per farlo serve
%%esprimere un punto del triangolo rispetto a 6 punti e non tre. (a meno di
%%spezzare l'integrale.

% relabeling
indexes = [vertex_indexes,midpoint_indexes];

values(indexes(1)) = coeff * (integral(integrand_11,-pi/4,pi/4)+integral(integrand_12,pi/4,3*pi/4));
values(indexes(2)) = coeff * (integral(integrand_21,-pi/4,pi/4)+integral(integrand_22,pi/4,3*pi/4));
values(indexes(3)) = coeff * (integral(integrand_31,-pi/4,pi/4)+integral(integrand_32,pi/4,3*pi/4));
values(indexes(4)) = coeff * (integral(integrand_41,-pi/4,pi/4)+integral(integrand_42,pi/4,3*pi/4));
values(indexes(5)) = coeff * (integral(integrand_51,-pi/4,pi/4)+integral(integrand_52,pi/4,3*pi/4));
values(indexes(6)) = coeff * (integral(integrand_61,-pi/4,pi/4)+integral(integrand_62,pi/4,3*pi/4));

%% qui mancano ovviamente tree assegnazioni, perche' il single layer potential va cambiato.

end