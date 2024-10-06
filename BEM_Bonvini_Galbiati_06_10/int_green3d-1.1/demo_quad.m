clear variables, close all

%% Calculate the integral of 1/r over the square [0 0 0; 1 0 0; 1 1 0; 0 1 0].
% The integration points belongs to a line orthogonal to the square,
% passing through the square barycenter

% definition of square: note that the vertice coordinates are in 3d
V = [0 0 0; 1 0 0; 1 1 0; 0 1 0];

% get the triangle barycenter
B = mean(V,1);

% number of field points. Let's keep big to measure the CPU time
npt = 5001;

% define the varying z coordinate
z = linspace(-1,1,npt);

% define the field points
Q = [repmat(B(1:2),npt,1) z(:)];

% plot the triangle and the field points
figure; hold all
patch('Faces',[1 2 3 4],'Vertices',V,'FaceColor','r','EdgeColor','k','FaceAlpha',.5);
plot3(Q(:,1),Q(:,2),Q(:,3),'bo','MarkerFaceColor','b');
view([1 1 1])
axis equal
title('geometry')

% Now start the integration. INT_GREEN3D_TRI is used twice exploiting the
% integral additivity

% integration over two triangle
t0 = tic;
[Itri1,Igradtri1] = int_green3d_tri(Q,V([1 2 3],:));
[Itri2,Igradtri2] = int_green3d_tri(Q,V([1 3 4],:));
Itri = Itri1+Itri2;
Igradtri = Igradtri1+Igradtri2;
ttri = toc(t0);

% integration over a polygon
t0 = tic;
[Ipoly,Igradpoly] = int_green3d_poly(Q,V);
tpoly = toc(t0);

% Use the Matlab built-in numerical quadrature QUAD2D
Inum = zeros(npt,1);
t0 = tic;
for i = 1:npt
    fun = @(x,y)1./sqrt((x-Q(i,1)).^2+(y-Q(i,2)).^2+Q(i,3).^2);
    Inum(i) = quad2d(fun,0,1,0,1);
end
tnum = toc(t0);

% Use the Matlab built-in numerical quadrature QUAD2D
Igradnum = zeros(npt,3);
t0 = tic;
for i = 1:npt
    funx = @(x,y)(x-Q(i,1))./((x-Q(i,1)).^2+(y-Q(i,2)).^2+Q(i,3).^2).^1.5;
    funy = @(x,y)(y-Q(i,2))./((x-Q(i,1)).^2+(y-Q(i,2)).^2+Q(i,3).^2).^1.5;
    funz = @(x,y)-Q(i,3)./((x-Q(i,1)).^2+(y-Q(i,2)).^2+Q(i,3).^2).^1.5;
    Igradnum(i,1) = quad2d(funx,0,1,0,1,'Singular',true);
    Igradnum(i,2) = quad2d(funy,0,1,0,1,'Singular',true);
    Igradnum(i,3) = quad2d(funz,0,1,0,1,'Singular',true);
end
tnum2 = toc(t0);

% Now compare results
figure, hold all
plot(z,Itri);
plot(z,Ipoly);
plot(z,Inum);
xlabel('z-coordinate')
ylabel('Integral of 1/r');
legend('Analytic triangle','Analytic polygon','Numerical');

% Now compare results
figure
subplot(3,1,1), hold all
plot(z,Igradtri(:,1));
plot(z,Igradpoly(:,1));
plot(z,Igradnum(:,1));
xlabel('z-coordinate')
ylabel('Integral of grad(1/r)_x');
legend('Analytic triangle','Analytic polygon','Numerical');

subplot(3,1,2), hold all
plot(z,Igradtri(:,2));
plot(z,Igradpoly(:,2));
plot(z,Igradnum(:,2));
xlabel('z-coordinate')
ylabel('Integral of grad(1/r)_y');
legend('Analytic triangle','Analytic polygon','Numerical');

subplot(3,1,3), hold all
plot(z,Igradtri(:,3));
plot(z,Igradpoly(:,3));
plot(z,Igradnum(:,3));
xlabel('z-coordinate')
ylabel('Integral of grad(1/r)_z');
legend('Analytic triangle','Analytic polygon','Numerical');

% print timing
fprintf('CPU time (%d points):\n',npt)
fprintf('* Analytic triangle: %f sec\n',ttri);
fprintf('* Analytic polygon: %f sec\n',tpoly);
fprintf('* Numerical quad2d 1/r: %f sec\n',tnum);
fprintf('* Numerical quad2d grad(1/r): %f sec\n',tnum2);