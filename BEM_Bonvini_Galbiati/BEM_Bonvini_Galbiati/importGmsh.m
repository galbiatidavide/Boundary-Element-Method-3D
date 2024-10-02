function mesh = importGmsh(file_name,i)
% Modified from "Sathyanarayan Rao (2021). Visualizing GMSH .msh file in
% MATLAB using TRIPLOT
% (https://www.mathworks.com/matlabcentral/fileexchange/59682-visualizing-gmsh-msh-file-in-matlab-using-triplot),
% MATLAB Central File Exchange. Retrieved June 30, 2021."
%------------------------------------------------------------------------%
%------ Gmsh to Matlab script: Import mesh to matlab---------------------%
%------------------------------------------------------------------------%
%-----------------------------------------------------------------------%
% dlmread(filename,delimiter,[R1 C1 R2 C2]) reads only the range
% bounded by row offsets R1 and R2 and column offsets C1 and C2.
%-----------------------------------------------------------------------%
file    =  (file_name);
if i==2
    % no of nodes is mentioned in 5th row and first column
    
    N_n      = dlmread(file,'',[5-1 1-1 5-1 1-1]);
    N_e      = dlmread(file,'',[7+N_n 0 7+N_n 0]);
    
    nodes       = dlmread(file,'',[5 1 4+N_n 3]);
    elements    = dlmread(file,'',[8+N_n 0 7+N_n+N_e 7]);
    
    %------- 2D Geometry
    
    two_d_nodes = nodes(:,1:2);
    elem_type   = elements(:,2);
    
    %--- find the starting indices of 2D elements
    two_ind=sum(elem_type~=2)+1;
    %----------------------------------------------
    
    two_d_elements(1:N_e-two_ind,1:3) = 0;
    k = 1;
    for i = two_ind:N_e
        two_d_elements(k,1:3) = elements(i,6:8);
        k = k+1;
    end
    
    boundary_temp=elements(elements(:,2)==1,5:7);
    boundary_length=zeros(max(boundary_temp(:,1)),1);
    for n=1:max(boundary_temp(:,1))
        boundary=reshape(boundary_temp(boundary_temp(:,1)==n,2:3)',[],1)';
        boundary=unique(boundary,'stable');
        boundary_length(n,1)=size(boundary,2);
    end
    
    boundary=zeros(max(boundary_temp(:,1)),max(boundary_length));
    
    for n=1:max(boundary_temp(:,1))
        boundary(n,3:boundary_length(n,1)+2)=unique(reshape(boundary_temp(boundary_temp(:,1)==n,2:3)',[],1)','stable');
        boundary(n,1)=n;
        boundary(n,2)=boundary_length(n,1);
    end
    
    mesh.Nodes=two_d_nodes;
    mesh.Elements=two_d_elements;
    mesh.Boundary=boundary;
    
    %---- visualize in matlab ---------------------
    
    figure(1)
    triplot(two_d_elements,two_d_nodes(:,1),two_d_nodes(:,2))
    xlabel('X','fontsize',14)
    ylabel('Y','fontsize',14)
    title('Gmsh to MATLAB import','fontsize',14)
    fh = figure(1);
    set(fh, 'color', 'white');
    axis equal
    
    %-------------------------------------------------------------------------
elseif i==1  %quadratic element
    % no of nodes is mentioned in 5th row and first column
    
    N_n      = dlmread(file,'',[5-1 1-1 5-1 1-1]);
    N_e      = dlmread(file,'',[7+N_n 0 7+N_n 0]);
    
    nodes       = dlmread(file,'',[5 1 4+N_n 3]);
    elements    = dlmread(file,'',[8+N_n 0 7+N_n+N_e 10]);
    
    %------- 2D Geometry
    
    two_d_nodes = nodes(:,1:2);
    elem_type   = elements(:,2);
    
    %--- find the starting indices of 2D elements
    two_ind=sum(elem_type~=9)+1;
    %----------------------------------------------
    
    two_d_elements(1:N_e-two_ind,1:6) = 0;
    k = 1;
    for i = two_ind:N_e
        two_d_elements(k,1:6) = elements(i,6:11);
        k = k+1;
    end
    
    %two_d_elements=two_d_elements(:,[1 4 2 5 3 6]);
    
    boundary_temp=elements(elements(:,2)==8,5:8);
    boundary_temp=boundary_temp(:,[1 2 4 3]);
    boundary_length=zeros(max(boundary_temp(:,1)),1);
    for n=1:max(boundary_temp(:,1))
        boundary=reshape(boundary_temp(boundary_temp(:,1)==n,2:4)',[],1)';
        boundary=unique(boundary,'stable');
        boundary_length(n,1)=size(boundary,2);
    end
    
    boundary=zeros(max(boundary_temp(:,1)),max(boundary_length));
    
    for n=1:max(boundary_temp(:,1))
        boundary(n,3:boundary_length(n,1)+2)=unique(reshape(boundary_temp(boundary_temp(:,1)==n,2:4)',[],1)','stable');
        boundary(n,1)=n;
        boundary(n,2)=boundary_length(n,1);
    end
    
    mesh.Nodes=two_d_nodes;
    mesh.Elements=two_d_elements;
    mesh.Boundary=boundary;
    
    %---- visualize in matlab ---------------------
    
    figure(1)
    trimesh(two_d_elements,two_d_nodes(:,1),two_d_nodes(:,2),zeros(size(two_d_nodes(:,2),1),1))
    xlabel('X','fontsize',14)
    ylabel('Y','fontsize',14)
    title('Gmsh to MATLAB import','fontsize',14)
    fh = figure(1);
    set(fh, 'color', 'white');
    view(2)
    axis equal
    
    %-------------------------------------------------------------------------
    
    
else
    fprintf('Please specify whether the element is linear (2) or quadratic (1)')
end
