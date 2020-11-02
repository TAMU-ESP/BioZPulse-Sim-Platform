function plot_3d_grid_sides_fn(A,opt,clims)
% specifies all the vertices that comprises the object you want to draw

[Nx,Ny,Nz]=size(A);
vert = [-0.5 -0.5 -0.5;  ...
    -0.5 0.5 -0.5;  ...
    0.5 0.5 -0.5;  ...
    0.5 -0.5 -0.5; ...
    -0.5 -0.5 0.5; ...
    -0.5 0.5 0.5;  ...
    0.5 0.5 0.5; ...
    0.5 -0.5 0.5];
% define the arbitrary polygon(patch) using the vertice number(index) you defined above.
fac = [1 2 3 4; ...
    2 6 7 3; ...
    4 3 7 8; ...
    1 5 8 4; ...
    1 2 6 5; ...
    5 6 7 8];
% specify patch (polygons) in patch() function
% just call the patch function multiple times to draw multiple cubes
if isempty(clims)
    cmin = min(A(~isnan(A)&~isinf(A)));
    cmax = max(A(~isnan(A)&~isinf(A)));
    clims=[cmin,cmax];
    if cmin==cmax
        cmin=cmin/2;
        cmax=cmax*2;
        clims=[cmin,cmax];
    end
else
    cmin=clims(1);
    cmax=clims(2);
end

colormap jet

cmap = colormap;
m = length(cmap);
x_arr=1:Nx;y_arr=1:Ny;z_arr=1:Nz;

A2d=squeeze(A(x_arr(end),y_arr,z_arr));
im=imagesc(A2d,clims);
set(im,'Visible','off')
colorbar
index = fix((A2d-cmin)/(cmax-cmin)*m)+1; %A
RGB = ind2rgb(index,cmap);
RGB=rgb_edit_fn(RGB,A2d);
for x=x_arr(end)
    for y=y_arr
        for z=z_arr
            shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
            fig=patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',RGB(y,z,:));  % draw the blue cube
        end
    end
end
hold on;

A2d=squeeze(A(x_arr(1),y_arr,z_arr));
index = fix((A2d-cmin)/(cmax-cmin)*m)+1; %A
RGB = ind2rgb(index,cmap);
RGB=rgb_edit_fn(RGB,A2d);
for x=x_arr(1)
    for y=y_arr
        for z=z_arr
            shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
            patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',RGB(y,z,:));  % draw the blue cube
        end
    end
end

A2d=squeeze(A(x_arr,y_arr(1),z_arr));
index = fix((A2d-cmin)/(cmax-cmin)*m)+1; %A
RGB = ind2rgb(index,cmap);
RGB=rgb_edit_fn(RGB,A2d);
for x=x_arr
    for y=y_arr(1)
        for z=z_arr
            shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
            patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',RGB(x,z,:));  % draw the blue cube
        end
    end
end

A2d=squeeze(A(x_arr,y_arr(end),z_arr));
index = fix((A2d-cmin)/(cmax-cmin)*m)+1; %A
RGB = ind2rgb(index,cmap);
RGB=rgb_edit_fn(RGB,A2d);
for x=x_arr
    for y=y_arr(end)
        for z=z_arr
            shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
            patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',RGB(x,z,:));  % draw the blue cube
        end
    end
end

A2d=squeeze(A(x_arr,y_arr,z_arr(end)));
index = fix((A2d-cmin)/(cmax-cmin)*m)+1; %A
RGB = ind2rgb(index,cmap);
RGB=rgb_edit_fn(RGB,A2d);
for x=x_arr
    for y=y_arr
        for z=z_arr(end)
            shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
            patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',RGB(x,y,:));  % draw the blue cube
        end
    end
end

A2d=squeeze(A(x_arr,y_arr,z_arr(1)));
index = fix((A2d-cmin)/(cmax-cmin)*m)+1; %A
RGB = ind2rgb(index,cmap);
RGB=rgb_edit_fn(RGB,A2d);
for x=x_arr
    for y=y_arr
        for z=z_arr(1)
            shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
            patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',RGB(x,y,:));  % draw the blue cube
        end
    end
end

axis([0 max(x_arr)+1 0 max(y_arr)+1 0 max(z_arr)+1]);
grid();
view(30,30);

%set(gca,'DataAspectRatio',[Nx Nz Ny])
%material shiny;
%alpha('color');
%alphamap('rampdown');


