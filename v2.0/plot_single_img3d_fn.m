function plot_single_img3d_fn(net,opt,clims)
%opt=1 --> cond

%A=permute(net.img_arr,[3,1,2]);
%A=(net.img_arr);
A=(net.ri_img);
figure
if opt==1
    inf_ind=(A==inf);
    inf_neg_ind=(A==-inf);
    A=1./A;
    A(inf_ind)=inf;
    A(inf_neg_ind)=-inf;
end

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
if size(A2d,1)>1 && size(A2d,2)>1
    for x=x_arr(end)
        for y=y_arr
            for z=z_arr
                shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
                fig=patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',RGB(y,z,:));  % draw the blue cube
            end
        end
    end
end
hold on;

A2d=squeeze(A(x_arr(1),y_arr,z_arr));
index = fix((A2d-cmin)/(cmax-cmin)*m)+1; %A
RGB = ind2rgb(index,cmap);
RGB=rgb_edit_fn(RGB,A2d);
if size(A2d,1)>1 && size(A2d,2)>1
    for x=x_arr(1)
        for y=y_arr
            for z=z_arr
                shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
                patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',RGB(y,z,:));  % draw the blue cube
            end
        end
    end
end

A2d=squeeze(A(x_arr,y_arr(1),z_arr));
index = fix((A2d-cmin)/(cmax-cmin)*m)+1; %A
RGB = ind2rgb(index,cmap);
RGB=rgb_edit_fn(RGB,A2d);
if size(A2d,1)>1 && size(A2d,2)>1
    for x=x_arr
        for y=y_arr(1)
            for z=z_arr
                shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
                patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',RGB(x,z,:));  % draw the blue cube
            end
        end
    end
end

A2d=squeeze(A(x_arr,y_arr(end),z_arr));
index = fix((A2d-cmin)/(cmax-cmin)*m)+1; %A
RGB = ind2rgb(index,cmap);
RGB=rgb_edit_fn(RGB,A2d);
if size(A2d,1)>1 && size(A2d,2)>1
    for x=x_arr
        for y=y_arr(end)
            for z=z_arr
                shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
                patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',RGB(x,z,:));  % draw the blue cube
            end
        end
    end
end

A2d=squeeze(A(x_arr,y_arr,z_arr(end)));
index = fix((A2d-cmin)/(cmax-cmin)*m)+1; %A
RGB = ind2rgb(index,cmap);
RGB=rgb_edit_fn(RGB,A2d);
if size(A2d,1)>1 && size(A2d,2)>1
    for x=x_arr
        for y=y_arr
            for z=z_arr(end)
                shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
                patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',RGB(x,y,:));  % draw the blue cube
            end
        end
    end
end

A2d=squeeze(A(x_arr,y_arr,z_arr(1)));
index = fix((A2d-cmin)/(cmax-cmin)*m)+1; %A
RGB = ind2rgb(index,cmap);
RGB=rgb_edit_fn(RGB,A2d);
if size(A2d,1)>1 && size(A2d,2)>1
    for x=x_arr
        for y=y_arr
            for z=z_arr(1)
                shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
                patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',RGB(x,y,:));  % draw the blue cube
            end
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
xlabel('Y');ylabel('X');zlabel('Z');
camorbit(90,0,'data',[0 1 0])
camorbit(180,0,'data',[0 1 0])

for i=1:length(net.s)
    x=net.s(i).n1i;y=net.s(i).n1j;z=net.s(i).n1k;
    shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
    %patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',net.param.color_arr(mod(i-1,10)+1,:));  % draw the blue cube
    patch('Faces',fac,'Vertices',(vert+shift),'FaceColor','k');  % draw the blue cube
end
for i=1:length(net.s)
    x=net.s(i).n2i;y=net.s(i).n2j;z=net.s(i).n2k;
    shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
    %patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',net.param.color_arr(mod(i-1,10)+1,:));  % draw the blue cube
    patch('Faces',fac,'Vertices',(vert+shift),'FaceColor','k');  % draw the blue cube
end

for i=1:length(net.curr)
    x=net.curr(i).n1i;y=net.curr(i).n1j;z=net.curr(i).n1k;
    shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
    patch('Faces',fac,'Vertices',(vert+shift),'FaceColor','r');  % draw the blue cube
    x=net.curr(i).n2i;y=net.curr(i).n2j;z=net.curr(i).n2k;
    shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
    patch('Faces',fac,'Vertices',(vert+shift),'FaceColor','r');  % draw the blue cube
end

