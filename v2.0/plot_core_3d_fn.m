function plot_core_3d_fn(A,net,opt,clims,colormap_in)
%opt=0 => no elec, opt=1 => plot electrodes, 3=>log

% if opt==1
%     inf_ind=(A==inf);
%     inf_neg_ind=(A==-inf);
%     A=1./A;
%     A(inf_ind)=inf;
%     A(inf_neg_ind)=-inf;
% end
A=A(:,:,:,1);
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

if colormap_in==0
    colormap jet
else
    colormap hot
end
cmap = colormap;
m = length(cmap);
x_arr=1:Nx;y_arr=1:Ny;z_arr=1:Nz;

A2d=squeeze(A(x_arr(end),y_arr,z_arr));
im=imagesc(A2d,clims);
colorbar
if colormap_in==3
    A2d=log(A2d);
    cmin=log(cmin);
    cmax=log(cmax);
    clims=[cmin cmax];
end
im=imagesc(A2d,clims);
set(im,'Visible','off')
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
if colormap_in==3
    A2d=log(A2d);
end
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
if colormap_in==3
    A2d=log(A2d);
end
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
if colormap_in==3
    A2d=log(A2d);
end
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
if colormap_in==3
    A2d=log(A2d);
end
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
if colormap_in==3
    A2d=log(A2d);
end
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
xlabel('Y (mm)');ylabel('X (mm)');zlabel('Z (mm)');
%view(45,30);
camorbit(270,0,'data',[0 0 1])
camorbit(180,0,'data',[1 0 0])
set(gca,'Ydir','normal')
camorbit(-30,0,'data',[0 1 0])
camorbit(10,0,'data',[1 0 0])
%set(gca,'DataAspectRatio',[Nx Ny Nz])
%material shiny;
%alpha('color');
%alphamap('rampdown');

% camorbit(90,0,'data',[0 1 0])
% camorbit(180,0,'data',[0 1 0])
% camorbit(180,0,'data',[0 0 1])
% %[a,b,c]=get(gca,'DataAspectRatio');
%set(gca,'DataAspectRatio',[Nz*a Ny Nx])
set(gca,'DataAspectRatio',[1 1 1])
if opt==1
    if net.param.se_en==1
        for i=1:length(net.s)
            n1imax=max(net.s(i).n1i);n1jmax=max(net.s(i).n1j);
            for j=1:length(net.s(i).n1i)
                if net.s(i).n1i(j)<n1imax && net.s(i).n1j(j)<n1jmax
                    x=net.s(i).n1i(j);y=net.s(i).n1j(j);z=net.s(i).n1k(j)-1;
                    shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
                    %patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',net.param.color_arr(mod(i-1,10)+1,:));  % draw the blue cube
                    patch('Faces',fac,'Vertices',(vert+shift),'FaceColor','b');  % draw the blue cube
                end
            end
        end
    else
        for i=1:length(net.s)
            n1imax=max(net.s(i).n1i);n1jmax=max(net.s(i).n1j);
            for j=1:length(net.s(i).n1i)
                if net.s(i).n1i(j)<n1imax && net.s(i).n1j(j)<n1jmax
                    x=net.s(i).n1i(j);y=net.s(i).n1j(j);z=net.s(i).n1k(j);
                    shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
                    %patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',net.param.color_arr(mod(i-1,10)+1,:));  % draw the blue cube
                    patch('Faces',fac,'Vertices',(vert+shift),'FaceColor','b');  % draw the blue cube
                end
            end
        end
    end
    
    if net.param.se_en==1
        for i=1:length(net.s)
            n2imax=max(net.s(i).n2i);n2jmax=max(net.s(i).n2j);
            for j=1:length(net.s(i).n2i)
                if net.s(i).n2i(j)<n2imax && net.s(i).n2j(j)<n2jmax
                    x=net.s(i).n2i(j);y=net.s(i).n2j(j);z=net.s(i).n2k(j)-1;
                    shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
                    %patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',net.param.color_arr(mod(i-1,10)+1,:));  % draw the blue cube
                    patch('Faces',fac,'Vertices',(vert+shift),'FaceColor','b');  % draw the blue cube
                end
            end
        end
    else
        for i=1:length(net.s)
            n2imax=max(net.s(i).n2i);n2jmax=max(net.s(i).n2j);
            for j=1:length(net.s(i).n2i)
                if net.s(i).n2i(j)<n2imax && net.s(i).n2j(j)<n2jmax
                    x=net.s(i).n2i(j);y=net.s(i).n2j(j);z=net.s(i).n2k(j);
                    shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
                    %patch('Faces',fac,'Vertices',(vert+shift),'FaceColor',net.param.color_arr(mod(i-1,10)+1,:));  % draw the blue cube
                    patch('Faces',fac,'Vertices',(vert+shift),'FaceColor','b');  % draw the blue cube
                end
            end
        end
    end
    
    if net.param.se_en==1
        for i=1:length(net.curr)
            n1imax=max(net.curr(i).n1i);n1jmax=max(net.curr(i).n1j);
            for j=1:length(net.curr(i).n1i)
                if net.curr(i).n1i(j)<n1imax && net.curr(i).n1j(j)<n1jmax
                    x=net.curr(i).n1i(j);y=net.curr(i).n1j(j);z=net.curr(i).n1k(j)-1;
                    shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
                    patch('Faces',fac,'Vertices',(vert+shift),'FaceColor','r');  % draw the blue cube
                    x=net.curr(i).n2i(j);y=net.curr(i).n2j(j);z=net.curr(i).n2k(j)-1;
                    shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
                    patch('Faces',fac,'Vertices',(vert+shift),'FaceColor','r');  % draw the blue cube
                end
            end
        end
    else
        for i=1:length(net.curr)
            n1imax=max(net.curr(i).n1i);n1jmax=max(net.curr(i).n1j);
            for j=1:length(net.curr(i).n1i)
                if net.curr(i).n1i(j)<n1imax && net.curr(i).n1j(j)<n1jmax
                x=net.curr(i).n1i(j);y=net.curr(i).n1j(j);z=net.curr(i).n1k(j);
                shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
                patch('Faces',fac,'Vertices',(vert+shift),'FaceColor','r');  % draw the blue cube
                x=net.curr(i).n2i(j);y=net.curr(i).n2j(j);z=net.curr(i).n2k(j);
                shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
                patch('Faces',fac,'Vertices',(vert+shift),'FaceColor','r');  % draw the blue cube
                end
            end
        end
    end
%     if net.param.re_var_max(1)>0
%         for i=1:length(net.r_var)
%             for j=1:length(net.r_var(i).n1ibase)
%                 x=net.r_var(i).n1ibase(j);y=net.r_var(i).n1jbase(j);z=net.r_var(i).n1kbase(j);
%                 shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
%                 patch('Faces',fac,'Vertices',(vert+shift),'FaceColor','r');  % draw the blue cube
%             end
%             x=net.r_var(i).n1ibase(end)+1;y=net.r_var(i).n1jbase(end);z=net.r_var(i).n1kbase(end);
%             shift=[ones(8,1)*x,ones(8,1)*y,ones(8,1)*z];
%             patch('Faces',fac,'Vertices',(vert+shift),'FaceColor','r');  % draw the blue cube
%         end
%     end
end
