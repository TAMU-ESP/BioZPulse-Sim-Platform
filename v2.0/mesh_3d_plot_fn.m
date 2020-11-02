function mesh_3d_plot_fn(A,opt,clims)
% A = rand(50); %random data
% for i = 1:10 %turn some random points into NaNs
%    A(randi(50),randi(50),:) = NaN;
% end
% opt=0 -> No color, opt=1 -> With color

%A=flipud(A);

if opt==1
    %[Aim]=imshow(A);
    %colormap(gray);
    if isempty(clims)
        cmin = min(A(:));
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
    %imagesc(A,clims)
    hold on;
    axP = get(gca,'Position');
    colorbar;
    set(gca, 'Position', axP);
    
    cmap = colormap;
    
    m = length(cmap);    
    index = fix((A-cmin)/(cmax-cmin)*m)+1; %A
    RGB = ind2rgb(index,cmap);
else
    RGB=ones(size(A,1),size(A,2),3);
end

R = RGB(:,:,1); % turn your data into a "pseudo-gray" rgb image.
G = RGB(:,:,2);
B = RGB(:,:,3);
R(isnan(A)) = 1; %turn NaNs(I) to Red
G(isnan(A)) = 0;
B(isnan(A)) = 0;

R(A==inf) = 0; %turn infs(V) to Black
G(A==inf) = 0;
B(A==inf) = 0;

R(A==-inf) = 0; %turn infs(V) to Black
G(A==-inf) = 1;
B(A==-inf) = 0;

RGB2(:,:,1) = R; %combine color slices into A
RGB2(:,:,2) = G;
RGB2(:,:,3) = B;
%imagesc(A)
im=imagesc(RGB2,'CDataMapping','direct');
im.CDataMapping='direct';
impixelinfo;

dcm_obj = datacursormode();
set(dcm_obj,'UpdateFcn',@(hObject, event_obj) myupdatefcn(hObject, event_obj, A) );

hold on;
for i = 1:size(A,1)
    plot([i-.5,i-.5],[.5,size(A,1)+.5],'k-');
end
for i = 1:size(A,2)
    plot([.5,size(A,2)+.5],[i-.5,i-.5],'k-');
end
set(gca,'YDir','normal')
axis equal;
xlim([0.5 size(A,2)+0.5])
ylim([0.5 size(A,1)+0.5])


