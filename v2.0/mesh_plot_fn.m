function mesh_plot_fn(A,opt,clims,colormap_in)
% A = rand(50); %random data
% for i = 1:10 %turn some random points into NaNs
%    A(randi(50),randi(50),:) = NaN;
% end
% opt=0 -> No color, opt=1 -> With color, opt=2 -> invert, opt=3 ->log

%A=flipud(A);
%set(gca,'ColorScale','log')
if opt==2
    inf_ind=(A==inf);
    inf_neg_ind=(A==-inf);
    A=1./A;
    A(inf_ind)=inf;
    A(inf_neg_ind)=-inf;
end


if opt>0
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
    imagesc(A,clims)
    hold on;
    axP = get(gca,'Position');
    colorbar;
    set(gca, 'Position', axP);
    if colormap_in==0
        colormap jet
    else
        colormap hot
    end
    cmap = colormap;
    
    m = length(cmap);
    if colormap_in==3
        A=log(A);
        cmin=log(cmin);
        cmax=log(cmax);
    end
    index = fix((A-cmin)/(cmax-cmin)*m)+1; %A
    RGB = ind2rgb(index,cmap);
else
    RGB=ones(size(A,1),size(A,2),3);
end
if colormap_in==3
    set(gca,'ColorScale','log')
    % preallocate Ticks and TickLabels
    num_of_ticks = 5;
    Ticks      = [];
    TickLabels = [];
    % distribute Ticks and TickLabels
    for n = 1:1:num_of_ticks
        
        Ticks(n)      = 2^floor(log2((exp(cmax)/1000/(2^(num_of_ticks-1)))))*1000*(2^(n-1));
        TickLabels{n} = num2str(Ticks(n));
    end
    %Ticks=[100,500,1000,2000,4000];
    %TickLabels={'100';'500';'1000';'2000';'4000'};
    colorbar('Ticks',Ticks,'TickLabels',TickLabels)
end
R = RGB(:,:,1); % turn your data into a "pseudo-gray" rgb image.
G = RGB(:,:,2);
B = RGB(:,:,3);
if opt~=0
    R(isnan(A)) = 1; %turn NaNs(I) to Red
    G(isnan(A)) = 0;
    B(isnan(A)) = 0;
    
    R(A==inf) = 0; %turn infs(V) to blue
    G(A==inf) = 0;
    B(A==inf) = 1;
    
    R(A==-inf) = 1; %turn -infs(A) to green
    G(A==-inf) = 0;
    B(A==-inf) = 0;
else
    gray_scale=0.75;
    R(isnan(A)) = gray_scale; %turn NaNs(I) to Red
    G(isnan(A)) = gray_scale;
    B(isnan(A)) = gray_scale;
    
    R(A==inf) = gray_scale; %turn infs(V) to blue
    G(A==inf) = gray_scale;
    B(A==inf) = gray_scale;
    
    R(A==-inf) = gray_scale; %turn -infs(A) to green
    G(A==-inf) = gray_scale;
    B(A==-inf) = gray_scale;
end


RGB2(:,:,1) = R; %combine color slices into A
RGB2(:,:,2) = G;
RGB2(:,:,3) = B;
%imagesc(A)
im=imagesc(RGB2,'CDataMapping','direct');
%set(gca,'ColorScale','log')
im.CDataMapping='direct';
impixelinfo;

dcm_obj = datacursormode();
set(dcm_obj,'UpdateFcn',@(hObject, event_obj) myupdatefcn(hObject, event_obj, A) );

hold on;
for i = 1:max(size(A))
    plot([i-.5,i-.5],[.5,max(size(A))+.5],'k-');
end
for i = 1:max(size(A))
    plot([.5,max(size(A))+.5],[i-.5,i-.5],'k-');
end
set(gca,'YDir','normal')
axis equal;
xlim([0.5 size(A,2)+0.5])
ylim([0.5 size(A,1)+0.5])


