function plot_single_img_fn(net,opt,clims_in,colormap_in)
% Plot Impedance/Conductivity Image
% Inputs:
% net: EIT Structure
% opt: Plot options: 0 -> impedance, 1-> conductivity, 2->conductivity+electrodes, 3->impedance+electrodes
% clims_in: colormap limits : [min max] or [] for automatic colormap
% colormap_in 0=>jet, 1=>hot
%r_img=net.img_arr_nonan;
r_img=net.img_arr;
r_img_nan=net.img_arr;
if opt==0||opt==3
    A=r_img;
else
    A=1./r_img;
end
if isempty(clims_in)
    clims=[min(min(min(min(A)))) max(max(max(max(A))))];
    cmin=clims(1);cmax=clims(2);
    if cmin==cmax
        cmin=cmin/2;
        cmax=cmax*2;
        clims=[cmin,cmax];
    end
else
    clims=clims_in;
end
figure
set(gcf,'position',[159.40        117.00       1260.00        608.80])
y=round(net.Ni*3/4)+1;
%y=27;
subplot(2,2,1)
if opt==0|opt==1
    img2d=squeeze(r_img(y,:,:))';
else
    img2d=squeeze(r_img_nan(y,:,:))';
end
if opt==2||opt==1
    optmesh=0;
else
    optmesh=1;
end
img2d_extend=[img2d];
img2d_extend(1,:)=mean(img2d_extend(1,img2d_extend(1,:)~=inf&~isnan(img2d_extend(1,:))));
img2d_extend=[img2d(1,:);img2d_extend];
img2d_extend(1,img2d(1,:)~=inf&~isnan(img2d_extend(1,:)))=1e6;

mesh_plot_fn(img2d_extend,optmesh,clims,colormap_in);
% if opt==1
%     Vs_arrow_plot_fn(net);
% end
if opt==2|opt==1
    title(strcat('Conductivity 2D Image at Y=',num2str(y),' mm (S/m)'));
else
    title(strcat('Impedance 2D Image at Y=',num2str(y),' mm (\Omega)'));
end
xlabel('X (mm)');ylabel('Z (mm)');
axis equal;
xlim([0.5 net.Nj+1.5])
ylim([0.5 net.Nk+1.5+1])
ax = gca; ax.FontSize=13;
if net.param.L==2
    for i=1:length(ax.XTickLabel)
        ax.XTickLabel{i}=num2str(str2num(ax.XTickLabel{i})*2);
    end
    for i=1:length(ax.YTickLabel)
        ax.YTickLabel{i}=num2str(str2num(ax.YTickLabel{i})*2-2);
    end
end
set(gca, 'FontName', 'Times New Roman');  

subplot(2,2,2)
x=round(net.curr(1).n1j(1))+round(net.param.i_dia_x/2)-1;
img2d=((squeeze(r_img_nan(:,x,:)))');
if size(net.img_arr,1)==1
    img2d=img2d';
end

img2d_extend=[img2d];
img2d_extend(1,:)=mean(img2d_extend(1,img2d_extend(1,:)~=inf&~isnan(img2d_extend(1,:))));
img2d_extend=[img2d(1,:);img2d_extend];
img2d_extend(1,img2d(1,:)~=inf&~isnan(img2d_extend(1,:)))=1e6;
mesh_plot_fn(img2d_extend,optmesh,clims,colormap_in);
set(gca, 'FontName', 'Times New Roman');  
% if opt==1
%     Vs_arrow_plot_fn(net);
% end
if opt==2|opt==1
    title(strcat('Conductivity 2D Image at X=',num2str(x),' mm (S/m)'));
else
    title(strcat('Impedance 2D Image at X=',num2str(x),' mm (\Omega)'));
end
xlabel('Y (mm)');ylabel('Z (mm)');
axis equal;
ylim([0.5 net.Nk+1.5+1])
xlim([0.5 net.Ni+1.5])
ax = gca; ax.FontSize=12;
set(gca,'Xdir','reverse')
set(gca,'Ydir','reverse')
if net.param.L==2
    for i=1:length(ax.XTickLabel)
        ax.XTickLabel{i}=num2str(str2num(ax.XTickLabel{i})*2);
    end
    for i=1:length(ax.YTickLabel)
        ax.YTickLabel{i}=num2str(str2num(ax.YTickLabel{i})*2-2);
    end
end


subplot(2,2,3)
img2d=squeeze(r_img_nan(:,:,1));
mesh_plot_fn(img2d,optmesh,clims,colormap_in);
if opt==2|opt==1
    title(strcat('Conductivity 2D Image at Z=',num2str(1),' mm (S/m)'));
else
    title(strcat('Impedance 2D Image at Z=',num2str(1),' mm (\Omega)'));
end
xlabel('X (mm)');ylabel('Y (mm)');
axis equal;
xlim([0.5 net.Nj+1.5])
ylim([0.5 net.Ni+1.5])
ax = gca; ax.FontSize=13;
if net.param.L==2
    for i=1:length(ax.XTickLabel)
        ax.XTickLabel{i}=num2str(str2num(ax.XTickLabel{i})*2);
    end
    for i=1:length(ax.YTickLabel)
        ax.YTickLabel{i}=num2str(str2num(ax.YTickLabel{i})*2);
    end
end
set(gca, 'FontName', 'Times New Roman');  

subplot(2,2,4)
plot_core_3d_fn(A,net,1,clims,colormap_in);
set(gca, 'FontName', 'Times New Roman');  
ax = gca; ax.FontSize=10;
if net.param.L==2
    for i=1:length(ax.XTickLabel)
        ax.XTickLabel{i}=num2str(str2num(ax.XTickLabel{i})*2);
    end
    for i=1:length(ax.YTickLabel)
        ax.YTickLabel{i}=num2str(str2num(ax.YTickLabel{i})*2);
    end
    for i=1:length(ax.ZTickLabel)
        ax.ZTickLabel{i}=num2str(str2num(ax.ZTickLabel{i})*2);
    end
end



end


