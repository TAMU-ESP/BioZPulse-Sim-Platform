function net=plot_spice_pulse_v_i_2d_fn(net,opt,y,V_node_min,V_node_max)

%   Bassem Ibrahim
%   2019.08.22
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation
%   Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling,
%   IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan.
%
%   function net=spice_netlist_3d_plot_fn(net)
%   Plot Output Results

V_base=min(abs(net.V_node(:,y,:,1)),[],'all');
figure
set(gcf,'position',[0.2238    0.1406    1.1644    0.6384]*1e3)
subplot(2,2,1)
mesh_plot_fn(squeeze(abs(net.V_node(:,y,:,1))-V_base)'*1e3,1,([V_node_min V_node_max]-V_base)*1e3,0);
title(strcat('2D Node and Electrode Voltage(Ve) (mV) at Y=',num2str(y)));
%plot_velec_fn(net,opt);
xlabel('Y (mm)');ylabel('Z (mm)')
set(gca, 'FontName', 'Times New Roman');  
ax = gca; ax.FontSize=14;
set(gca,'xdir','reverse','ydir','reverse')
if net.param.L==2
    for i=1:length(ax.XTickLabel)
        ax.XTickLabel{i}=num2str(str2num(ax.XTickLabel{i})*2);
    end
    for i=1:length(ax.YTickLabel)
        ax.YTickLabel{i}=num2str(str2num(ax.YTickLabel{i})*2);
    end
end


if net.param.icalc_en==1
    subplot(2,2,2)
    img2d=squeeze(net.img_arr(:,y,:))';
    img2d_extend=[img2d];
    img2d_extend(1,:)=mean(img2d_extend(1,img2d_extend(1,:)~=inf&~isnan(img2d_extend(1,:))));
    img2d_extend=[img2d(1,:);img2d_extend];
    img2d_extend(1,img2d(1,:)~=inf&~isnan(img2d_extend(1,:)))=1e6;

    mesh_plot_fn(img2d_extend,0,[],1);
     set(gca, 'FontName', 'Times New Roman');  
    %Vs_arrow_plot_fn(net)
    xlabel('Y (mm)');ylabel('Z (mm)')
    title(strcat('2D Current Distribution at Y=',num2str(round(net.Ni/2)+1)))
    iy=squeeze(net.i_vector(:,y,:,2,1))';
    iz=squeeze(net.i_vector(:,y,:,3,1))';
    quiver([0:net.Ni]+1,[0:net.Nk]+2,abs(iy).*sign(real(iy))*-1,abs(iz).*sign(real(iz))*-1,'linewidth',1,'color','r');
    ax = gca; ax.FontSize=14;
    set(gca,'xdir','reverse','ydir','reverse')
    if net.param.L==2
        for i=1:length(ax.XTickLabel)
            ax.XTickLabel{i}=num2str(str2num(ax.XTickLabel{i})*2);
        end
        for i=1:length(ax.YTickLabel)
            ax.YTickLabel{i}=num2str(str2num(ax.YTickLabel{i})*2-2);
        end
    end
   
end

subplot(2,2,3)
mesh_plot_fn(squeeze(abs(net.V_node(:,:,1,1))-V_base)*1e3,1,([V_node_min V_node_max]-V_base)*1e3,0);
title(strcat('2D Node and Electrode Voltage(Ve) (mV) at Z=0'));
%plot_velec_fn(net,opt);
xlabel('X (mm)');ylabel('Y (mm)')
ax = gca; ax.FontSize=14;
set(gca, 'FontName', 'Times New Roman');  
%set(gca,'xdir','reverse','ydir','reverse')
if net.param.L==2
    for i=1:length(ax.XTickLabel)
        ax.XTickLabel{i}=num2str(str2num(ax.XTickLabel{i})*2);
    end
    for i=1:length(ax.YTickLabel)
        ax.YTickLabel{i}=num2str(str2num(ax.YTickLabel{i})*2);
    end
end


if net.param.icalc_en==1
    subplot(2,2,4)
    set(gca,'xdir','reverse','ydir','reverse')
    z=4;
    mesh_plot_fn(squeeze(net.img_arr(:,:,z)),0,[],1);
     set(gca, 'FontName', 'Times New Roman');  
    %Vs_arrow_plot_fn(net)
    xlabel('X (mm)');ylabel('Y (mm)')
    title(strcat('2D Current Distribution at Z=',num2str(z)))
    ix=squeeze(net.i_vector(:,:,z,2,1));
    iy=squeeze(net.i_vector(:,:,z,1,1));
    quiver([0:net.Nj]+2,[0:net.Ni]+2,abs(iy).*sign(real(iy)),abs(ix).*sign(real(ix)),'linewidth',1,'color','r');
    ax = gca; ax.FontSize=14;    
    if net.param.L==2
        for i=1:length(ax.XTickLabel)
            ax.XTickLabel{i}=num2str(str2num(ax.XTickLabel{i})*2);
        end
        for i=1:length(ax.YTickLabel)
            ax.YTickLabel{i}=num2str(str2num(ax.YTickLabel{i})*2-2);
        end
    end
   
end
