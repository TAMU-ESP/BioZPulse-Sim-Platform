function net=plot_spice_v_i_3d_fn(net,opt,V_node_min,V_node_max)

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


s_arr=[net.s];
Ve_arr=abs([s_arr.Ve]);
Ve_max=max(Ve_arr,[],'all');
figure
subplot(1,2,1)
V_base=min(abs(net.V_node),[],'all');
A=(abs(net.V_node)-V_base)*1e3;
clims=([V_node_min V_node_max]-V_base)*1e3;
plot_core_3d_fn(A,net,0,clims,0);
title('3D Node Voltage Map (mV)');
xlabel('Y (mm)');ylabel('X (mm)');zlabel('Z (mm)');
ax = gca; ax.FontSize=14;
set(gca, 'FontName', 'Times New Roman');  
set(gcf,'position',[0.2238    0.1406    1.1644    0.6384]*1e3)
%
% if net.param.icalc_en==1
%     subplot(1,2,2)
%     hold on;
%     Z=repmat([1:net.Ni+1]',1,net.Nk+1,net.Nj+1);
%     scale=abs(real(squeeze(max(max(max(net.i_vector(:,:,:,:),[],2),[],3),[],4))));
%     scale_norm=scale/max(scale);
%     for y=1:round(net.Ni)+1%round(net.Ni/2)+1:round(net.Ni)+1
%         ix=squeeze(net.i_vector(y,:,:,1))';
%         iy=squeeze(net.i_vector(y,:,:,2))';
%         iz=squeeze(net.i_vector(y,:,:,3))';
%         U=abs(ix).*sign(real(ix));
%         V=abs(iz).*sign(real(iz));
%         W=abs(iy).*sign(real(iy));
%         quiver3(squeeze(Z(y,:,:)),U,V,W,scale_norm(y),'linewidth',1,'color','r');
%         scale=scale/2;
%     end
%     grid on;
%     %title(strcat('3D Current Distribution at Y=',num2str(round(net.Ni/2)+1)))
%     title(strcat('3D Current Distribution'))
%     xlabel('X (mm)');ylabel('Z (mm)');zlabel('Y (mm)');
%     view(-30,30)
%     ax = gca; ax.FontSize=13;
%     if net.param.L==2
%         for i=1:length(ax.XTickLabel)
%             ax.XTickLabel{i}=num2str(str2num(ax.XTickLabel{i})*2);
%         end
%         for i=1:length(ax.YTickLabel)
%             ax.YTickLabel{i}=num2str(str2num(ax.YTickLabel{i})*2);
%         end
%         for i=1:length(ax.ZTickLabel)
%             ax.ZTickLabel{i}=num2str(str2num(ax.ZTickLabel{i})*2);
%         end
%     end
% end
