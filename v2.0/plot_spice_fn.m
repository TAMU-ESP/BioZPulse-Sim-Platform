function net=plot_spice_fn(net,opt,V_node_min,V_node_max,Z_img_max,Z_img_min)

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

L=net.L;

spacing_arr=net.spacing_arr;
if sum(net.param.re_var_max>0)~=0
    r_var_t=net.r_var.r_var_t;
end

s_arr=[net.s];
Vs_arr=abs([s_arr.Vs]);
Vs_max=max(Vs_arr,[],'all');

fi=0;
% 1 All V Heatmaps
fi=fi+1;
if any(fi==net.param.fi_req)
    %     V_sns_v_max=max(max(max((squeeze(net.V_sns_v_NaN(:,:,1,:))))));
    %     V_sns_h_max=max(max(max((squeeze(net.V_sns_h_NaN(:,:,1,:))))));
    %     V_sns_v_min=min(min(min((squeeze(net.V_sns_v_NaN(:,:,1,:))))));
    %     V_sns_h_min=min(min(min((squeeze(net.V_sns_h_NaN(:,:,1,:))))));
    %     V_sns_max=max(V_sns_v_max);
    %     V_sns_min=min(V_sns_v_min);
    %     V_node_max=max(max(max((squeeze(net.V_node_NaN(:,:,1,:))))));
    %     V_node_min=min(min(min((squeeze(net.V_node_NaN(:,:,1,:))))));
    fig_h(fi)=figure;
    %     if net.Nk>0
    %         subplot(3,2,1)
    %         %heatmap([0:net.Nj],[0:net.Ni],squeeze(net.node_arr(:,:,1)),'Colormap',white,'MissingDataLabel','Elec.');
    %         mesh_plot_fn(squeeze(net.img_arr(:,:,1)));
    %         xlabel('X');ylabel('Y')
    %         title('Electrodes Location (Black square), Z=0')
    %         subplot(3,2,3)
    %         heatmap([0:net.Nj],[0:net.Nk],squeeze(net.node_arr(1,:,:))','Colormap',white,'MissingDataLabel','Artery');
    %         xlabel('X');ylabel('Z')
    %         title('Artery Location (Black square), Z=0')
    %
    %         subplot(3,2,5)
    %         heatmap([0:net.Nj],[0:net.Ni],net.V_node(:,:,1,1),'Colormap',parula,'ColorLimits',[V_node_min V_node_max]);
    %         %mesh_plot_fn(normalize(net.V_node(:,:,1,1)));
    %         xlabel('X');ylabel('Y')
    %         title('Node V (mV)')
    %
    %         k_index=1;
    %         subplot(3,2,2)
    %         heatmap([0:net.Nj],[0:net.Ni-1],net.V_sns_v(:,:,k_index,1),'Colormap',parula,'ColorLimits',[V_sns_min V_sns_max]);
    %         %mesh_plot_fn(net.V_sns_v(:,:,k_index,1));
    %         xlabel('X');ylabel('Y')
    %         title(strcat('V_{DC}(mV), Z=',num2str(k_index-1),',spacing=',num2str(net.spacing_single)));
    %
    %         k_index=2;
    %         subplot(3,2,4)
    %         heatmap([0:net.Nj],[0:net.Ni-1],net.V_sns_v_NaN(:,:,k_index,1),'Colormap',parula,'ColorLimits',[V_sns_min V_sns_max]);
    %         xlabel('X');ylabel('Y')
    %         title(strcat('V_{DC}(mV), Z=',num2str(k_index-1)))
    %         k_index=3;
    %         subplot(3,2,6)
    %         heatmap([0:net.Nj],[0:net.Ni-1],net.V_sns_v_NaN(:,:,k_index,1),'Colormap',parula,'ColorLimits',[V_sns_min V_sns_max]);
    %         xlabel('X');ylabel('Y')
    %         title(strcat('V_{DC}(mV), Z=',num2str(k_index-1)))
    %
    %     else
    
    if opt==1
        subplot(1,2,1)
    else
        subplot(2,2,1)
    end
    
    if Vs_max<10
        mesh_plot_fn(squeeze(abs(net.V_node(round(end/2),:,:)))'*1e3,1,[V_node_min V_node_max]*1e3);
        title('Node and Electrode Voltage(Ve) (uV)');
    else
        mesh_plot_fn(squeeze(abs(net.V_node(round(end/2),:,:)))',1,[V_node_min V_node_max]);
        title('Node and Electrode Voltage(Ve) (mV)');
    end
    plot_velec_fn(net,opt);
    xlabel('X');ylabel('Z')
    ax = gca; ax.FontSize=13;
    
    %         subplot(2,2,3)
    %         mesh_plot_fn(squeeze(net.img_arr(:,:,1)),0);
    %         %Vs_arrow_plot_fn(net)
    %         xlabel('X');ylabel('Y')
    %         title({'Electrode Voltage Difference V_e(mV)',''})
    %         plot_velec_fn(net);
    if opt==1
        subplot(1,2,2)
    else
        subplot(2,2,2)
    end
    
    
        
        
    
    hold on;
    if net.Ni>1
        %mesh_plot_fn(squeeze(net.img_arr(round(y),:,:)),0);
        %Vs_arrow_plot_fn(net)
        xlabel('X');ylabel('Z')
        title('Current Distribution')
         Z=repmat([1:net.Ni+1]',1,net.Nj+1,net.Nk+1);
        for y=7%round(net.Ni/2)+1%1:net.Ni+1
            ix=squeeze(net.i_vector(y,:,:,1))';
            iy=squeeze(net.i_vector(y,:,:,2))';
            iz=squeeze(net.i_vector(y,:,:,3))';           
            U=abs(ix).*sign(real(ix));
            V=abs(iy).*sign(real(iy));
            W=abs(iz).*sign(real(iz));
            quiver3(squeeze(Z(y,:,:)),U,V,W,'linewidth',1,'color','r');
        end
        grid on;
        xlabel('Y');ylabel('X');zlabel('Z');
        view(-35,45)
        ax = gca; ax.FontSize=13;
        
%         mesh_plot_fn(squeeze(net.img_arr(round(net.Ni/2)+1,:,:))',0);
%         %Vs_arrow_plot_fn(net)
%         xlabel('X');ylabel('Z')
%         title('Current Distribution')
%         ix=squeeze(net.i_vector(round(net.Ni/2)+1,:,:,1))';
%         iz=squeeze(net.i_vector(round(net.Ni/2)+1,:,:,3))';
%         quiver([0:net.Nj]+1,[0:net.Nk]+1,abs(ix).*sign(real(ix)),abs(iz).*sign(real(iz)),'linewidth',1,'color','r');
%         ax = gca; ax.FontSize=13;
        
    else
        mesh_plot_fn(squeeze(net.img_arr(1,:,:))',0);
        %Vs_arrow_plot_fn(net)
        xlabel('X');ylabel('Z')
        title('Current Distribution')
        ix=squeeze(net.i_vector(1,:,:,1))';
        iz=squeeze(net.i_vector(1,:,:,3))';
        quiver([0:net.Nj]+1,[0:net.Nk]+1,abs(ix).*sign(real(ix)),abs(iz).*sign(real(iz)),'linewidth',1,'color','r');
        ax = gca; ax.FontSize=13;
    end
    
    
    if opt==2
        subplot(2,2,3)
        if net.Ni>1
            mesh_plot_fn(abs(squeeze(net.Y_img(round(end/2)+1,:,:)))',1,[Z_img_min,Z_img_max]);
        else
            mesh_plot_fn(abs(squeeze(net.Y_img(1,:,:)))',1,[Z_img_min,Z_img_max]);
        end
        xlabel('X');ylabel('Z')
        title('Conductivity Image')
        ax = gca; ax.FontSize=13;
        
        %             subplot(2,2,4)
        %             grid on;hold on;
        %             plot([1:length(net.plot.freq_arr)],net.plot.freq_arr/1e3,'linewidth',2)
        %             plot(net.plot.freq_ind,net.plot.freq_arr(net.plot.freq_ind)/1e3,'or','linewidth',2)
        %             xlim([1 length(net.plot.freq_arr)])
        %             xlabel('Iteration');ylabel('Frequency(kHz)')
        %             title('Frequency Step')
        %             ax = gca; ax.FontSize=13;
    end
    
    set(gcf,'position',[0.0010    0.0410    1.2800    0.6073]*1e3);
    
    if net.param.write_figures==1 && net.param.debug==0
        %             if ~exist(net.param.plot_path,'dir')
        %                 mkdir(net.param.plot_path);
        %             end
        file_str=strcat(net.param.plot_path);
        print(strcat(file_str,'eit_',num2str(net.param.index),'.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
        savefig(strcat(file_str,'eit_',num2str(net.param.index),'.fig')); %<-Save as PNG with 300 DP
    end
    %  end
    
end
