function net=spice_netlist_3d_plot_fn(net)

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
S=net.param.S;
spacing_arr=net.spacing_arr;
if sum(net.param.r_var_max>0)~=0
    r_var_t=net.r_var.r_var_t;
end


fi=0;
% 1 All V Heatmaps
fi=fi+1;
if any(fi==net.param.fi_req)
    V_sns_v_max=max(max(max((squeeze(net.V_sns_v_NaN(:,:,1,:))))));
    V_sns_h_max=max(max(max((squeeze(net.V_sns_h_NaN(:,:,1,:))))));
    V_sns_v_min=min(min(min((squeeze(net.V_sns_v_NaN(:,:,1,:))))));
    V_sns_h_min=min(min(min((squeeze(net.V_sns_h_NaN(:,:,1,:))))));
    V_sns_max=max(V_sns_v_max);
    V_sns_min=min(V_sns_v_min);
    V_node_max=max(max(max((squeeze(net.V_node_NaN(:,:,1,:))))));
    V_node_min=min(min(min((squeeze(net.V_node_NaN(:,:,1,:))))));
    fig_h(fi)=figure;
    if net.Nk>0
        subplot(3,2,1)
        %heatmap([0:net.Nj],[0:net.Ni],squeeze(net.node_arr(:,:,1)),'Colormap',white,'MissingDataLabel','Elec.');
        mesh_plot_fn(squeeze(net.img_arr(:,:,1)));
        xlabel('X');ylabel('Y')
        title('Electrodes Location (Black square), Z=0')
        subplot(3,2,3)
        heatmap([0:net.Nj],[0:net.Nk],squeeze(net.node_arr(1,:,:))','Colormap',white,'MissingDataLabel','Artery');
        xlabel('X');ylabel('Z')
        title('Artery Location (Black square), Z=0')
        
        subplot(3,2,5)
        heatmap([0:net.Nj],[0:net.Ni],net.V_node(:,:,1,1),'Colormap',parula,'ColorLimits',[V_node_min V_node_max]);
        %mesh_plot_fn(normalize(net.V_node(:,:,1,1)));
        xlabel('X');ylabel('Y')
        title('Node V (mV)')
        
        k_index=1;
        subplot(3,2,2)
        heatmap([0:net.Nj],[0:net.Ni-1],net.V_sns_v(:,:,k_index,1),'Colormap',parula,'ColorLimits',[V_sns_min V_sns_max]);
        %mesh_plot_fn(net.V_sns_v(:,:,k_index,1));
        xlabel('X');ylabel('Y')
        title(strcat('V_{DC}(mV), Z=',num2str(k_index-1),',spacing=',num2str(net.spacing_single)));
        
        k_index=2;
        subplot(3,2,4)
        heatmap([0:net.Nj],[0:net.Ni-1],net.V_sns_v_NaN(:,:,k_index,1),'Colormap',parula,'ColorLimits',[V_sns_min V_sns_max]);
        xlabel('X');ylabel('Y')
        title(strcat('V_{DC}(mV), Z=',num2str(k_index-1)))
        k_index=3;
        subplot(3,2,6)
        heatmap([0:net.Nj],[0:net.Ni-1],net.V_sns_v_NaN(:,:,k_index,1),'Colormap',parula,'ColorLimits',[V_sns_min V_sns_max]);
        xlabel('X');ylabel('Y')
        title(strcat('V_{DC}(mV), Z=',num2str(k_index-1)))
        
    else
        subplot(1,2,1)
        %heatmap([0:net.Nj],[0:net.Ni],squeeze(net.node_arr(:,:,1)),'Colormap',white,'MissingDataLabel','Elec.');
        mesh_plot_fn(squeeze(net.img_arr(:,:,1)));
        xlabel('X');ylabel('Y')
        title('Impedance Image Top View, I(Red), V(Black), Z=0')
        
        subplot(1,2,2)
        h=heatmap([0:net.Nj]+1,[0:net.Ni]+1,net.V_node(:,:,1,1),'Colormap',parula,'ColorLimits',[V_node_min V_node_max]);
        idx=[0:net.Ni]+1;        
        h.XDisplayLabels((mod(idx,5))~=0) = {''};
        idy=[0:net.Nj]+1;        
        h.YDisplayLabels((mod(idy,5))~=0) = {''};
        %mesh_plot_fn((net.V_node(:,:,1,1)));
        xlabel('X');ylabel('Y')
        title('Node V (mV)');
        if net.param.write_figures==1
            file_str=strcat(net.param.path,net.param.write_netlist);
            print(strcat(file_str,'_eit_',num2str(net.param.index),'.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
            savefig(strcat(file_str,'_eit_',num2str(net.param.index),'.fig')); %<-Save as PNG with 300 DP
        end

                
        %         k_index=1;
        %         subplot(3,2,2)
        %         heatmap([0:net.Nj],[0:net.Ni-1],net.V_sns_v(:,:,k_index,1),'Colormap',parula,'ColorLimits',[V_sns_min V_sns_max]);
        %         %mesh_plot_fn(net.V_sns_v(:,:,k_index,1));
        %         xlabel('X');ylabel('Y')
        %         title(strcat('V_{DC}(mV), Z=',num2str(k_index-1),',spacing=',num2str(net.spacing_single)));
        %
    end
    
    %set(gcf, 'Position',[ 93.8000   20.2000  746.4000  742.4000])
end

% Location change (Z=0)
V_sns_v_loc=squeeze(net.V_sns_v(net.curr.n1i(1)+1:net.curr.n2i(1)-2,net.curr.n1j(1),1,1));
for p=1:length(spacing_arr)
    V_sns_v_spacing_p=net.V_sns_v_space{p}(net.curr.n1i(1)+1+net.curr.dia_y:net.curr.n2i(1)-1-spacing_arr(p),net.curr.n1j(1)-1,1,1);
    V_sns_v_spacing(p)=V_sns_v_spacing_p(round(length(V_sns_v_spacing_p)/2));
end
% 2
fi=fi+1;
if any(fi==net.param.fi_req)
    fig_h(fi)=figure;
    subplot(2,3,1)
    hold on;grid on;
    plot([0:net.Ni-1],squeeze(net.V_sns_v_NaN(:,net.curr.n1j(1)-1,1,1)),'linewidth',2)
    plot([0:net.Ni-1],squeeze(net.V_sns_v_NaN(:,round(net.Nj/4),1,1)),'linewidth',2)
    plot([0:net.Ni-1],squeeze(net.V_sns_v_NaN(:,1,1,1)),'linewidth',2)
    xlabel('Y')
    ylabel('V_{DC}(mV)')
    ax = gca; ax.FontSize=13;
    subplot(2,3,4);
    hold on;grid on;
    plot([0:net.Ni-1],squeeze(net.V_sns_v_diff(:,net.curr.n1j(1)-1,1,1)),'linewidth',2)
    xlabel('Y')
    ylabel('\DeltaV_{PP}(mV)')
    ax = gca; ax.FontSize=13;
    subplot(2,3,2)
    hold on;grid on;
    plot([1:net.Nj+1]-1,squeeze([net.V_sns_v(min(round(mean([max(net.curr.n1i)+1])),net.Ni-1),:,1,1)]),'linewidth',2)
    if ~isempty(net.param.s_n1i)
        for i=1:length(net.s)
            plot([1:net.Nj+1]-1,squeeze([net.V_sns_v_NaN(min(round(mean([mean(net.s(i).n1i),mean(net.s(i).n2i)])),net.Ni-1),:,1,1)]),'linewidth',2)
        end
        str2=strcat('Y=',num2str(-1+round(mean([mean(net.s(i).n1i),mean(net.s(i).n2i)]))));
    else
        plot([1:net.Nj+1]-1,squeeze([net.V_sns_v_NaN(min(round(mean([mean(net.curr.n1i),mean(net.curr.n2i)])),net.Ni-1),:,1,1)]),'linewidth',2)
        str2=strcat('Y=',num2str(-1+min(round(mean([mean(net.curr.n1i),mean(net.curr.n2i)])),net.Ni-1)));
    end
    legend(strcat('Y=',num2str(mean([max(net.curr.n1i)]))),str2);
    xlabel('X')
    ylabel('V_{DC}(mV)')
    ax = gca; ax.FontSize=13;
    subplot(2,3,5);
    hold on;grid on;
    if ~isempty(net.param.s_n1i)
        for i=1:length(net.s)
            plot([1:net.Nj+1]-1,squeeze(net.V_sns_v_diff(min(round(mean([mean(net.s(i).n1i),mean(net.s(i).n2i)])),net.Ni-1),:,1)),'linewidth',2)
        end
    else
        plot([1:net.Nj+1]-1,squeeze(net.V_sns_v_diff(min(round(mean([net.curr.n1i(1),net.curr.n2i(1)]),net.Ni-1)),:,1)),'linewidth',2)
    end
    xlabel('X')
    ylabel('\DeltaV_{PP}(mV)')
    ax = gca; ax.FontSize=13;
    subplot(2,3,3)
    hold on;grid on;
    plot([1:size(net.V_sns_v,3)]-1,squeeze([net.V_sns_v(min(max(net.curr.n1i)+1,net.Ni-1),min(net.curr.n1j(1),net.Nj-1),:,1)]),'linewidth',2)
    plot([1:size(net.V_sns_v,3)]-1,squeeze([net.V_sns_v_NaN(round(mean([net.curr.n1i(1),net.curr.n2i(1)])),min(net.curr.n1j(1),net.Nj-1),:,1)]),'linewidth',2)
    legend(strcat('Y=',num2str(net.curr.n1i)),strcat('Y=',num2str(-1+round(mean([net.curr.n1i(1),net.curr.n2i(1)])))));
    xlabel('Z')
    ylabel('V_{DC}(mV)')
    ax = gca; ax.FontSize=13;
    
    subplot(2,3,6)
    hold on;grid on;
    plot([1:size(net.V_sns_v_diff,3)]-1,squeeze([net.V_sns_v_diff(min(max(net.curr.n1i)+1,net.Ni-1),min(net.curr.n1j(1),net.Nj-1),:)]),'linewidth',2)
    plot([1:size(net.V_sns_v_diff,3)]-1,squeeze([net.V_sns_v_diff(round(mean([net.curr.n1i(1),net.curr.n2i(1)])),min(net.curr.n1j(1),net.Nj-1),:)]),'linewidth',2)
    legend(strcat('Y=',num2str(net.curr.n1i)),strcat('Y=',num2str(-1+round(mean([net.curr.n1i(1),net.curr.n2i(1)])))));
    xlabel('Z')
    ylabel('V_{DC}(mV)')
    ax = gca; ax.FontSize=13;
    set(gcf, 'Position',[ 0.0794    0.0010    1.4576    0.7888]*1e3)
end

% 3 Delta V SNS Z change Heatmaps
if net.param.r_var_max>0
    fi=fi+1;
    if any(fi==net.param.fi_req)
        V_sns_v_diff_abs=abs(squeeze(net.V_sns_v_diff(:,:,:)));
        V_sns_h_diff_abs=abs(squeeze(net.V_sns_h_diff(:,:,:)));
        V_node_diff_abs=abs(squeeze(net.V_node_diff_NaN(:,:,:)));
        V_sns_diff_max=max(max(V_sns_v_diff_abs(:)),max(V_sns_h_diff_abs(:)));
        V_node_diff_max=max(net.V_node_diff_NaN(:));
        V_node_diff_min=min(net.V_node_diff_NaN(:));
        V_sns_v_diff_max=max(net.V_sns_v_diff(:));
        V_sns_v_diff_min=min(net.V_sns_v_diff(:));
        V_sns_h_diff_max=max(net.V_sns_h_diff(:));
        V_sns_h_diff_min=min(net.V_sns_h_diff(:));
        fig_h(fi)=figure;
        k_index=1;
        subplot(3,3,1)
        heatmap([0:net.Nj],[0:net.Ni],net.V_node_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[V_node_diff_min V_node_diff_max]);
        title(strcat('Node Volatge Delta V(mV), Z=',num2str(k_index-1)))
        subplot(3,3,2)
        heatmap([0:net.Nj],[0:net.Ni-1],net.V_sns_v_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[V_sns_v_diff_min V_sns_v_diff_max]);
        title(strcat('Vertical Sensing Delta V(mV), Z=',num2str(k_index-1)))
        subplot(3,3,3)
        heatmap([0:net.Nj-1],[0:net.Ni],net.V_sns_h_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[V_sns_h_diff_min V_sns_h_diff_max]);
        title(strcat('Horizontal Sensing Delta V(mV), Z=',num2str(k_index-1)))
        if net.Nk>2
            k_index=2;
            subplot(3,3,4)
            heatmap([0:net.Nj],[0:net.Ni],net.V_node_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[V_node_diff_min V_node_diff_max]);
            title(strcat('Node Volatge Delta V(mV), Z=',num2str(k_index-1)))
            subplot(3,3,5)
            heatmap([0:net.Nj],[0:net.Ni-1],net.V_sns_v_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[V_sns_v_diff_min V_sns_v_diff_max]);
            title(strcat('Vertical Sensing Delta V(mV), Z=',num2str(k_index-1)))
            subplot(3,3,6)
            heatmap([0:net.Nj-1],[0:net.Ni],net.V_sns_h_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[V_sns_h_diff_min V_sns_h_diff_max]);
            title(strcat('Horizontal Sensing Delta V(mV), Z=',num2str(k_index-1)))
            k_index=3;
            subplot(3,3,7)
            heatmap([0:net.Nj],[0:net.Ni],net.V_node_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[V_node_diff_min V_node_diff_max]);
            title(strcat('Node Volatge Delta V(mV), Z=',num2str(k_index-1)))
            subplot(3,3,8)
            heatmap([0:net.Nj],[0:net.Ni-1],net.V_sns_v_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[V_sns_v_diff_min V_sns_v_diff_max]);
            title(strcat('Vertical Sensing Delta V(mV), Z=',num2str(k_index-1)))
            subplot(3,3,9)
            heatmap([0:net.Nj-1],[0:net.Ni],net.V_sns_h_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[V_sns_h_diff_min V_sns_h_diff_max]);
            title(strcat('Horizontal Sensing Delta V(mV), Z=',num2str(k_index-1)))
        end
        set(gcf, 'Position',[ 0.0794    0.0010    1.4576    0.7888]*1e3)
    end
end
% 4
fi=fi+1;
if any(fi==net.param.fi_req)
    fig_h(fi)=figure;
    subplot1=subplot(2,2,1);
    surf([0:net.Nj],[0:net.Ni-1],net.V_sns_v_NaN(:,:,1,1))
    title(strcat('Z=',num2str(0)))
    view(subplot1,[-50.6199999999998 41.3600000000001]);
    xlabel('X');ylabel('Y');zlabel('V_{DC}(mV)');
    set(gcf, 'Position',[372.2000    5.0000  690.4000  768.0000])
    ax = gca; ax.FontSize=13;
    subplot2=subplot(2,2,2);
    surf([0:net.Nj],[0:net.Ni-1],net.V_sns_v_NaN(:,:,end,1))
    title(strcat('Z=',num2str(net.Nk)))
    view(subplot2,[-50.6199999999998 41.3600000000001]);
    xlabel('X');ylabel('Y');zlabel('V_{DC}(mV)');
    ax = gca; ax.FontSize=13;
    subplot3=subplot(2,2,3);
    surf([0:net.Nj],[0:net.Ni-1],net.V_sns_v_diff(:,:,1))
    view(subplot3,[-50.6199999999998 41.3600000000001]);
    title(strcat('Z=',num2str(0)))
    xlabel('X');ylabel('Y');zlabel('\DeltaV_{PP}(mV)');
    ax = gca; ax.FontSize=13;
    subplot4=subplot(2,2,4);
    surf([0:net.Nj],[0:net.Ni-1],net.V_sns_v_diff(:,:,end))
    title(strcat('Z=',num2str(net.Nk)))
    xlabel('X');ylabel('Y');zlabel('\DeltaV_{PP}(mV)');
    view(subplot4,[-50.6199999999998 41.3600000000001]);
    ax = gca; ax.FontSize=13;
    set(gcf, 'Position',[0.1906    0.0146    1.1200    0.7680]*1e3)
end
% Delta V Location change (Z=0)

% 5 R and V vs time
if net.param.r_var_max>0
    fi=fi+1;
    if any(fi==net.param.fi_req)
        fig_h(fi)=figure;
        te=net.t_var(end)*1e3;
        subplot(2,4,1)
        hold on;grid on;
        plot(net.t_var*1e3,squeeze(r_var_t(1,:,:)),'linewidth',2)
        title('Artery1 R')
        xlabel('Time (ms)');ylabel('R(\Omega)');xlim([0 te]);
        ax = gca; ax.FontSize=13;
        subplot(2,4,2)
        hold on;grid on;
        plot(net.t_var*1e3,squeeze(net.V_sns_v_n(1:net.Ni,1,1,:)),'linewidth',2)
        title(strcat('\DeltaV at X=',num2str(0),',Z=0'));%xlim([0 te]);ylim([-1 0.5]);
        xlabel('Time (ms)');ylabel('\DeltaV')
        ax = gca; ax.FontSize=13;
        subplot(2,4,3)
        hold on;grid on;
        plot(net.t_var*1e3,squeeze(net.V_sns_v_n(1:net.Ni,net.param.r_var_n1j(1),1,:)),'linewidth',2)
        title(strcat('\DeltaV at X=',num2str(net.param.r_var_n1j(1)-1),',Z=0 (Artery1)'))
        xlabel('Time (ms)');ylabel('\DeltaV');%xlim([0 te]);ylim([-1 0.5]);
        ax = gca; ax.FontSize=13;
        subplot(2,4,4)
        hold on;grid on;
        plot(net.t_var*1e3,squeeze(net.V_sns_v_n(1:net.Ni,end,1,:)),'linewidth',2)
        title(strcat('\DeltaV at X=',num2str(net.Nj),',Z=0'));%xlim([0 te]);ylim([-1 0.5]);
        xlabel('Time (ms)');ylabel('\DeltaV')
        ax = gca; ax.FontSize=13;
        if S>1
            subplot(2,3,4)
            hold on;grid on;
            plot(net.t_var*1e3,squeeze(r_var_t(2,:,:)),'linewidth',2)
            title('Artery2 R')
            xlabel('Time (ms)');ylabel('R(\Omega)');%xlim([0 te]);
            ax = gca; ax.FontSize=13;
            subplot(2,3,5)
            hold on;grid on;
            plot(net.t_var*1e3,squeeze(net.V_sns_v_n(1:net.Ni,net.param.r_var_n1j(2),1,:)),'linewidth',2)
            title(strcat('\DeltaV at X=',num2str(net.param.r_var_n1j(2)),',Z=0 (Artery2)'))
            xlabel('Time (ms)');ylabel('\DeltaV');%xlim([0 te]);ylim([-1 0.5]);
            ax = gca; ax.FontSize=13;
        elseif net.Nk>2
            subplot(2,4,5)
            hold on;grid on;
            plot(net.t_var*1e3,squeeze(net.V_sns_v_n(1:net.Ni,1,2,:)),'linewidth',2)
            title(strcat('\DeltaV at X=',num2str(0),',Z=1'));%xlim([0 te]);ylim([-1 0.5]);
            xlabel('Time (ms)');ylabel('\DeltaV')
            ax = gca; ax.FontSize=13;
            subplot(2,4,6)
            hold on;grid on;
            plot(net.t_var*1e3,squeeze(net.V_sns_v_n(1:net.Ni,net.param.r_var_n1j(1)-1,2,:)),'linewidth',2)
            title(strcat('\DeltaV at X=',num2str(net.param.r_var_n1j(1)-1),',Z=1'));%xlim([0 te]);ylim([-1 0.5]);
            xlabel('Time (ms)');ylabel('\DeltaV')
            ax = gca; ax.FontSize=13;
            subplot(2,4,7)
            hold on;grid on;
            plot(net.t_var*1e3,squeeze(net.V_sns_v_n(1:net.Ni,net.param.r_var_n1j(1),2,:)),'linewidth',2)
            title(strcat('\DeltaV at X=',num2str(net.param.r_var_n1j(1)-1),',Z=1 (Artery1)'))
            xlabel('Time (ms)');ylabel('\DeltaV');%xlim([0 te]);ylim([-1 0.5]);
            ax = gca; ax.FontSize=13;
            subplot(2,4,8)
            hold on;grid on;
            plot(net.t_var*1e3,squeeze(net.V_sns_v([1 2 end-1 end],net.param.r_var_n1j(1),2,:)),'linewidth',2)
            title(strcat('\DeltaV at X=',num2str(net.Nj),',Z=1'));%xlim([0 te]);ylim([-1 0.5]);
            xlabel('Time (ms)');ylabel('\DeltaV')
            ax = gca; ax.FontSize=13;
        end
        set(gcf, 'Position',[ 0.0794    0.0010    1.4576    0.7888]*1e3)
    end
    
    %6
    fi=fi+1;
    if any(fi==net.param.fi_req)
        if net.param.r_var_dim==3
            fig_h(fi)=figure;
            te=net.t_var(end)*1e3;
            subplot(2,4,1)
            hold on;grid on;
            plot(net.t_var*1e3,squeeze(r_var_t(1,:,:)),'linewidth',2)
            title('Artery1 R')
            xlabel('Time (ms)');ylabel('R(\Omega)');xlim([0 te]);
            ax = gca; ax.FontSize=13;
            subplot(2,4,2)
            hold on;grid on;
            plot(net.t_var*1e3,squeeze(net.V_sns_v_n(net.curr.n1i+1:net.curr.n2i-2,1,1,:)),'linewidth',2)
            title(strcat('\DeltaV at X=',num2str(0),',Z=0'));%xlim([0 te]);ylim([-1 0.5]);
            xlabel('Time (ms)');ylabel('\DeltaV')
            ax = gca; ax.FontSize=13;
            subplot(2,4,3)
            hold on;grid on;
            plot(net.t_var*1e3,squeeze(net.V_sns_v_n(net.curr.n1i+1:net.curr.n2i-2,net.param.r_var_n1j(1),1,:)),'linewidth',2)
            title(strcat('\DeltaV at X=',num2str(net.param.r_var_n1j(1)-1),',Z=0 (Artery1)'))
            xlabel('Time (ms)');ylabel('\DeltaV');%xlim([0 te]);ylim([-1 0.5]);
            ax = gca; ax.FontSize=13;
            subplot(2,4,4)
            hold on;grid on;
            plot(net.t_var*1e3,squeeze(net.V_sns_v_n(net.curr.n1i+1:net.curr.n2i-2,end,1,:)),'linewidth',2)
            title(strcat('\DeltaV at X=',num2str(net.Nj),',Z=0'));%xlim([0 te]);ylim([-1 0.5]);
            xlabel('Time (ms)');ylabel('\DeltaV')
            ax = gca; ax.FontSize=13;
            subplot(2,4,5)
            hold on;grid on;
            plot(net.t_var*1e3,squeeze(net.V_sns_v_n(net.curr.n1i+1:net.curr.n2i-2,1,2,:)),'linewidth',2)
            title(strcat('\DeltaV at X=',num2str(0),',Z=1'));%xlim([0 te]);ylim([-1 0.5]);
            xlabel('Time (ms)');ylabel('\DeltaV')
            ax = gca; ax.FontSize=13;
            subplot(2,4,6)
            hold on;grid on;
            plot(net.t_var*1e3,squeeze(net.V_sns_v_n(net.curr.n1i+1:net.curr.n2i-2,net.param.r_var_n1j(1)-1,2,:)),'linewidth',2)
            title(strcat('\DeltaV at X=',num2str(net.param.r_var_n1j(1)-1),',Z=1'));%xlim([0 te]);ylim([-1 0.5]);
            xlabel('Time (ms)');ylabel('\DeltaV')
            ax = gca; ax.FontSize=13;
            subplot(2,4,7)
            hold on;grid on;
            plot(net.t_var*1e3,squeeze(net.V_sns_v_n(net.curr.n1i+1:net.curr.n2i-2,net.param.r_var_n1j(1),2,:)),'linewidth',2)
            title(strcat('\DeltaV at X=',num2str(net.param.r_var_n1j(1)-1),',Z=1 (Artery1)'))
            xlabel('Time (ms)');ylabel('\DeltaV');%xlim([0 te]);ylim([-1 0.5]);
            ax = gca; ax.FontSize=13;
            subplot(2,4,8)
            hold on;grid on;
            %plot(t_var*1e3,squeeze(V_sns_v_n(1:net.Ni,end,2,:)),'linewidth',2)
            plot(net.t_var*1e3,squeeze(net.V_sns_v([1 2 end-1 end],net.param.r_var_n1j(1),2,:)),'linewidth',2)
            title(strcat('\DeltaV at X=',num2str(net.Nj),',Z=1'));%xlim([0 te]);ylim([-1 0.5]);
            xlabel('Time (ms)');ylabel('\DeltaV')
            ax = gca; ax.FontSize=13;
            set(gcf, 'Position',[ 0.0794    0.0010    1.4576    0.7888]*1e3)
        end
    end
    
    if S==1 && net.param.r_var_dim==3
        % 7 PAT vs Y, PTT vs Y, PAT between the current electrodes
        fi=fi+1;
        if any(fi==net.param.fi_req)
            fig_h(fi)=figure;
            subplot(2,2,1)
            grid on;hold on;
            plot([net.curr.n1i+1:net.curr.n2i-2]-1,net.PAT_r_var(net.curr.n1i+1:net.curr.n2i-2)*1e3,'k','linewidth',2)
            plot([net.curr.n1i+1:net.curr.n2i-2]-1,squeeze(net.PAT_sns_v_n_k(net.curr.n1i+1:net.curr.n2i-2,net.param.r_var_n1j:3:end,2))*1e3,'linewidth',2)
            title('Z=1,spacing=1');
            xlabel('Y');ylabel('PAT (ms)')
            legend(["Artery";strcat(repmat("X=",floor((net.Ni+1-net.param.r_var_n1j)/3)+1,1),num2str([net.param.r_var_n1j:3:net.Ni+1]'))],'location','best')
            ax = gca; ax.FontSize=13;
            subplot(2,2,2)
            grid on;hold on;
            plot([net.curr.n1i+1:net.curr.n2i-2]-1,net.PAT_r_var(net.curr.n1i+1:net.curr.n2i-2)*1e3,'k','linewidth',2)
            plot([net.curr.n1i+1:net.curr.n2i-2]-1,squeeze(net.PAT_sns_v_n_k(net.curr.n1i+1:net.curr.n2i-2,net.param.r_var_n1j:3:end,1))*1e3,'linewidth',2)
            legend(["Artery";strcat(repmat("X=",floor((net.Ni+1-net.param.r_var_n1j)/3)+1,1),num2str([net.param.r_var_n1j:3:net.Ni+1]'))],'location','best')
            title('Z=0, spacing=1');
            xlabel('Y');ylabel('PAT (ms)')
            ax = gca; ax.FontSize=13;
            if spacing_arr>2
                subplot(2,2,3)
                p=3;
                grid on;hold on;
                plot([net.curr.n1i+1:net.curr.n2i-1-spacing_arr(p)]-1+floor(spacing_arr(p)/2),net.PAT_r_var([net.curr.n1i+1:net.curr.n2i-1-spacing_arr(p)]+floor(spacing_arr(p)/2))*1e3,'k','linewidth',2)
                plot([net.curr.n1i+1:net.curr.n2i-1-spacing_arr(p)]-1+floor(spacing_arr(p)/2),net.PAT_sns_v_space_n_k1{p}(net.curr.n1i+1:net.curr.n2i-1-spacing_arr(p),net.param.r_var_n1j:3:end)*1e3,'linewidth',2)
                title(strcat('Z=0, spacing=',num2str(spacing_arr(p))));
                xlabel('Y');ylabel('PAT (ms)')
                legend(["Artery";strcat(repmat("X=",floor((net.Ni+1-net.param.r_var_n1j)/3)+1,1),num2str([net.param.r_var_n1j:3:net.Ni+1]'))],'location','best')
                ax = gca; ax.FontSize=13;
            end
            subplot(2,2,4)
            grid on;hold on;
            plot([0 net.Nj],[net.PTT_r_var net.PTT_r_var]*1e3,'k','linewidth',2)
            plot([0:net.Nj],net.PTT_sns_v_space*1e3,'linewidth',2)
            title('Z=0');
            xlabel('X');ylabel('PTT between current (ms)')
            legend(["Artery";strcat(repmat("spacing=",length(spacing_arr),1),num2str(spacing_arr'))],'location','best')
            ax = gca; ax.FontSize=13;
            set(gcf, 'Position',[ 120 2  1086 682])
        end
        
        % 8
        fi=fi+1;
        if any(fi==net.param.fi_req)
            fig_h(fi)=figure;
            subplot(2,2,1)
            grid on;hold on;
            plot([net.curr.n1i+1:net.curr.n2i-2]-1,net.PAT_r_var(net.curr.n1i+1:net.curr.n2i-2)*1e3,'k','linewidth',2)
            plot([net.curr.n1i+1:net.curr.n2i-2]-1,squeeze(net.PAT_sns_v_n_k(net.curr.n1i+1:net.curr.n2i-2,net.param.r_var_n1j,:))*1e3,'linewidth',2)
            title(strcat('Variable Z, X=',num2str(net.curr.n1j),',spacing=1'));
            xlabel('Y');ylabel('PAT (ms)')
            legend(["Artery";strcat(repmat("Z=",net.Nk+1,1),num2str([0:net.Nk]'))],'location','best')
            ax = gca; ax.FontSize=13;
            subplot(2,2,2)
            grid on;hold on;
            plot([net.curr.n1i+1:net.curr.n2i-2]-1,net.PAT_r_var(net.curr.n1i+1:net.curr.n2i-2)*1e3,'k','linewidth',2)
            plot([net.curr.n1i+1:net.curr.n2i-2]-1,squeeze(net.PAT_sns_v_n_k(net.curr.n1i+1:net.curr.n2i-2,net.param.r_var_n1j:2:end,1))*1e3,'linewidth',2)
            legend(["Artery";strcat(repmat("X=",floor((net.Ni+1-net.param.r_var_n1j)/2)+1,1),num2str([net.param.r_var_n1j:2:net.Ni+1]'))],'location','best')
            title('Variable X, Z=0, spacing=1');
            xlabel('Y');ylabel('PAT (ms)')
            ax = gca; ax.FontSize=13;
            subplot(2,2,3)
            grid on;hold on;
            p=1;
            plot([net.curr.n1i+1:net.curr.n2i-1-net.spacing_arr(p)]-1+floor(net.spacing_arr(p)/2),net.PAT_r_var([net.curr.n1i+1:net.curr.n2i-1-net.spacing_arr(p)]+floor(net.spacing_arr(p)/2))*1e3,'k','linewidth',2)
            for p=1:4:length(net.spacing_arr)
                plot([net.curr.n1i+1:net.curr.n2i-1-net.spacing_arr(p)]-1+floor(net.spacing_arr(p)/2),net.PAT_sns_v_space_n_k1{p}(net.curr.n1i+1:net.curr.n2i-1-net.spacing_arr(p),net.param.r_var_n1j)*1e3,'linewidth',2)
            end
            title(strcat('Variable Spacing, Z=0, X=',num2str(net.curr.n1j)));
            xlabel('Y');ylabel('PAT (ms)')
            legend(["Artery";strcat(repmat("Spacing=",floor(length(1:3:length(net.spacing_arr))),1),num2str([1:3:length(net.spacing_arr)]'))],'location','best')
            ax = gca; ax.FontSize=13;
            set(gcf, 'Position',[ 120 2  1086 682])
        end
        
        % 9 PAT over XY plan
        fi=fi+1;
        if any(fi==net.param.fi_req)
            fig_h(fi)=figure;
            subplot1=subplot(2,1,1);
            surf(squeeze(net.PAT_sns_v_n_k(net.curr.n1i+1:net.curr.n2i-2,:,2))*1e3)
            xlabel('X');ylabel('Y');zlabel('PAT (ms)');
            title('Z=1');
            view(subplot1,[-49.6599999999997 22.8000000000001]);
            ax = gca; ax.FontSize=13;
            subplot2=subplot(2,1,2);
            surf(squeeze(net.PAT_sns_v_n_k(net.curr.n1i+1:net.curr.n2i-2,:,1))*1e3)
            xlabel('X');ylabel('Y');zlabel('PAT (ms)');
            title('Z=0');
            view(subplot2,[-60.8599999999998 17.6800000000001]);
            set(gcf, 'Position',[372.2000    5.0000  690.4000  768.0000])
            ax = gca; ax.FontSize=13;
        end
    end
end

% if net.param.write_figures==1
%     file_str=strcat(res_folder,'N',num2str(Ni));
%     file_str=strcat(file_str,'_i',num2str(net.param.i_src_n1i),'_',num2str(net.param.i_src_n2i),'_',num2str(net.param.i_src_n1j),'_num',num2str(net.param.i_num),'_ip',num2str(net.param.i_opp),'_a',num2str(net.curr(1).area),'_rj',num2str(net.param.r_var_n1j),'_rk',num2str(net.param.r_var_n1k),'_rs',num2str(net.param.r_var_seg));
%     for i=1:length(fi)
%         if any(fi==net.param.fi_req)
%             figure(fi(i))
%             print(strcat(file_str,'_',num2str(i),'.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
%         end
%     end
% end
