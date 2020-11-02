function net=spice_netlist_3d_plot_ac_fn(net)

%   Bassem Ibrahim
%   2019.08.22
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation
%   Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling,
%   IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan.
%
%   function net=spice_netlist_3d_plot_ac_fn(net)
%   Plot Output Results in Frequency domain


S=net.param.S;
spacing_arr=net.spacing_arr;


V_sns_v_max=max(max(max(abs(squeeze(net.V_sns_v_NaN(:,:,1,:))))));
V_sns_h_max=max(max(max(abs(squeeze(net.V_sns_h_NaN(:,:,1,:))))));
V_sns_max=max(V_sns_h_max,V_sns_v_max);
V_node_max=max(max(max(abs(squeeze(net.V_node_NaN(:,:,1,:))))));

fi=0;

% 1 All V Heatmaps
fi=fi+1;
if any(fi==net.param.fi_req)
    fig_h(fi)=figure;
    %suptitle(strcat('Figure ',num2str(fi)))
    subplot(3,2,1)
    heatmap([0:net.Nj],[0:net.Ni],squeeze(net.node_arr(:,:,1)),'Colormap',white,'MissingDataLabel','Elec.');
    xlabel('X');ylabel('Y')
    title('Electrodes Location (Black squares), Z=0')
    if net.Nk>0
        subplot(3,2,3)
        heatmap([0:net.Nj],[0:net.Nk],squeeze(net.node_arr(1,:,:))','Colormap',white,'MissingDataLabel','Artery');
        xlabel('X');ylabel('Z')
        title('Artery Location (Black squares), Z=0')
    end
    subplot(3,2,5)
    heatmap([0:net.Nj],[0:net.Ni],real(net.V_node(:,:,1,1)),'Colormap',parula,'ColorLimits',[-1*V_node_max V_node_max]);
    xlabel('X');ylabel('Y')
    title('Node V (mV), Z=0')
    
    k_index=1;
    subplot(3,2,2)
    heatmap([0:net.Nj],[0:net.Ni-1],real(net.V_sns_v_NaN(:,:,k_index,1)),'Colormap',parula,'ColorLimits',[-1*V_sns_max V_sns_max]);
    xlabel('X');ylabel('Y')
    title(strcat('V_{DC}(mV), Z=',num2str(k_index-1)));
    if net.Nk>2
        k_index=2;
        subplot(3,2,4)
        heatmap([0:net.Nj],[0:net.Ni-1],real(net.V_sns_v_NaN(:,:,k_index,1)),'Colormap',parula,'ColorLimits',[-1*V_sns_max V_sns_max]);
        xlabel('X');ylabel('Y')
        title(strcat('V_{DC}(mV), Z=',num2str(k_index-1)))
        k_index=3;
        subplot(3,2,6)
        heatmap([0:net.Nj],[0:net.Ni-1],real(net.V_sns_v_NaN(:,:,k_index,1)),'Colormap',parula,'ColorLimits',[-1*V_sns_max V_sns_max]);
        xlabel('X');ylabel('Y')
        title(strcat('V_{DC}(mV), Z=',num2str(k_index-1)))
    end
    set(gcf, 'Position',[93.8000   20.2000  746.4000  742.4000])
end


% Location change (Z=0)
V_sns_v_loc=squeeze(net.V_sns_v(net.curr.n1i(1)+1:net.curr.n2i(1)-2,net.curr.n1j(1),1,1));
for p=1:length(net.spacing_arr)
    V_sns_v_spacing_p=net.V_sns_v_space{p}(net.curr.n1i(1)+1:net.curr.n2i(1)-1-net.spacing_arr(p),net.curr.n1j(1),1,1);
    V_sns_v_spacing(p)=V_sns_v_spacing_p(round(length(V_sns_v_spacing_p)/2));
end
% 2
fi=fi+1;
if any(fi==net.param.fi_req)
    fig_h(fi)=figure;
    subplot(3,3,1)
    hold on;grid on;
    plot([0:net.Ni-1],squeeze(net.V_sns_v_NaN(:,net.curr.n1j(1)-1,1,1)),'linewidth',2)
    plot([0:net.Ni-1],squeeze(net.V_sns_v_NaN(:,round(net.Nj/4),1,1)),'linewidth',2)
    plot([0:net.Ni-1],squeeze(net.V_sns_v_NaN(:,1,1,1)),'linewidth',2)
    xlabel('Y')
    ylabel('Vertical Sensing V(mV)')
    ax = gca; ax.FontSize=13;
    
    subplot(3,3,2)
    hold on;grid on;
    plot([1:net.Nj+1]-1,squeeze([net.V_sns_v(min(round(mean([mean(net.curr.n1i),mean(net.curr.n2i)])),net.Ni-1),:,1,1)]),'linewidth',2)
    plot([1:net.Nj+1]-1,squeeze([net.V_sns_v(min(round(mean([max(net.curr.n1i)+1])),net.Ni-1),:,1,1)]),'linewidth',2)
    xlabel('X')
    ylabel('Vertical Sensing V(mV)')
    ax = gca; ax.FontSize=13;
    
    subplot(3,3,3)
    hold on;grid on;
    plot([1:size(net.V_sns_v,3)]-1,squeeze([net.V_sns_v(min(max(net.curr.n1i)+1,net.Ni-1),min(net.curr.n1j(1),net.Nj-1),:,1)]),'linewidth',2)
    plot([1:size(net.V_sns_v,3)]-1,squeeze([net.V_sns_v_NaN(round(mean([net.curr.n1i(1),net.curr.n2i(1)])),min(net.curr.n1j(1),net.Nj-1),:,1)]),'linewidth',2)
    xlabel('Z')
    ylabel('Vertical Sensing V(mV)')
    ax = gca; ax.FontSize=13;
    if ~isempty(net.spacing_arr)
        subplot(3,3,6)
        hold on;grid on;
        plot(net.spacing_arr,V_sns_v_spacing,'linewidth',2)
        xlabel('spacing of Sensing Electrodes')
        ylabel('Vertical Sensing V(mV)')
        ax = gca; ax.FontSize=13;
    end
    if length(net.f_var(:,1))>1
        subplot(3,3,7)
        hold on;grid on;
        if isempty(net.param.s_n1i)==1
            plot(net.f_var(:,1),squeeze([real(net.V_sns_v(min(max(net.curr.n1i)+1,net.Ni-1),min(max(net.curr.n1j),net.Nj-1),1,:))]),'linewidth',2)
        else
            plot(net.f_var(:,1),real(net.V_s),'linewidth',2)
        end
        xlabel('Frequency(Hz)')
        ylabel('Vertical Sensing V(mV) [Real]')
        ax = gca; ax.FontSize=13;
        subplot(3,3,8)
        hold on;grid on;
        if isempty(net.param.s_n1i)==1
            plot(net.f_var(:,1),squeeze([imag(net.V_sns_v(min(max(net.curr.n1i)+1,net.Ni-1),min(max(net.curr.n1j),net.Nj-1),1,:))]),'linewidth',2)
        else
            plot(net.f_var(:,1),imag(net.V_s),'linewidth',2)
        end
        xlabel('Frequency(Hz)')
        ylabel('Vertical Sensing V(mV) [Imag.]')
        ax = gca; ax.FontSize=13;
    end
    subplot(3,3,9)
    hold on;grid on;
    if isempty(net.param.s_n1i)==1
        plot(squeeze([real(net.V_sns_v(min(max(net.curr.n1i)+1,net.Ni-1),min(max(net.curr.n1j),net.Nj-1),1,:))]),-1*squeeze([(imag(net.V_sns_v(min(max(net.curr.n1i)+1,net.Ni-1),min(max(net.curr.n1j),net.Nj-1),1,:)))]),'-o','linewidth',2)
    else
        plot(real(net.V_s),imag(net.V_s)*-1,'-o','linewidth',2)
    end
    xlabel('Vertical Sensing V(mV) [Real]')
    ylabel('Vertical Sensing V(mV) [Imag.]')
    ax = gca; ax.FontSize=13;
    set(gcf, 'Position',[ 0.0794    0.0010    1.4576    0.7888]*1e3)
end

% 3 Delta V SNS Z change Heatmaps
if sum(net.param.r_var_max>0)~=0 && net.param.mode==0
    fi=fi+1;
    if any(fi==net.param.fi_req)
        fig_h(fi)=figure;
        k_index=1;
        subplot(3,3,1)
        heatmap([0:net.Nj],[0:net.Ni],net.V_node_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[-1*V_node_diff_max V_node_diff_max]);
        title(strcat('Node Volatge Delta V(mV), Z=',num2str(k_index-1)))
        subplot(3,3,2)
        heatmap([0:net.Nj],[0:net.Ni-1],net.V_sns_v_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[-1*V_sns_diff_max V_sns_diff_max]);
        title(strcat('Vertical Sensing Delta V(mV), Z=',num2str(k_index-1)))
        subplot(3,3,3)
        heatmap([0:net.Nj-1],[0:net.Ni],net.V_sns_h_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[-1*V_sns_diff_max V_sns_diff_max]);
        title(strcat('Horizontal Sensing Delta V(mV), Z=',num2str(k_index-1)))
        if net.Nk>2
            k_index=2;
            subplot(3,3,4)
            heatmap([0:net.Nj],[0:net.Ni],net.V_node_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[-1*V_node_diff_max V_node_diff_max]);
            title(strcat('Node Volatge Delta V(mV), Z=',num2str(k_index-1)))
            subplot(3,3,5)
            heatmap([0:net.Nj],[0:net.Ni-1],net.V_sns_v_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[-1*V_sns_diff_max V_sns_diff_max]);
            title(strcat('Vertical Sensing Delta V(mV), Z=',num2str(k_index-1)))
            subplot(3,3,6)
            heatmap([0:net.Nj-1],[0:net.Ni],net.V_sns_h_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[-1*V_sns_diff_max V_sns_diff_max]);
            title(strcat('Horizontal Sensing Delta V(mV), Z=',num2str(k_index-1)))
            k_index=3;
            subplot(3,3,7)
            heatmap([0:net.Nj],[0:net.Ni],net.V_node_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[-1*V_node_diff_max V_node_diff_max]);
            title(strcat('Node Volatge Delta V(mV), Z=',num2str(k_index-1)))
            subplot(3,3,8)
            heatmap([0:net.Nj],[0:net.Ni-1],net.V_sns_v_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[-1*V_sns_diff_max V_sns_diff_max]);
            title(strcat('Vertical Sensing Delta V(mV), Z=',num2str(k_index-1)))
            subplot(3,3,9)
            heatmap([0:net.Nj-1],[0:net.Ni],net.V_sns_h_diff(:,:,k_index),'Colormap',parula,'ColorLimits',[-1*V_sns_diff_max V_sns_diff_max]);
            title(strcat('Horizontal Sensing Delta V(mV), Z=',num2str(k_index-1)))
        end
        set(gcf, 'Position',[ 0.0794    0.0010    1.4576    0.7888]*1e3)
    end
end
%4
fi=fi+1;
if any(fi==net.param.fi_req)
    fig_h(fi)=figure;
    subplot1=subplot(1,2,1);
    surf([0:net.Nj],[0:net.Ni-1],real(net.V_sns_v_NaN(:,:,1,1)))
    title(strcat('Z=',num2str(0)))
    view(subplot1,[-50.6199999999998 41.3600000000001]);
    xlabel('X');ylabel('Y');zlabel('V_{DC}(mV)');
    set(gcf, 'Position',[372.2000    5.0000  690.4000  768.0000])
    ax = gca; ax.FontSize=13;
    subplot2=subplot(1,2,2);
    surf([0:net.Nj],[0:net.Ni-1],real(net.V_sns_v_NaN(:,:,end,1)))
    title(strcat('Z=',num2str(net.Nk)))
    view(subplot2,[-50.6199999999998 41.3600000000001]);
    xlabel('X');ylabel('Y');zlabel('V_{DC}(mV)');
    ax = gca; ax.FontSize=13;
    set(gcf, 'Position',[793.8000   25.0000  730.4000  261.6000])
end

if net.param.write_figures==1
    file_str=strcat(res_folder,'N',num2str(Ni));
    file_str=strcat(file_str,'_i',num2str(net.param.i_src_n1i),'_',num2str(net.param.i_src_n2i),'_',num2str(net.param.i_src_n1j),'_num',num2str(net.param.i_num),'_ip',num2str(net.param.i_opp),'_a',num2str(net.curr(1).area),'_rj',num2str(net.param.r_var_n1j),'_rk',num2str(net.param.r_var_n1k),'_rs',num2str(net.param.r_var_seg));
    for i=1:length(fi)
        if any(fi==net.param.fi_req)
            figure(fi(i))
            print(strcat(file_str,'_',num2str(i),'.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
        end
    end
end
