function plot_freq_compact_fn(net_arr)


Av=[1];

fig1=openfig('G:\My Drive\ESP lab\projects\TBioCAS-Bio-Z-Sim\bassem\bioz-sim-freq-Nov2019\results\2020-09-18\freq_r_compact.fig');
fig2=openfig('G:\My Drive\ESP lab\projects\TBioCAS-Bio-Z-Sim\bassem\bioz-sim-freq-Nov2019\results\2020-09-18\freq_i_compact.fig');

n=length(net_arr);


for j=1:length(net_arr)
    k=length(net_arr{j});
    curr_arr(j,:)=cat(1,net_arr{j}.curr);
    s_arr(j,:,:)=cat(1,net_arr{j}.s);
end

% 3D: 1=Paramter, 2=frequency, 3=sensing voltage
Vdc_curr_i=reshape(cat(1,curr_arr.Vdc_e_i),n,k);
dVpp_curr_i=reshape(cat(1,curr_arr.dVpp_e_i),n,k);
Vdc_curr_r=reshape(cat(1,curr_arr.Vdc_e_r),n,k);
dVpp_curr_r=reshape(cat(1,curr_arr.dVpp_e_r),n,k);
Vdc_curr=reshape(cat(1,curr_arr.Vdc_e),n,k);
dVpp_curr=reshape(cat(1,curr_arr.dVpp_e),n,k);

Vdc_s_i=reshape(cat(1,s_arr.Vdc_e_i),n,k,length(net_arr{1}(1).s));
dVpp_s_i=reshape(cat(1,s_arr.dVpp_e_i),n,k,length(net_arr{1}(1).s));
Vdc_s_r=reshape(cat(1,s_arr.Vdc_e_r),n,k,length(net_arr{1}(1).s));
dVpp_s_r=reshape(cat(1,s_arr.dVpp_e_r),n,k,length(net_arr{1}(1).s));
Vdc_s=reshape(cat(1,s_arr.Vdc_e),n,k,length(net_arr{1}(1).s));
dVpp_s=reshape(cat(1,s_arr.dVpp_e),n,k,length(net_arr{1}(1).s));
PAT_s=reshape(cat(1,s_arr.PAT_e),n,k,length(net_arr{1}(1).s));




%for i=1:size(Vdc_s_r,3)
    %mode_color=net_arr{1}(1).param.color_arr(7,:);
    mode_color=[0 0 0];
    xlim_freq=[net_arr{1}(1).f_var(1,1)*1e-3 net_arr{1}(1).f_var(end,1)*1e-3*1.01];
    figure(fig1)
    subplot(2,2,1)    
    yyaxis right
    grid on;hold on;
    h2=plot(net_arr{1}(1).f_var(:,1)*1e-3,squeeze(Vdc_s_r(1,:,1))*Av*1e3,'-','linewidth',2,'Color',mode_color,'DisplayName','Model');
    xlabel({'Frequency (kHz)'})
    ylabel('{\it{V}}_{dc} (mV)')
    ax = gca; ax.FontSize=15;
    set(gca, 'FontName', 'Times New Roman');   
    ylim([0 140]);xlim(xlim_freq)
    AX=gca;set(AX,'ycolor',mode_color)

    subplot(2,2,2)    
    yyaxis right
    grid on;hold on;
    h2=plot(net_arr{1}(1).f_var(:,1)*1e-3,squeeze(Vdc_s_r(1,:,2))*Av*1e3,'-','linewidth',2,'Color',mode_color,'DisplayName','Model');
    xlabel({'Frequency (kHz)'})
    ylabel('{\it{V}}_{dc} (mV)')
    ax = gca; ax.FontSize=15;
    set(gca, 'FontName', 'Times New Roman');   
    ylim([0 140]);xlim(xlim_freq)
    AX=gca;set(AX,'ycolor',mode_color)
    
    subplot(2,2,3)
    yyaxis right
    grid on;hold on;
    plot(net_arr{1}(1).f_var(:,1)*1e-3,squeeze(dVpp_s_r(1,:,1))*Av*1e6,'-','linewidth',2,'Color',mode_color,'DisplayName','Model');
    xlabel({'Frequency (kHz)'})
    ylabel('{\it{\DeltaV}}_{pp} (\muV)')
    ax = gca; ax.FontSize=15;
    set(gca, 'FontName', 'Times New Roman');    
    ylim([0 80]);xlim(xlim_freq)
    AX=gca;set(AX,'ycolor',mode_color)
    subplot(2,2,4)
    yyaxis right
    grid on;hold on;
    plot(net_arr{1}(1).f_var(:,1)*1e-3,squeeze(dVpp_s_r(1,:,2))*Av*1e6,'-','linewidth',2,'Color',mode_color,'DisplayName','Model');
    xlabel({'Frequency (kHz)'})
    ylabel('{\it{\DeltaV}}_{pp} (\muV)')
    ax = gca; ax.FontSize=15;
    set(gca, 'FontName', 'Times New Roman');    
    ylim([0 80]);xlim(xlim_freq)
    AX=gca;set(AX,'ycolor',mode_color)
    
    figure(fig2)
    subplot(2,2,1)
    yyaxis right
    grid on;hold on;
    plot(net_arr{1}(1).f_var(:,1)*1e-3,squeeze((Vdc_s_i(1,:,1)))*Av*1e3,'-','linewidth',2,'Color',mode_color,'DisplayName','Model');
    xlabel({'Frequency (kHz)'})
    ylabel('{\it{V}}_{dc} (mV)')
    ax = gca; ax.FontSize=15;
    set(gca, 'FontName', 'Times New Roman');
    ylim([-5 0]);xlim(xlim_freq)
    AX=gca;set(AX,'ycolor',mode_color)
    subplot(2,2,2)
    yyaxis right
    grid on;hold on;
    plot(net_arr{1}(1).f_var(:,1)*1e-3,squeeze((Vdc_s_i(1,:,2)))*Av*1e3,'-','linewidth',2,'Color',mode_color,'DisplayName','Model');
    xlabel({'Frequency (kHz)'})
    ylabel('{\it{V}}_{dc} (mV)')
    ax = gca; ax.FontSize=15;
    set(gca, 'FontName', 'Times New Roman');    
    ylim([-5 0]);xlim(xlim_freq)
    AX=gca;set(AX,'ycolor',mode_color)
    
    subplot(2,2,3)
    yyaxis right
    grid on;hold on;
    plot(net_arr{1}(1).f_var(:,1)*1e-3,squeeze((dVpp_s_i(1,:,1)))*Av*1e6*-1,'-','linewidth',2,'Color',mode_color,'DisplayName','Model');
    xlabel({'Frequency (kHz)'})
    ylabel('{\it{\DeltaV}}_{pp} (\muV)')
    ax = gca; ax.FontSize=15;
    set(gca, 'FontName', 'Times New Roman');    
    ylim([0 4]);xlim(xlim_freq)    
    AX=gca;set(AX,'ycolor',mode_color)
    subplot(2,2,4)    
    yyaxis right
    grid on;hold on;
    plot(net_arr{1}(1).f_var(:,1)*1e-3,squeeze((dVpp_s_i(1,:,2)))*Av*1e6*-1,'-','linewidth',2,'Color',mode_color,'DisplayName','Model');
    xlabel({'Frequency (kHz)'})
    ylabel('{\it{\DeltaV}}_{pp} (\muV)')
    ax = gca; ax.FontSize=15;
    set(gca, 'FontName', 'Times New Roman');   
    ylim([0 4]);xlim(xlim_freq)
    AX=gca;set(AX,'ycolor',mode_color)
    set(gcf,'position',[1.0000    1.0000  998.4000  641.2000]) 
%end

if net_arr{1}(1).param.write_figures==1 && net_arr{1}(1).param.debug==0
    file_str=strcat(net_arr{1}(1).param.path);
    figure(fig1)
    print(strcat(file_str,'freq_r_compact.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
    %savefig(strcat(file_str,'freq_r.fig')); %<-Save as PNG with 300 DP
    figure(fig2)
    print(strcat(file_str,'freq_i_compact.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
    %savefig(strcat(file_str,'freq_i.fig')); %<-Save as PNG with 300 DP
end
