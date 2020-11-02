function plot_freq_fn(net_arr)


Av=[1];

fig1=openfig('G:\My Drive\ESP lab\projects\TBioCAS-Bio-Z-Sim\bassem\bioz-sim-freq-Nov2019\results\2020-09-11\freq_r.fig');
fig2=openfig('G:\My Drive\ESP lab\projects\TBioCAS-Bio-Z-Sim\bassem\bioz-sim-freq-Nov2019\results\2020-09-11\freq_i.fig');

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




for i=1:size(Vdc_s_r,3)
    figure(fig1)
    subplot(2,3,3)
    grid on;hold on;
    plot(net_arr{1}(1).f_var(:,1)*1e-3,squeeze(Vdc_s_r(1,:,i))*Av*1e3,'-','linewidth',2,'Color',net_arr{1}(1).param.color_arr(i,:));
    xlabel({'Frequency (kHz)'})
    ylabel('{\it{V}}_{dc} (mV)')
    ax = gca; ax.FontSize=15;
    set(gca, 'FontName', 'Times New Roman');
    legend('Re({\it{V}}_{1}) Model','Re({\it{V}}_{2}) Model','fontsize',12)
    ylim([0 90])
    
    subplot(2,3,6)
    grid on;hold on;
    plot(net_arr{1}(1).f_var(:,1)*1e-3,squeeze(dVpp_s_r(1,:,i))*Av*1e6,'-','linewidth',2,'Color',net_arr{1}(1).param.color_arr(i,:));
    xlabel({'Frequency (kHz)'})
    ylabel('{\it{\DeltaV}}_{pp} (\muV)')
    ax = gca; ax.FontSize=15;
    set(gca, 'FontName', 'Times New Roman');
    legend('Re({\it{V}}_{1}) Model','Re({\it{V}}_{2}) Model','fontsize',12)
    ylim([0 60])
    
    figure(fig2)
    subplot(2,3,3)
    grid on;hold on;
    plot(net_arr{1}(1).f_var(:,1)*1e-3,squeeze((Vdc_s_i(1,:,i)))*Av*1e3,'-','linewidth',2,'Color',net_arr{1}(1).param.color_arr(i,:));
    xlabel({'Frequency (kHz)'})
    ylabel('{\it{V}}_{dc} (mV)')
    ax = gca; ax.FontSize=15;
    set(gca, 'FontName', 'Times New Roman');
    legend('Im({\it{V}}_{1}) Model','Im({\it{V}}_{2}) Model','fontsize',12)
    ylim([-2.5 0.5])
    
    subplot(2,3,6)
    grid on;hold on;
    plot(net_arr{1}(1).f_var(:,1)*1e-3,squeeze((dVpp_s_i(1,:,i)))*Av*1e6*-1,'-','linewidth',2,'Color',net_arr{1}(1).param.color_arr(i,:));
    xlabel({'Frequency (kHz)'})
    ylabel('{\it{\DeltaV}}_{pp} (\muV)')
    ax = gca; ax.FontSize=15;
    set(gca, 'FontName', 'Times New Roman');
    legend('Im({\it{V}}_{1}) Model','Im({\it{V}}_{2}) Model','fontsize',12)
    ylim([0 4])
   %set(gca,'position',[0.6916    0.1185    0.2134    0.3326]) 
end

if net_arr{1}(1).param.write_figures==1 && net_arr{1}(1).param.debug==0
    file_str=strcat(net_arr{1}(1).param.path);
    figure(fig1)
    print(strcat(file_str,'freq_r.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
    %savefig(strcat(file_str,'freq_r.fig')); %<-Save as PNG with 300 DP
    figure(fig2)
    print(strcat(file_str,'freq_i.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
    %savefig(strcat(file_str,'freq_i.fig')); %<-Save as PNG with 300 DP
end
