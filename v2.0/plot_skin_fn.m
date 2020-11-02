function plot_skin_fn(net_arr)


Av=[1];

fig1=openfig('.\skin-elec-model\results\2020-09-10\all.fig');

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




for i=size(Vdc_s_r,3):-1:1
    figure(fig1)
    subplot(2,1,1)
    [old_legend]=findobj(gcf, 'Type', 'Legend');
    h1 = findobj(gca,'Type', 'line');
    grid on;hold on;
    h2=plot(net_arr{1}(1).f_var(:,1)*1e-3,squeeze(Vdc_curr_r)'*Av,'--','linewidth',2,'Color',net_arr{1}(1).param.color_arr(1,:));
    xlabel({'Frequency (kHz)'})
    ylabel('Re({\it{V}}_{se}) (V)')
    ax = gca; ax.FontSize=12;
    set(gca, 'FontName', 'Times New Roman');
    legend([h1([2,1]);h2],[old_legend(1).String,'Model'],'fontsize',12)
    %ylim([10 80])
    
    figure(fig1)
    subplot(2,1,2)
    [old_legend]=findobj(gcf, 'Type', 'Legend');
    h1 = findobj(gca,'Type', 'line');
    grid on;hold on;
    h2=plot(net_arr{1}(1).f_var(:,1)*1e-3,squeeze(Vdc_curr_i)'*Av,'--','linewidth',2,'Color',net_arr{1}(1).param.color_arr(2,:));
    xlabel({'Frequency (kHz)'})
    ylabel('Im({\it{\DeltaV}}_{se}) (V)')
    ax = gca; ax.FontSize=12;
    set(gca, 'FontName', 'Times New Roman');
    legend([h1([2,1]);h2],[old_legend(2).String],'fontsize',12,'location','SouthEast')
    %ylim([20 50])
    set(gca,'position',[0.1300    0.1382    0.7750    0.3129])
    
end

if net_arr{1}(1).param.write_figures==1 && net_arr{1}(1).param.debug==0
    file_str=strcat(net_arr{1}(1).param.path);
    figure(fig1)
    print(strcat(file_str,'skin.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
    %savefig(strcat(file_str,'freq_r.fig')); %<-Save as PNG with 300 DP   
end

figure(fig2)
