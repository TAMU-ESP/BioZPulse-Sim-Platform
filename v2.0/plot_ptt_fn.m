function plot_ptt_fn(net_arr)




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



% if net_arr{1}(1).param.write_figures==1 && net_arr{1}(1).param.debug==0
%     file_str=strcat(net_arr{1}(1).param.path);
%     print(strcat(file_str,'Paramter.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
%     %savefig(strcat(file_str,'freq.fig')); %<-Save as PNG with 300 DP
% end
% 
figure
x=[1:8]*8;
xfit=[x(1):0.1:x(end)];
p=spline(x,PAT_s(:,1,1));
PAT_s_1_fit=ppval(p,xfit);
p=spline(x,PAT_s(:,1,2));
PAT_s_2_fit=ppval(p,xfit);


    
subplot(2,1,1)
grid on;hold on;
plot(xfit,PAT_s_1_fit*1e3,'-','linewidth',2);
plot(xfit,PAT_s_2_fit*1e3,'-','linewidth',2);
xlabel('Horizontal Sensing Location {\it{PI}}_X (mm)')
ylabel('{\it{T}}_d (ms)')
ax = gca; ax.FontSize=13;xlim([8 64])
set(gca, 'FontName', 'Times New Roman');  
legend('{\it{V1}}','{\it{V2}}','fontsize',13,'location','northeast')
ylim([0 20])

subplot(2,1,2)
grid on;hold on;
plot(xfit,(PAT_s_2_fit-PAT_s_1_fit)*1e3,'-','linewidth',2);
xlabel('Horizontal Sensing Location {\it{PI}}_X (mm)')
ylabel('PTT (ms)')
ax = gca; ax.FontSize=13;
set(gca, 'FontName', 'Times New Roman');
ylim([0 15]);xlim([8 64])
set(gcf,'position',[488.2000  307.8000  534.8000  454.4000])

if net_arr{1}(1).param.write_figures==1 && net_arr{1}(1).param.debug==0
    file_str=strcat(net_arr{1}(1).param.path);
    print(strcat(file_str,'ptt.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
    %savefig(strcat(file_str,'freq.fig')); %<-Save as PNG with 300 DP
end