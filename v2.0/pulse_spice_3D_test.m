clc
clear
close all
fclose('all');
set(0,'DefaultFigureWindowStyle','docked')

net.param.write_netlist='eit_3d_test';
net.param.write_figures=1;
net.param.write_mat=1;
net.param.plot_eit=1;
net.param.debug=1;
net.param.date=datestr(now,'yyyy-mm-dd');
net.param.path=strcat('./results/',net.param.date,'_',net.param.write_netlist,'/');

%n_arr=[1:1:21];
%m_arr=[1:1:21];

%Tissue  Conduct(S/m)Permittivity
%Skin    2.04E-4     1.13E+3
%Bone    2.0E-2      5.0E+2
%Fat     4.3E-2      9.12E+2
%Muscles 3.4E-1      2.59E+4
%Blood   7.0E-1      5.25E+3

% Tissue objects location and RC
% bone,muscle,muscle,muscle,muscle,blood,blood
center_z_arr=[0.5,0.20,0.725,0.625,0.65,0.425,0.425];
center_x_arr=[0.5,0.5,0.6,0.725,0.2,0.275,0.10];
center_y_arr=[0.5,0.5,0.5,0.5,0.5,0.5,0.5];
len_z_arr=[0.15,0.275,0.15,0.35,0.200,0.00125,0.1];
len_x_arr=[0.15,0.725,0.3,0.15,0.2,0.00125,0.1];
len_y_arr=[1,1,1,1,1,1,1];
cond_input=[2e-1,3.4e-1,3.4e-1,3.4e-1,3.4e-1,7e-1,7e-1];
permit_input=[5e2,2.6e4,2.6e4,2.6e4,2.6e4,5.2e3,5.2e3];
art_en=[0 0 0 0 0 1 1];
cond_input_rel=cond_input./4.3e-2;
r_ratio_arr=1./cond_input_rel;  %R=L/cond*A
ri_dc=10;
re_dc=10;

permit_input_rel=permit_input./9.1e2;
c_ratio_arr=permit_input_rel;   %c=eps*A/d
cm_dc=100e-9;

% Electrode configuration
Nelec_arr=20;
EIT_spacing_arr=4;
direction_arr=0; %0=> Adjacent, 1=>Opposite

% Simulation Paramters (Frequency, Time, R vs RC)
ac_freq=[10 1e3 100e3]; %npts fstart fend
z_en=1;

name_base={'Inhomo Sens. '};
net.param.EIT_spacing=EIT_spacing_arr;
net.param.Nelec=Nelec_arr;
net.param.dia=0;
net.param.direction=direction_arr;
net.param.ri_dc=ri_dc;
net.param.re_dc=re_dc;
net.param.cm_dc=cm_dc;
net.param.ac_freq=ac_freq;
net.param.z_en=z_en;

for i=1:length(center_x_arr)
    net.param.obj(i).r_ratio=r_ratio_arr(i);
    net.param.obj(i).c_ratio=c_ratio_arr(i);
    net.param.obj(i).center.x=center_x_arr(i);
    net.param.obj(i).center.y=center_y_arr(i);
    net.param.obj(i).center.z=center_z_arr(i);
    net.param.obj(i).len.x=len_x_arr(i);
    net.param.obj(i).len.y=len_y_arr(i);
    net.param.obj(i).len.z=len_z_arr(i);
    net.param.obj(i).art_en=art_en(i);
end


% Reference
%Create electrodes
[net]=eit_elec_array_fn(net);
%Create Grid
net=img_gen_fn(net);
net=grid_gen_fn(net);
% Simulate SPICE Netlist
net_ref{1}=run_sim_fn(net);
plot_single_img_fn(net_ref{1}(5),0,0);
plot_single_img3d_fn(net_ref{1}(2),1,[]);


if net.param.plot_eit==1 && net.param.debug==0
    plot_spice_arr_fn(net_ref,[1:1],[1:1:10],2)
end


if net.param.write_figures==1 && net.param.debug==0
    if ~exist(net.param.path,'dir')
        mkdir(net.param.path);
    end
    file_str=strcat(net.param.path);
    print(strcat(file_str,'img_ref.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
    savefig(strcat(file_str,'img_ref.fig')); %<-Save as PNG with 300 DP
end
% 
% for n=1:length(n_arr)
%     for m=1:length(m_arr)
%         %Create electrodes
%         %[net]=eit_elec_array_fn(net);
%         %Create Grid
%         %net=imp_grid_gen_fn(net);
%         % Update Grid for Senstivity Analysis
%         %Create Grid
%         i=n_arr(n);
%         j=m_arr(m);
%         net=img_gen_fn(net);
%         net.param.update.center_x=(i);
%         net.param.update.center_y=(j);
%         net.param.update.delta_r_dc=0.1;  % percentage
%         net=img_update_fn(net);
%         net=grid_gen_fn(net);
%         % Simulate SPICE Netlist
%         net_arr{n,m}=EIT_SPICE_fn(net);
%     end
% end
% 
% %name_arr=strcat(repmat(name_base,size(net_arr,1)*size(net_arr,2),1),num2str(reshape([round(center_y_arr*(net_arr{1,1}(1).Ni+1))],1,[])'));
% plot_img_fn(net_arr,1,0);
% 
% 
% if net.param.write_figures==1 && net.param.debug==0
%     if ~exist(net.param.path,'dir')
%         mkdir(net.param.path);
%     end
%     file_str=strcat(net.param.path);
%     print(strcat(file_str,'img_step.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
%     savefig(strcat(file_str,'img_step.fig')); %<-Save as PNG with 300 DP
% end

if net.param.write_mat==1 && net.param.debug==0
    if ~exist(net.param.path,'dir')
        mkdir(net.param.path);
    end
    save(strcat(net.param.path,'netlist.mat'));
end


%%

if net.param.debug==0
    
    % single senestivity map per electrode config
    curr_idx_arr=1:length(net_arr{1,1});
    s_idx_arr=1:length(net_arr{1,1}(1).s);
    for k=1:length(curr_idx_arr)
        for s=1:length(s_idx_arr)
            curr_idx=curr_idx_arr(k);
            s_idx=s_idx_arr(s);
            n=size(net_arr,1);
            m=size(net_arr,2);
            for i=1:n
                for j=1:m
                    s_arr=cat(1,net_arr{i,j}(curr_idx).s);
                    dVe_arr_all(k,s,j,i)=s_arr(s_idx).Vs-net_ref{1}(curr_idx).s(s_idx).Vs;
                end
            end
        end
    end
    
    % Plotting
    plot_curr_idx_arr=[1,5];
    plot_s_idx_arr=[2:4:15];
    figure
    a=ceil(sqrt(length(plot_curr_idx_arr)*length(plot_s_idx_arr)));
    b=ceil(length(plot_curr_idx_arr)*length(plot_s_idx_arr)/a);
    clims_max=max(abs(max(dVe_arr_all(plot_curr_idx_arr,plot_s_idx_arr,:,:),[],'all')),abs(min(dVe_arr_all(plot_curr_idx_arr,plot_s_idx_arr,:,:),[],'all')));
    clims=[-clims_max,clims_max]/2;
    p=0;
    for kp=1:length(plot_curr_idx_arr)
        for sp=1:length(plot_s_idx_arr)
            k=plot_curr_idx_arr(kp);
            s=plot_s_idx_arr(sp);
            curr_idx=curr_idx_arr(k);
            s_idx=s_idx_arr(s);
            p=p+1;
            subplot(a,b,p)
            grid on;hold on;
            dVe_arr=squeeze(dVe_arr_all(k,s,:,:));
            dVe_arr(sub2ind(size(dVe_arr),[net_ref{1}(curr_idx).s(s_idx).n1i,net_ref{1}(curr_idx).s(s_idx).n2i],[net_ref{1}(curr_idx).s(s_idx).n1j,net_ref{1}(curr_idx).s(s_idx).n2j]))=inf;
            dVe_arr(sub2ind(size(dVe_arr),[net_ref{1}(curr_idx).curr.n1i,net_ref{1}(curr_idx).curr.n2i],[net_ref{1}(curr_idx).curr.n1j,net_ref{1}(curr_idx).curr.n2j]))=NaN;
            mesh_plot_fn(dVe_arr,1,clims);
            title(strcat('Impedance Senestivity Map k=',num2str(curr_idx),', s=',num2str(s_idx),' (mV)'));
            xlabel('X');ylabel('Y');
        end
    end
    if net.param.write_figures==1 && net.param.debug==0
        if ~exist(net.param.path,'dir')
            mkdir(net.param.path);
        end
        file_str=strcat(net.param.path);
        print(strcat(file_str,'eit_sens_single.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
        savefig(strcat(file_str,'eit_sens_single.fig')); %<-Save as PNG with 300 DP
    end
    
    % Combined senestivity map
    dVe_arr_all_re=reshape(dVe_arr_all,[],n,m);
    dVenorm_arr_all=squeeze(vecnorm(dVe_arr_all_re,2,1));
    figure
    grid on;hold on;
    clims=[min(dVenorm_arr_all,[],'all'),max(dVenorm_arr_all,[],'all')];
    mesh_plot_fn(dVenorm_arr_all',1,clims/2.5);
    xlim([1.5,size(dVenorm_arr_all,1)-0.5]);
    ylim([1.5,size(dVenorm_arr_all,2)-0.5]);
    title(strcat('Combined Impedance Sensitivity Map (mV)'));
    xlabel('X');ylabel('Y');
    ax = gca; ax.FontSize=13;
    
    if net.param.write_figures==1 && net.param.debug==0
        if ~exist(net.param.path,'dir')
            mkdir(net.param.path);
        end
        file_str=strcat(net.param.path);
        print(strcat(file_str,'eit_sens_comb.png'), '-dpng', '-r600'); %<-Save as PNG with 300 DP
        savefig(strcat(file_str,'eit_sens_comb.fig')); %<-Save as PNG with 300 DP
    end
end