function net_arr_t1=run_sim_fn(net)

load(net.param.skin_elec_data_path)
Ni=net.Ni;
Nj=net.Nj;
Nk=net.Nk;

i_src_dist=0;
fi_req=[1];
i_num=1;
i_opp=0;
i_area=1;
ival=0.5e-3;
small_edges=net.param.small_edges;
i_dia_x=net.param.i_dia_x;
s_dia_x=net.param.s_dia_x;
i_dia_y=net.param.i_dia_y;
s_dia_y=net.param.s_dia_y;
nstep=8;
tindex=0;
sim_speed=3;

% Configurations
z_en=net.param.z_en;
se_en=net.param.se_en;
i_freq=100;
% mode
% 0=> Transient Analysis      
% 1=> Single DC Bio-Z        re_var_max=0,t=0
% 2=> Multi DC Bio-Z (Pulse) re_var_max>0,t=time_arr
if net.param.re_var_max>0
    mode=2;  
else
    mode=1;
end
ac_freq=net.param.ac_freq;
L=1;
reltol=0.001;
abstol=1e-6;
gmin=0;
volttol=1e-5;
a_se=net.param.se_ratio;
r1_se=data.R1*a_se; %1.12e5
c1_se=data.C1/a_se; %1.97e-9
r2_se=data.R2*a_se; %4.8e3
c2_se=data.C2/a_se; %1.97e-9
% r1_se=0;
% r2_se=1;
% c_se=100e-9;

%time_en=1; %1=> time in R,C var
if mode==1
    time_arr=0;
    iac=ival;
    idc=0;
else
    iac=ival;
    idc=0;
    if net.param.debug==0
        r_var_time_arr=1./min(net.param.r_var_freq,[],2)*1;
        time_arr=[0:r_var_time_arr/nstep:r_var_time_arr/nstep*(nstep-1)]';
    else
        time_arr=0;
    end
end
net.param.color_arr=[0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    1.0000         0         0;
    0.2500    0.2500    0.2500;
    0    0.5000         0];

for t1=1:length(net.elec_array.i_src_n1i_arr)
    clear net_tarr;
    if exist('net.param.r_var_freq')
        tstep=1./min([net.param.r_var_freq],[],2)*2*nstep;
    else
        tstep=1;
    end
    
    for t0=1:length(time_arr)
        %% Generate and Simualte Spice        
        
        net.param.i_src_n1i=round(net.elec_array.i_src_n1i_arr(t1,:)/L);
        net.param.i_src_n2i=round(net.elec_array.i_src_n2i_arr(t1,:)/L);
        net.param.i_src_n1j=round(net.elec_array.i_src_n1j_arr(t1,:)/L);
        net.param.i_src_n2j=round(net.elec_array.i_src_n2j_arr(t1,:)/L);
        net.param.i_src_n1k=round(net.elec_array.i_src_n1k_arr(t1,:)/L);
        net.param.i_src_n2k=round(net.elec_array.i_src_n2k_arr(t1,:)/L);
        net.param.i_src_dist=i_src_dist;
        net.param.i_num=i_num;
        net.param.i_opp=i_opp;
        net.param.tstep=tstep;
        % Current Source 1
        net.curr(1).dia_x=round(i_dia_x/L);
        net.curr(1).dia_y=round(i_dia_y/L);
        net.curr(1).iac=iac;
        net.curr(1).idc=idc;
        
        net.curr(1).freq=i_freq;
        net.param.s_n1i=round(net.elec_array.s_n1i_arr(t1,:)/L);
        net.param.s_n2i=round(net.elec_array.s_n2i_arr(t1,:)/L);
        net.param.s_n1j=round(net.elec_array.s_n1j_arr(t1,:)/L);
        net.param.s_n2j=round(net.elec_array.s_n2j_arr(t1,:)/L);
        net.param.s_n1k=round(net.elec_array.s_n1k_arr(t1,:)/L);
        net.param.s_n2k=round(net.elec_array.s_n2k_arr(t1,:)/L);
        net.param.s_dia_x=round(s_dia_x/L);
        net.param.s_dia_y=round(s_dia_y/L);        
        
        net.param.fi_req=fi_req;
        net.param.sim_speed=sim_speed;
        net.param.mode=mode;
        net.param.z_en=z_en;
        net.param.time=time_arr(t0);
        net.param.ac_freq=ac_freq;
        net.param.r1_se=r1_se;
        net.param.r2_se=r2_se;
        net.param.c1_se=c1_se;
        net.param.c2_se=c2_se;
        net.param.se_en=se_en;
        net.param.reltol=reltol;
        net.param.abstol=abstol;
        net.param.gmin=gmin;
        net.param.volttol=volttol;
        
        tindex=tindex+1;
        fprintf('Time Step %2d, curr n1=(%3d,%3d), n2=(%3d,%3d)',tindex,net.elec_array.i_src_n1i_arr(t1,:),net.elec_array.i_src_n1j_arr(t1,:),net.elec_array.i_src_n2i_arr(t1,:),net.elec_array.i_src_n2j_arr(t1,:));
        if exist('net.param.r_var_freq')
            fprintf(', artery nk=%d',net.param.r_var_n1k);
        end
        
        % Generate the impedance(resistance) image
        
        % Run
        net=spice_netlist_3d_run_fn(net);
        if t0==1
            net_tarr(1)=net;
        else
            net_tarr(1).param=net.param;
        end
        % Sparse Spice Output
        if net.param.debug==0
            net_tarr=read_spice_out_fast_gnd_fn(net_tarr);
        end
        if t0==1
            for i=1:length(net_tarr)
                for j=1:length(net_tarr(i).curr)
                    net_tarr(i).V_node=[];
                    net_tarr(i).curr(j).Ve_n1=[];
                    net_tarr(i).curr(j).Ve_n2=[];
                    net_tarr(i).curr(j).Ve_a=[];
                    net_tarr(i).curr(j).Ve=[];
                end
                for j=1:length(net_tarr(i).s)
                    net_tarr(i).s(j).Ve_n1=[];
                    net_tarr(i).s(j).Ve_n2=[];
                    net_tarr(i).s(j).Ve_a=[];
                    net_tarr(i).s(j).Ve=[];
                end
            end
        end
        net_tarr(1).param.index=t1;
        if net.param.debug==0
            for i=1:length(net_tarr)
                net_tarr(i).V_node=cat(4,net_tarr(i).V_node,net_tarr(i).V_node_x);
                for j=1:length(net_tarr(i).curr)
                    net_tarr(i).curr(j).Ve_n1=cat(2,net_tarr(i).curr(j).Ve_n1,net_tarr(i).curr(j).Ve_n1_x);
                    net_tarr(i).curr(j).Ve_n2=cat(2,net_tarr(i).curr(j).Ve_n2,net_tarr(i).curr(j).Ve_n2_x);
                    net_tarr(i).curr(j).Ve_a=cat(2,net_tarr(i).curr(j).Ve_a,net_tarr(i).curr(j).Ve_a_x);
                    net_tarr(i).curr(j).Ve=cat(2,net_tarr(i).curr(j).Ve,net_tarr(i).curr(j).Ve_x);
                end
                for j=1:length(net_tarr(i).s)
                    net_tarr(i).s(j).Ve_n1=cat(2,net_tarr(i).s(j).Ve_n1,net_tarr(i).s(j).Ve_n1_x);
                    net_tarr(i).s(j).Ve_n2=cat(2,net_tarr(i).s(j).Ve_n2,net_tarr(i).s(j).Ve_n2_x);
                    net_tarr(i).s(j).Ve_a=cat(2,net_tarr(i).s(j).Ve_a,net_tarr(i).s(j).Ve_a_x);
                    net_tarr(i).s(j).Ve=cat(2,net_tarr(i).s(j).Ve,net_tarr(i).s(j).Ve_x);
                end
            end
        end
        %clear net
    end
    if net.param.debug==0
        %if net.param.mode==1 && length(time_arr)>6
        %r_var_arr=cat(1,[net_tarr.r_var]);
        %r_art_all_arr_arr=[r_var_arr.r_art_all_arr];
        %         V_node_arr=cat(5,net_tarr.V_node);
        %         V_node_arr=permute(V_node_arr,[1 2 3 5 4]);
        %
        %         V_node_arr=cat(5,net_tarr.V_node);
        %         V_node_arr=permute(V_node_arr,[1 2 3 5 4]);
        
        
        %net_t_r.param.mode=0;
        fprintf('Post processing')
        net_temp=net_tarr(1);
        for i=1:length(net_tarr)
            fprintf('.') 
            % Post Processing - Pulse Bio-Z
            net_temp.V_node=net_tarr(i).V_node;
            net_temp.I_node=net_tarr(i).I_node;
            for j=1:length(net_tarr(i).curr)
                net_temp.curr(j).Ve_n1=net_tarr(i).curr(j).Ve_n1;
                net_temp.curr(j).Ve_n2=net_tarr(i).curr(j).Ve_n2;
                net_temp.curr(j).Ve_a=net_tarr(i).curr(j).Ve_a;
                net_temp.curr(j).Ve=net_tarr(i).curr(j).Ve;
            end
            for j=1:length(net_tarr(i).s)
                net_temp.s(j).Ve_n1=net_tarr(i).s(j).Ve_n1;
                net_temp.s(j).Ve_n2=net_tarr(i).s(j).Ve_n2;
                net_temp.s(j).Ve_a=net_tarr(i).s(j).Ve_a;
                net_temp.s(j).Ve=net_tarr(i).s(j).Ve;
            end
            net_temp.t_var=time_arr;
            net_temp=spice_netlist_3d_pp_fn(net_temp); 
            net_arr(i)=net_temp;
        end
        fprintf('\n')    
        net_arr_t1(t1,:)=net_arr;
    else
        net_arr_t1(t1,:)=net_tarr;
    end
end


