function net=spice_netlist_3d_pp_fn(net)

%   Bassem Ibrahim
%   2019.08.22
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation
%   Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling,
%   IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan.
%
%   function net=spice_netlist_3d_pp_fn(net)
%   Post-processing for voltage and PTT calculations
%   net.V_node: Node voltage
%   net.V_node_NaN: Node voltage with NaN at electrode nodes
%   net.V_sns_v: Vertical sensing voltage (Voltage difference between each 2 adjacent vertical nodes)
%   net.V_sns_h: Horizontal sensing voltage (Voltage difference between each 2 adjacent horizontal nodes)
%   net.V_sns_h_NaN: Horizontal sensing voltage with NaN at electrode nodes
%   net.V_sns_v_NaN: Vertical sensing voltage with NaN at electrode nodes
%   net.V_node_diff: Node delta voltage
%   net.V_node_diff_NaN: Node delta voltage
%   net.V_sns_v_diff: Vertical sensing delta voltage
%   net.V_sns_h_diff: Horizontal sensing delta voltage
%   net.PAT_r_var: Artery's pulse arrival time (PAT)
%   net.PAT_sns_v_space_n_k1: Sensed pulse arrival time (PAT) for different spacing between electrodes and depths
%   net.PAT_sns_v_n_k: Sensed pulse arrival time (PAT) for different depths
%   net.PTT_sns_v_space: Sensed pulse arrival time (PAT) for different spacing between electrodes
%   net.PTT_sns_v_n_k1: Sensed pulse arrival time (PAT) at skin surface (Z=0)
%   net.V_s: DC Voltage at sensing electrodes 
%   net.V_s_diff: Delta voltage at sensing electrodes
%   net.PAT_V_s: Pulse arrival time (PAT) at sensing electrodes

V_node=net.V_node;
t_var=net.t_var;
L=length(t_var);
node_arr=net.node_arr;
S=net.param.S;
V_node_NaN=V_node;

for i=1:length(net.curr)
    V_node_NaN(net.curr(i).n1i,net.curr(i).n1j,1,:)=NaN;
    V_node_NaN(net.curr(i).n2i,net.curr(i).n2j,1,:)=NaN;
end
if net.param.r_var_max>0 && ~isempty(strmatch(net.param.write_netlist,'pene'))
    V_node_NaN(net.r_var.n1i,net.r_var.n1j,net.r_var.n1k,:)=NaN;
    V_node_NaN(net.r_var.n2i,net.r_var.n2j,net.r_var.n2k,:)=NaN;
end
if ~isempty(net.param.s_n1i)
    for i=1:length(net.s)
        V_node_NaN(net.s(i).n1i,net.s(i).n1j,1,:)=NaN;
        V_node_NaN(net.s(i).n2i,net.s(i).n2j,1,:)=NaN;
    end
end
spacing_single=1;
V_sns_v_NaN=zeros(size(V_node_NaN)-[1 zeros(1,length(size(V_node_NaN))-1)]);
V_sns_v=zeros(size(V_node)-[1 zeros(1,length(size(V_node))-1)]);
for k=1:size(V_node,1)-spacing_single
    V_sns_v_NaN(k,:,:,:)=V_node_NaN(k,:,:,:)-V_node_NaN(k+spacing_single,:,:,:);
    V_sns_v(k,:,:,:)=V_node(k,:,:,:)-V_node(k+spacing_single,:,:,:);
end
V_sns_h_NaN=zeros(size(V_node_NaN)-[0 1 zeros(1,length(size(V_node_NaN))-2)]);
V_sns_h=zeros(size(V_node_NaN)-[0 1 zeros(1,length(size(V_node_NaN))-2)]);
for k=1:size(V_node,2)-spacing_single
    V_sns_h_NaN(:,k,:,:)=squeeze(V_node_NaN(:,k,:,:)-V_node_NaN(:,k+spacing_single,:,:));
    V_sns_h(:,k,:,:)=squeeze(V_node(:,k,:,:)-V_node(:,k+spacing_single,:,:));
end

e=net.param.EIT_spacing;
if length(size(V_node_NaN))==2
    V_node_NaN(:,:,1)=V_node_NaN;
end
if length(size(V_node_NaN))<3
    net.V_node_EIT=[squeeze(V_node_NaN(1,round(e/2)+1:e:end-round(e/2),:)),squeeze(permute(V_node_NaN(round(e/2)+1:e:end-round(e/2),end,:),[2 1 3 4])),squeeze(V_node_NaN(end,end-round(e/2):-e:round(e/2)+1,:)),squeeze(permute(V_node_NaN(end-round(e/2):-e:round(e/2)+1,1,:),[2 1 3 4]))];
else
    net.V_node_EIT=[squeeze(V_node_NaN(1,round(e/2)+1:e:end-round(e/2),:));squeeze(permute(V_node_NaN(round(e/2)+1:e:end-round(e/2),end,:),[2 1 3 4]));squeeze(V_node_NaN(end,end-round(e/2):-e:round(e/2)+1,:));squeeze(permute(V_node_NaN(end-round(e/2):-e:round(e/2)+1,1,:),[2 1 3 4]))];
end

net.V_EIT=circshift(net.V_node_EIT,-1)-net.V_node_EIT;
if net.param.mode==0
    net.V_EIT_mean=(net.V_EIT(:,round(L/2)));
    net.V_EIT_max=(net.V_EIT(:,round(L/4)));
end
if net.param.mode==0
    V_node_diff=max(V_node,[],4)-min(V_node,[],4);
    V_node_diff_NaN=max(V_node_NaN,[],4)-min(V_node_NaN,[],4);
    V_sns_v_diff=V_sns_v(:,:,:,1)-V_sns_v(:,:,:,dsearchn(net.t_var,1/net.param.r_var_freq/2));
    V_sns_h_diff=V_sns_h(:,:,:,1)-V_sns_h(:,:,:,dsearchn(net.t_var,1/net.param.r_var_freq/2));
end

max_spacing=(net.param.i_src_n2i-net.param.i_src_n1i-2-net.curr.dia_y);
if net.param.sim_speed==1
    spacing_arr=1:max_spacing;
    dist_step=[1 1 1];
elseif net.param.sim_speed==2
    spacing_arr=round(linspace(1,max_spacing,4));
    dist_step=[ceil(net.Ni/10) ceil(net.Nj/10) ceil(net.Nk/10)];
elseif net.param.sim_speed==3
    if ((net.param.i_src_n2i-net.param.i_src_n1i-2)>0)&&((net.param.i_src_n2j-net.param.i_src_n1j-2)>0) 
        spacing_arr=1;
    else
        spacing_arr=[];
    end
    dist_step=[ceil(net.Ni/3) ceil(net.Nj/3) ceil(net.Nk/3)];
elseif net.param.sim_speed>3
    spacing_arr=net.param.spacing_arr;
    dist_step=[1 1 ceil(net.Nk/net.param.dist_step_count)];
end

V_sns_v_space=[];
V_sns_v_space_NaN=[];
for p=1:length(spacing_arr)
    spacing=spacing_arr(p);
    for k=1:size(V_node,1)-spacing
        V_sns_v_space_NaN{p}(k,:,:,:)=squeeze(V_node_NaN(k,:,:,:)-V_node_NaN(k+spacing,:,:,:));
        V_sns_v_space{p}(k,:,:,:)=squeeze(V_node(k,:,:,:)-V_node(k+spacing,:,:,:));
    end
end


% PTT

if net.param.r_var_max>0 && net.param.mode==0
    for s=1:S
        r_var_t(s,:,:)=-net.param.r_var_max(s)*sin(2*pi*net.param.r_var_freq(s)*(t_var+net.r_var.delayi(s,:))-pi/2);
    end
    
    
    
    if S==1 && net.param.r_var_dim==3 && net.param.mode==0
        for i=1:size(r_var_t,3)
            [PAT_r_var(i),amp_r_var(i),dc_r_var(i),r_var_fiterr(i)]=find_delay_fn(t_var,squeeze(r_var_t(s,:,i))',net.param.r_var_freq,0);
        end
        
        
        % PAT in X,Y,spacing at Z=0
        for p=1:length(V_sns_v_space)
            PAT_sns_v_space_n_k1{p}=NaN(net.Ni+1-spacing_arr(p),net.Nj+1);
            for i=1:dist_step(1):size(V_sns_v_space{p},1)
                for j=1:dist_step(2):size(V_sns_v_space{p},2)
                    [PAT_sns_v_space_n_k1{p}(i,j),amp_sns_v_space_n_k1{p}(i,j),dc_sns_v_space_n_k1{p}(i,j),V_sns_v_space_n_k1_fiterr{p}(i,j)]=find_delay_fn(t_var,squeeze(V_sns_v_space{p}(i,j,1,:)),net.param.r_var_freq(1),0);
                end
            end
        end
        PTT_sns_v_space=[];
        for p=1:length(V_sns_v_space)
            PTT_sns_v_space(p,:)=PAT_sns_v_space_n_k1{p}(max(net.curr.n2i)-1-spacing_arr(p),:)-PAT_sns_v_space_n_k1{p}(min(net.curr.n1i)+1,:);
        end
        
        % PAT in X,Y,Z at spacing=1
        if ~isempty(V_sns_v_space)
            p=1;
            PAT_sns_v_n_k=NaN(net.Ni,net.Nj+1,net.Nk+1);
            amp_sns_v_n_k=NaN(net.Ni,net.Nj+1,net.Nk+1);
            dc_sns_v_n_k=NaN(net.Ni,net.Nj+1,net.Nk+1);
            V_sns_v_n_k_fiterr=NaN(net.Ni,net.Nj+1,net.Nk+1);
            for i=1:dist_step(1):size(V_sns_v_space{p},1)
                for j=1:dist_step(2):size(V_sns_v_space{p},2)
                    for k=1:dist_step(3):size(V_sns_v_space{p},3)
                        [PAT_sns_v_n_k(i,j,k),amp_sns_v_n_k(i,j,k),dc_sns_v_n_k(i,j,k),V_sns_v_n_k_fiterr(i,j,k)]=find_delay_fn(t_var,squeeze(V_sns_v_space{p}(i,j,k,:)),net.param.r_var_freq(1),0);
                    end
                end
            end
            
            PTT_sns_v_n_k1=PAT_sns_v_n_k(end,:,1)-PAT_sns_v_n_k(1,:,1);
            PTT_r_var=PAT_r_var(end)-PAT_r_var(1);
            
            % Normalize r var for plotting
            V_sns_v_c=(V_sns_v-repmat(dc_sns_v_n_k,1,1,1,size(V_sns_v,4)));
            V_sns_v_n=-1*V_sns_v_c./repmat(amp_sns_v_n_k,1,1,1,size(V_sns_v,4));
        end
    end
end


net.V_node=V_node;
net.V_node_NaN=V_node_NaN;
net.V_sns_v=V_sns_v;
net.V_sns_h=V_sns_h;
net.V_sns_h_NaN=V_sns_h_NaN;
net.V_sns_v_NaN=V_sns_v_NaN;


net.V_sns_v_space=V_sns_v_space;
net.spacing_arr=spacing_arr;
net.spacing_single=spacing_single;
if S==1 && net.param.r_var_max>0    && net.param.mode==0
    
    net.V_node_diff=V_node_diff;
    net.V_node_diff_NaN=V_node_diff_NaN;
    net.V_sns_v_diff=V_sns_v_diff;
    net.V_sns_h_diff=V_sns_h_diff;
    net.r_var.r_var_t=r_var_t;
    if net.param.r_var_dim==3
        if ~isempty(V_sns_v_space)
            net.PAT_sns_v_space_n_k1=PAT_sns_v_space_n_k1;
            net.PAT_sns_v_n_k=PAT_sns_v_n_k;
            net.PTT_sns_v_space=PTT_sns_v_space;
            net.PTT_sns_v_n_k1=PTT_sns_v_n_k1;
            net.V_sns_v_n=V_sns_v_n;
            net.PAT_r_var=PAT_r_var;
            net.PTT_r_var=PTT_r_var;
        end        
    end
end
net.L=L;
if ~isempty(net.param.s_n1i)
    for i=1:length(net.s)
        net.V_s(i,:)=squeeze((net.V_node(net.s(i).n1i(1),net.s(i).n1j(1),1,:)))-squeeze(net.V_node(net.s(i).n2i(1),net.s(i).n2j(1),1,:));
    end
    net.V_s_diff=max(net.V_s,[],2)-min(net.V_s,[],2);
    if net.param.mode==0
        for i=1:length(net.s)
            %net.V_s(i,:)=squeeze((net.V_node(net.s(i).n1i(1),net.s(i).n1j(1),1,:)))-squeeze(net.V_node(net.s(i).n2i(1),net.s(i).n2j(1),1,:));
            [PAT_V_s(i),dVpp_s(i),Vdc_s(i)]=find_delay_fn(t_var,net.V_s(i,:),net.param.r_var_freq(1),0);
        end
        net.PAT_V_s=PAT_V_s;
        net.dVpp_s=dVpp_s;
        net.Vdc_s=Vdc_s;
    end
end

