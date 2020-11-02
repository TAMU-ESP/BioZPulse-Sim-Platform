function net=spice_netlist_3d_gen_fn(net)

%   Bassem Ibrahim
%   2020.11.02
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari
%
%   function net=spice_netlist_3d_gen_fn(net)
%   Generate SPICE netlist file (netlsit.cir)
%
%   Model Paramters
% 1     net.Ni					Number of nodes in Y-axis (YL)
% 2 	net.Nj					Number of nodes in X-axis (XL)
% 3     net.Nk					Number of nodes in Z-axis (ZL)
% 4     net.param.r_dc			Tissue's extra-cellular resistance (RE)
% 5     net.param.ri_dc			Tissue's intra-cellular resistance (RI)
% 6     net.param.ci_dc			Tissue's cell membrane capacitance (Cm)
% 7     net.param.r_var_num		Number of arteries
% 8     net.param.re_art_dc		Artery's extra-cellular resistance (RE)
% 9     net.param.ri_art_dc		Artery's intra-cellular resistance (RI)
% 10	net.param.cm_art_dc		Artery's cell membrane capacitance (Cm)
% 11	net.param.re_var_max		Artery's variable extra-cellular resistance (deltaRE)
% 12	net.param.ri_var_max	Artery's variable intra-cellular resistance (deltaRI)
% 13	net.param.cm_var_max	Artery's variable cell membrane capacitance (deltaCm)
% 14	net.param.r_var_steps	Artery's number of steps
% 15	net.param.r_var_n1i		Artery's location on Y-axis (YA)
% 16	net.param.r_var_n1j		Artery's location on X-axis (XA)
% 17	net.param.r_var_n1k		Artery's location on Z-axis (ZA)
% 18	net.param.r_var_v		Artery's direction: 1-> vertical, 0-> horizontal
% 19	net.param.r_var_dim		Artery's dimension (1D, 2D or 3D)
% 20	net.param.r_var_seg		Artery's number of voxels
% 21	net.param.r_var_dia		Artery's diameter (DA)
% 22	net.param.r_var_delay	Artery's pulse transit time (PTT)
% 23	net.param.r_var_freq	Artery's pulse rate (f)
% 24	net.param.s_n1i			Location on Y-axis of the first voltage electrode
% 25	net.param.s_n2i			Location on Y-axis of the second voltage electrode
% 26	net.param.s_n1j			Location on X-axis of the first voltage electrode
% 27	net.param.s_n2j			Location on X-axis of the second voltage electrode
% 28	net.param.s_dia_x		Size of voltage electrode in X-axis
% 29	net.param.s_dia_y		Size of voltage electrode in Y-axis
% 32	net.param.i_num			Number of current sources
% 33	net.param.i_src_n1i		Location on Y-axis of the first current electrode
% 34	net.param.i_src_n2i		Location on Y-axis of the second current electrode
% 35	net.param.i_src_n1j		Location on X-axis of the first current electrode
% 36	net.param.i_src_n2j		Location on X-axis of the second current electrode
% 37	net.param.i_opp			Direction of current
% 38	net.curr.ival			DC current amplitude
% 39	net.curr.ival			AC current amplitude
% 40	net.curr.freq			AC current frequency
% 30	net.curr.dia_x			Size of current electrode in X-axis
% 31	net.curr.dia_y			Size of current electrode in Y-axis
% 41	net.param.r_var_time	Simulation time
% 42	net.param.tstep			Simulation time step

tic
file_name=net.param.circuit_file;
if exist(file_name,'file')
    delete(file_name);
end


Ni=net.Ni;
Nj=net.Nj;
Nk=net.Nk;

tstep=net.param.tstep;

% Create Node Array
for i=1:Ni+1
    for j=1:Nj+1
        for k=1:Nk+1
            node_arr(i,j,k)=xyz2node_fn(i,j,k,Ni,Nj,Nk);
        end
    end
end
net.node_arr_raw=node_arr;
%img_arr=net.re_img;
img_arr_nonan=abs(1./((1./net.re_img)+1./(net.ri_img+1./(1i*(2*pi*10e3)*net.cm_img)))).^(1/1);
img_arr=img_arr_nonan;

% Variable Resistor nodes
if net.param.re_var_max(1)>0
    freq=net.param.r_var_freq;
    tend=net.param.r_var_time;
    for s=1:length(net.param.r_var_n1j)
        td=net.param.r_var_delay(s);
        delay_line=linspace(0,td,net.param.r_var_seg(s));
        S=floor((net.param.r_var_num(s))/net.param.r_var_seg(s));
        net.r_var(s).delayi=reshape(repmat(delay_line,S,1),1,[]);
        net.r_var(s).delayi=[net.r_var(s).delayi repmat(net.r_var(s).delayi(end),1,net.Ni-S*net.param.r_var_seg(s))];
        if net.param.r_var_dim(s)==1
            net.r_var.n1i=net.param.r_var_n1i;
            net.r_var.n1j=net.param.r_var_n1j;
            if net.param.r_var_v==1
                net.r_var.n2i=net.param.r_var_n1i+1;
                net.r_var.n2j=net.param.r_var_n1j;
            else
                net.r_var.n2i=net.param.r_var_n1i;
                net.r_var.n2j=net.param.r_var_n1j+1;
            end
            net.r_var.n1k=net.param.r_var_n1k;
            net.r_var.n2k=net.param.r_var_n1k;
        elseif net.param.r_var_dim(s)==2
            %             n1i=net.param.r_var_n1i+1+[0 0 1 0];
            %             n1j=net.param.r_var_n1j+1+[0 1 0 0];
            %             n2i=net.param.r_var_n1i+1+[0 1 1 1];
            %             n2j=net.param.r_var_n1j+1+[1 1 1 0];
            %             n1k=net.param.r_var_n1k+1+[0 0 0 0];
            %             n2k=net.param.r_var_n1k+1+[0 0 0 0];
            %             delay=[0,0,0,0];
            n1i=[];n1j=[];n1k=[];n2i=[];n2j=[];n2k=[];delay=[];
            for j=1:net.param.r_var_dia(s)
                n1ibase=net.param.r_var_n1i(s)+1;
                n1jbase=net.param.r_var_n1j(s)+1+j-1;
                n1kbase=1;
                
                n1i=[n1i n1ibase+[0]];
                n1j=[n1j n1jbase+[0]];
                n1k=[n1k n1kbase+[0]];
                
                n2i=[n2i n1ibase+[0]];
                n2j=[n2j n1jbase+[1]];
                n2k=[n2k n1kbase+[0]];
                delay=[delay repmat(net.r_var.delayi(s,1),1,4)];
                for i=1:net.param.r_var_num(s)
                    n1i=[n1i n1ibase+[0 0 1]];
                    n1j=[n1j n1jbase+[0 1 0]];
                    n1k=[n1k n1kbase+[0 0 0]];
                    
                    n2i=[n2i n1ibase+[1 1 1]];
                    n2j=[n2j n1jbase+[0 1 1]];
                    n2k=[n2k n1kbase+[0 0 0]];
                    
                    delay=[delay repmat(net.r_var.delayi(s,i),1,4)];
                    
                    n1ibase=n1ibase+1;
                end
            end
        elseif net.param.r_var_dim(s)==3
            n1i=[];n1j=[];n1k=[];n2i=[];n2j=[];n2k=[];delay=[];n1ibase_arr=[];n1jbase_arr=[];n1kbase_arr=[];
            for j=1:net.param.r_var_dia(s)
                for k=1:net.param.r_var_dia(s)
                    n1ibase=net.param.r_var_n1i(s)+1;
                    n1jbase=net.param.r_var_n1j(s)+1+j-1;
                    n1kbase=net.param.r_var_n1k(s)+1+k-1;                   
                                        
                    n1i=[n1i n1ibase+[0 0 0 0]];
                    n1j=[n1j n1jbase+[0 1 0 0]];
                    n1k=[n1k n1kbase+[0 0 1 0]];
                    
                    n2i=[n2i n1ibase+[0 0 0 0]];
                    n2j=[n2j n1jbase+[1 1 1 0]];
                    n2k=[n2k n1kbase+[0 1 1 1]];
                    
                    n1ibase_arr=[n1ibase_arr n1ibase];
                    n1jbase_arr=[n1jbase_arr n1jbase];
                    n1kbase_arr=[n1kbase_arr n1kbase];
                    
                    delay=[delay repmat(net.r_var(s).delayi(1),1,4)];
                    for i=1:net.param.r_var_num(s)
                        if net.param.r_var_tilt(s)>0 && i==round((net.param.r_var_num(s)+1)*(1/2))
                            n1kbase=n1kbase-net.param.r_var_tilt(s); % upper part (i>Ni/2) get closer to surface
                        %elseif net.param.r_var_tilt(s)>0 && i==round((net.param.r_var_num(s)+1)*(2/3))                
                        %    n1kbase=n1kbase-net.param.r_var_tilt(s); % upper part (i>Ni/2) get closer to surface
                        end
                        n1i=[n1i n1ibase+[0 0 0 0 1 1 1 1]];
                        n1j=[n1j n1jbase+[0 1 0 1 0 1 0 0]];
                        n1k=[n1k n1kbase+[0 0 1 1 0 0 1 0]];
                        
                        n2i=[n2i n1ibase+[1 1 1 1 1 1 1 1]];
                        n2j=[n2j n1jbase+[0 1 0 1 1 1 1 0]];
                        n2k=[n2k n1kbase+[0 0 1 1 0 1 1 1]];
                        
                        n1ibase_arr=[n1ibase_arr n1ibase];
                        n1jbase_arr=[n1jbase_arr n1jbase];
                        n1kbase_arr=[n1kbase_arr n1kbase];
                        
                        delay=[delay repmat(net.r_var(s).delayi(i),1,8)];
                        
                        n1ibase=n1ibase+1;
                    end
                end
            end
        end
        net.r_var(s).n1=xyz2node_fn(n1i,n1j,n1k,Ni,Nj,Nk);
        net.r_var(s).n2=xyz2node_fn(n2i,n2j,n2k,Ni,Nj,Nk);
        net.r_var(s).n1i=n1i;
        net.r_var(s).n1j=n1j;
        net.r_var(s).n1k=n1k;
        net.r_var(s).n2i=n2i;
        net.r_var(s).n2j=n2j;
        net.r_var(s).n2k=n2k;
        net.r_var(s).delay=delay;
        net.r_var(s).n1ibase=n1ibase_arr;
        net.r_var(s).n1jbase=n1jbase_arr;
        net.r_var(s).n1kbase=n1kbase_arr;        
    end
end

% Current Sources Nodes
if net.param.i_num==2 && net.param.i_opp==0
    net.curr(1).n1i=[round((Ni+1)/2-net.param.i_elec_spacing/2) round((Ni+1)/2-net.param.i_elec_spacing/2)];
    net.curr(1).n1j=[round((Nj+1)/2-net.param.i_src_dist/2) round((Nj+1)/2+net.param.i_src_dist/2)];
    net.curr(1).n1k=1;
    net.curr(1).n2i=[round((Ni+1)/2+net.param.i_elec_spacing/2) round((Ni+1)/2+net.param.i_elec_spacing/2)];
    net.curr(1).n2j=[round((Nj+1)/2-net.param.i_src_dist/2) round((Nj+1)/2+net.param.i_src_dist/2)];
    net.curr(1).n2k=1;
elseif net.param.i_num==2 && net.param.i_opp==1
    net.curr(1).n1i=net.param.i_src_n1i;
    net.curr(1).n1j=[round((Nj+1)/2-net.param.i_src_dist/2) round((Nj+1)/2+net.param.i_src_dist/2)];
    net.curr(1).n1k=1;
    net.curr(1).n2i=net.param.i_src_n2i;
    net.curr(1).n2j=[round((Nj+1)/2-net.param.i_src_dist/2) round((Nj+1)/2+net.param.i_src_dist/2)];
    net.curr(1).n2k=1;
else
    net.curr(1).n1i_1=net.param.i_src_n1i;
    net.curr(1).n1j_1=net.param.i_src_n1j;
    net.curr(1).n1k_1=net.param.i_src_n1k;
    
    net.curr(1).n2i_1=net.param.i_src_n2i;
    net.curr(1).n2j_1=net.param.i_src_n2j;
    net.curr(1).n2k_1=net.param.i_src_n2k;
end

for j=1:length(net.curr)
    net.curr(j).n1i=[];
    net.curr(j).n1j=[];
    net.curr(j).n2i=[];
    net.curr(j).n2j=[];
    net.curr(j).n1k=[];
    net.curr(j).n2k=[];
    for i=1:net.curr(j).dia_y+1
        for m=1:net.curr(j).dia_x+1
            net.curr(j).n1i=[net.curr(j).n1i net.curr(j).n1i_1+i-1];
            net.curr(j).n1j=[net.curr(j).n1j net.curr(j).n1j_1+m-1];
            net.curr(j).n1k=[net.curr(j).n1k net.curr(j).n1k_1];
            net.curr(j).n2i=[net.curr(j).n2i net.curr(j).n2i_1+i-1];
            net.curr(j).n2j=[net.curr(j).n2j net.curr(j).n2j_1+m-1];
            net.curr(j).n2k=[net.curr(j).n2k net.curr(j).n2k_1];
        end
    end
    net.curr(j).n1=xyz2node_fn(net.curr(j).n1i,net.curr(j).n1j,net.curr(j).n1k,Ni,Nj,Nk);
    net.curr(j).n2=xyz2node_fn(net.curr(j).n2i,net.curr(j).n2j,net.curr(j).n2k,Ni,Nj,Nk);
end

if ~isempty(net.param.s_n1i)
    for j=1:length(net.param.s_n1i)
        net.s(j).n1i=[];
        net.s(j).n1j=[];
        net.s(j).n2i=[];
        net.s(j).n2j=[];
        net.s(j).n1k=[];
        net.s(j).n2k=[];
        for i=1:net.param.s_dia_y+1
            for m=1:net.param.s_dia_x+1
                net.s(j).n1i=[net.s(j).n1i net.param.s_n1i(j)+i-1];
                net.s(j).n1j=[net.s(j).n1j net.param.s_n1j(j)+m-1];
                net.s(j).n2i=[net.s(j).n2i net.param.s_n2i(j)+i-1];
                net.s(j).n2j=[net.s(j).n2j net.param.s_n2j(j)+m-1];
                net.s(j).n1k=[net.s(j).n1k net.param.s_n1k(j)];
                net.s(j).n2k=[net.s(j).n2k net.param.s_n2k(j)];
            end
        end
        net.s(j).n1=xyz2node_fn(net.s(j).n1i,net.s(j).n1j,net.s(j).n1k,Ni,Nj,Nk);
        net.s(j).n2=xyz2node_fn(net.s(j).n2i,net.s(j).n2j,net.s(j).n2k,Ni,Nj,Nk);
    end
end
%fprintf('Nodes created:');toc;

% Node Array NaN
%tic
for j=1:length(net.curr)
    for i=1:length(net.curr(j).n1i)
        node_arr(net.curr(j).n1i(i),net.curr(j).n1j(i),net.curr(j).n1k(i))=NaN;
        node_arr(net.curr(j).n2i(i),net.curr(j).n2j(i),net.curr(j).n2k(i))=NaN;
    end
end
for j=1:length(net.curr)
    n1imax=max(net.curr(j).n1i);n1jmax=max(net.curr(j).n1j);
    for i=1:length(net.curr(j).n1i)
        if net.curr(j).n1i(i)<n1imax && net.curr(j).n1j(i)<n1jmax
            img_arr(net.curr(j).n1i(i),net.curr(j).n1j(i),net.curr(j).n1k(i))=NaN;
            img_arr(net.curr(j).n2i(i),net.curr(j).n2j(i),net.curr(j).n2k(i))=NaN;
        end
    end
end
%fprintf('Nodes NaN current:');toc;
%tic
if net.param.re_var_max>0
    for j=1:length(net.r_var)
        for i=1:length(net.r_var(j).n1i)           
            node_arr(net.r_var(j).n1i(i),net.r_var(j).n1j(i),net.r_var(j).n1k(i))=NaN;
            node_arr(net.r_var(j).n2i(i),net.r_var(j).n2j(i),net.r_var(j).n2k(i))=NaN;
        end
    end
end
r_var_node_arr=[];
if net.param.re_var_max>0
%     for j=1:length(net.r_var)
%         for i=1:length(net.r_var(j).n1i)
%             img_arr(net.r_var(j).n1i(i),net.r_var(j).n1j(i),net.r_var(j).n1k(i))=-inf;
%             img_arr(net.r_var(j).n2i(i),net.r_var(j).n2j(i),net.r_var(j).n2k(i))=-inf;
%         end
%     end
    for j=1:length(net.r_var)
        for i=1:length(net.r_var(j).n1ibase)           
            %img_arr(net.r_var(j).n1ibase(i),net.r_var(j).n1jbase(i),net.r_var(j).n1kbase(i))=-inf;         
            img_arr(net.r_var(j).n1ibase(i),net.r_var(j).n1jbase(i),net.r_var(j).n1kbase(i))=abs(1./((1./net.param.re_art_dc)+1./(net.param.ri_art_dc+1./(1i*(2*pi*10e3)*net.param.cm_art_dc))));             
        end
        img_arr(net.r_var(j).n1ibase(i)+1,net.r_var(j).n1jbase(i),net.r_var(j).n1kbase(i))=abs(1./((1./net.param.re_art_dc)+1./(net.param.ri_art_dc+1./(1i*(2*pi*10e3)*net.param.cm_art_dc))));
    end
end
%fprintf('Nodes NaN Art:');toc;
%tic
if ~isempty(net.param.s_n1i)
    for j=1:length(net.s)
        for i=1:length(net.s(j).n1i)
            node_arr(net.s(j).n1i(i),net.s(j).n1j(i),net.s(j).n1k(i))=NaN;
            node_arr(net.s(j).n2i(i),net.s(j).n2j(i),net.s(j).n2k(i))=NaN;
        end
    end
end
if ~isempty(net.param.s_n1i)
    for j=1:length(net.s)
        n1imax=max(net.s(j).n1i);n1jmax=max(net.s(j).n1j);
        for i=1:length(net.s(j).n1i)
            if net.s(j).n1i(i)<n1imax && net.s(j).n1j(i)<n1jmax
            img_arr(net.s(j).n1i(i),net.s(j).n1j(i),net.s(j).n1k(i))=inf;
            img_arr(net.s(j).n2i(i),net.s(j).n2j(i),net.s(j).n2k(i))=inf;
            end
        end
    end
end
net.node_arr=node_arr;
net.img_arr=img_arr;
net.img_arr_nonan=img_arr_nonan;


fprintf('Nodes Created:');toc;
tic
fileID=[];
if net.param.debug==0
    [fileID errmsg] = fopen(file_name,'w');
    if ~isempty(errmsg)
        fprintf(errmsg);
    end
    while isempty(fileID)||(fileID==-1)
        [fileID errmsg] = fopen(file_name,'w');
        if ~isempty(errmsg)
            fprintf(errmsg);
        end
        pause(0.1);
    end
    %disp(fileID);
    fprintf(fileID,'bioz netlist\n');
    re_art_all_arr=[];
    ri_art_all_arr=[];
    cm_art_all_arr=[];
    % Horizontal Resistors (x)
    nodes_horz=[];
    for i=1:Ni+1
        for j=1:Nj
            for k=1:Nk+1
                n1=xyz2node_fn(i,j,k,Ni,Nj,Nk);
                n2=n1+1;
                nodes_horz=[nodes_horz;[n1, n2]];
                found=0;
                if net.param.re_var_max(1)>0
                    for s=1:length(net.r_var)
                        l_all=find(n1==net.r_var(s).n1&n2==net.r_var(s).n2);
                        if length(l_all)>1
                            l_new=l_all(1);
                        else
                            l_new=l_all;
                        end
                        l=l_new;
                        if ~isempty(l)
                            found=1;
                            if net.param.re_var_max>0
                                if net.param.z_en==0
                                    if net.param.mode==0
                                        fprintf(fileID,'r_%g_%g %g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,net.param.re_art_dc,net.param.re_var_max(s),freq(s),net.r_var(s).delay(l));
                                    else
                                        re_art_all=net.param.re_art_dc-net.param.re_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var(s).delay(l))-pi/2);
                                        re_art_all_arr=[re_art_all_arr; re_art_all];
                                        fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,re_art_all);
                                    end
                                else
                                    if net.param.mode==0
                                        fprintf(fileID,'r_%g_%g_e %g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,net.param.re_art_dc,net.param.re_var_max(s),freq(s),net.r_var(s).delay(s));    %delta Re
                                        fprintf(fileID,'r_%g_%g_i %g_%g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,n2,net.param.ri_art_dc,net.param.ri_var_max(s),freq(s),net.r_var(s).delay(l)); %delta Ri
                                        fprintf(fileID,'c_%g_%g_i %g %g_%g Q=(%g+%g*sin(2*pi*%g*(time+%g)-pi/2))*x\n',n1,n2,n1,n1,n2,net.param.cm_art_dc,net.param.cm_var_max(s),freq(s),net.r_var(s).delay(l)); %delta Ci
                                    else
                                        re_art_all=net.param.re_art_dc-net.param.re_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var(s).delay(l))-pi/2);
                                        ri_art_all=net.param.ri_art_dc-net.param.ri_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var(s).delay(l))-pi/2);
                                        cm_art_all=net.param.cm_art_dc+net.param.cm_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var(s).delay(l))-pi/2);
                                        re_art_all_arr=[re_art_all_arr; re_art_all];
                                        ri_art_all_arr=[ri_art_all_arr; ri_art_all];
                                        cm_art_all_arr=[cm_art_all_arr; cm_art_all];
                                        fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,re_art_all);    %delta Re
                                        fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,ri_art_all); %delta Ri
                                        fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,cm_art_all); %delta Ci
                                    end
                                end
                            else
                                if net.param.z_en==0
                                    fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,net.re_art_dc);
                                else
                                    fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,net.re_art_dc);  % Re
                                    fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,net.ri_art_dc); % Ri
                                    fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,net.cm_art_dc); % Ci
                                end
                            end
                        end
                    end
                end
                if found==0
                    if net.param.z_en==0
                        fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,net.grid.x.re_dc_arr(i,j,k));
                    else
                        fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,net.grid.x.re_dc_arr(i,j,k));  % Re
                        fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,net.grid.x.ri_dc_arr(i,j,k)); % Ri
                        fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,net.grid.x.cm_dc_arr(i,j,k)); % Ci
                    end
                end
            end
        end
    end
    fprintf('Horizontal Resistors (x):');toc
    tic
    % Vertical Resistors (y)
    nodes_vert=[];
    for i=1:Ni
        for j=1:Nj+1
            for k=1:Nk+1
                n1=xyz2node_fn(i,j,k,Ni,Nj,Nk);
                n2=n1+Nj+1;
                nodes_vert=[nodes_vert;[n1, n2]];
                found=0;
                if net.param.re_var_max(1)>0
                    for s=1:length(net.r_var)
                        l_all=find(n1==net.r_var(s).n1&n2==net.r_var(s).n2);
                        if length(l_all)>1
                            l_new=l_all(1);
                        else
                            l_new=l_all;
                        end
                        l=l_new;
                        if ~isempty(l)
                            found=1;
                            if net.param.re_var_max>0
                                if net.param.z_en==0
                                    if net.param.mode==0
                                        fprintf(fileID,'r_%g_%g %g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,net.param.re_art_dc,net.param.re_var_max(s),freq(s),net.r_var(s).delay(l));
                                    else
                                        re_art_all=net.param.re_art_dc-net.param.re_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var(s).delay(l))-pi/2);
                                        re_art_all_arr=[re_art_all_arr; re_art_all];
                                        fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,re_art_all);
                                    end
                                else
                                    if net.param.mode==0
                                        fprintf(fileID,'r_%g_%g_e %g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,net.param.re_art_dc,net.param.re_var_max(s),freq(s),net.r_var(s).delay(l));    %delta Re
                                        fprintf(fileID,'r_%g_%g_i %g_%g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,n2,net.param.ri_art_dc,net.param.ri_var_max(s),freq(s),net.r_var(s).delay(l)); %delta Ri
                                        fprintf(fileID,'c_%g_%g_i %g %g_%g Q=(%g+%g*sin(2*pi*%g*(time+%g)-pi/2))*x\n',n1,n2,n1,n1,n2,net.param.cm_art_dc,net.param.cm_var_max(s),freq(s),net.r_var(s).delay(l)); %delta Ci
                                    else
                                        re_art_all=net.param.re_art_dc-net.param.re_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var(s).delay(l))-pi/2);
                                        ri_art_all=net.param.ri_art_dc-net.param.ri_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var(s).delay(l))-pi/2);
                                        cm_art_all=net.param.cm_art_dc+net.param.cm_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var(s).delay(l))-pi/2);
                                        re_art_all_arr=[re_art_all_arr; re_art_all];
                                        ri_art_all_arr=[ri_art_all_arr; ri_art_all];
                                        cm_art_all_arr=[cm_art_all_arr; cm_art_all];
                                        fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,re_art_all);    %delta Re
                                        fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,ri_art_all); %delta Ri
                                        fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,cm_art_all); %delta Ci
                                    end
                                end
                            else
                                if net.param.z_en==0
                                    fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,net.re_art_dc);
                                else
                                    fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,net.re_art_dc);  % Re
                                    fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,net.ri_art_dc); % Ri
                                    fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,net.cm_art_dc); % Ci
                                end
                            end
                        end
                    end
                end
                if found==0
                    if net.param.z_en==0
                        fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,net.grid.y.re_dc_arr(i,j,k));
                    else
                        fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,net.grid.y.re_dc_arr(i,j,k));  % Re
                        fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,net.grid.y.ri_dc_arr(i,j,k)); % Ri
                        fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,net.grid.y.cm_dc_arr(i,j,k)); % Ci
                    end
                end
            end
        end
    end
    fprintf('Vertical Resistors (y):');toc
    tic
    % Depth Resistors (z)
    nodes_depth=[];
    for i=1:Ni+1
        for j=1:Nj+1
            for k=1:Nk
                n1=xyz2node_fn(i,j,k,Ni,Nj,Nk);
                n2=n1+(Nj+1)*(Ni+1);
                nodes_depth=[nodes_depth;[n1, n2]];
                found=0;
                if net.param.re_var_max(1)>0
                    for s=1:length(net.r_var)
                        l_all=find(n1==net.r_var(s).n1&n2==net.r_var(s).n2);
                        if length(l_all)>1
                            l_new=l_all(1);
                        else
                            l_new=l_all;
                        end
                        l=l_new;
                        if ~isempty(l)
                            found=1;
                            if net.param.re_var_max>0
                                if net.param.z_en==0
                                    if net.param.mode==0
                                        fprintf(fileID,'r_%g_%g %g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,net.param.re_art_dc,net.param.re_var_max(s),freq(s),net.r_var(s).delay(l));
                                    else
                                        re_art_all=net.param.re_art_dc-net.param.re_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var(s).delay(l))-pi/2);
                                        re_art_all_arr=[re_art_all_arr; re_art_all];
                                        fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,re_art_all);
                                    end
                                else
                                    if net.param.mode==0
                                        fprintf(fileID,'r_%g_%g_e %g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,net.param.re_art_dc,net.param.re_var_max(s),freq(s),net.r_var(s).delay(l));    %delta Re
                                        fprintf(fileID,'r_%g_%g_i %g_%g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,n2,net.param.ri_art_dc,net.param.ri_var_max(s),freq(s),net.r_var(s).delay(l)); %delta Ri
                                        fprintf(fileID,'c_%g_%g_i %g %g_%g Q=(%g+%g*sin(2*pi*%g*(time+%g)-pi/2))*x\n',n1,n2,n1,n1,n2,net.param.cm_art_dc,net.param.cm_var_max(s),freq(s),net.r_var(s).delay(l)); %delta Ci
                                    else
                                        re_art_all=net.param.re_art_dc-net.param.re_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var(s).delay(l))-pi/2);
                                        ri_art_all=net.param.ri_art_dc-net.param.ri_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var(s).delay(l))-pi/2);
                                        cm_art_all=net.param.cm_art_dc+net.param.cm_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var(s).delay(l))-pi/2);
                                        re_art_all_arr=[re_art_all_arr; re_art_all];
                                        ri_art_all_arr=[ri_art_all_arr; ri_art_all];
                                        cm_art_all_arr=[cm_art_all_arr; cm_art_all];
                                        fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,re_art_all);    %delta Re
                                        fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,ri_art_all); %delta Ri
                                        fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,cm_art_all); %delta Ci
                                    end
                                end
                            else
                                if net.param.z_en==0
                                    fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,net.re_art_dc);
                                else
                                    fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,net.re_art_dc);  % Re
                                    fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,net.ri_art_dc); % Ri
                                    fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,net.cm_art_dc); % Ci
                                end
                            end
                        end
                    end
                end
                if found==0
                    if net.param.z_en==0
                        fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,net.grid.z.re_dc_arr(i,j,k));
                    else
                        fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,net.grid.z.re_dc_arr(i,j,k));  % Re
                        fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,net.grid.z.ri_dc_arr(i,j,k)); % Ri
                        fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,net.grid.z.cm_dc_arr(i,j,k)); % Ci
                    end
                end
            end
        end
    end
    if exist('net.r_var')
        net.r_var(s).re_art_all_arr=re_art_all_arr;
        net.r_var(s).ri_art_all_arr=ri_art_all_arr;
        net.r_var(s).cm_art_all_arr=cm_art_all_arr;
    end
    fprintf('Depth Resistors (z):');toc
    tic
    % voltage sensing electrodes
    if net.param.se_en==0
        if ~isempty(net.param.s_n1i)
            for j=1:length(net.param.s_n1i)
                for m=1:size(net.s(j).n2,2)-1
                    fprintf(fileID,'v_s%g_n1_%g %g %g dc 0 ac 0 0\n',j,m,net.s(j).n1(1),net.s(j).n1(m+1));
                    fprintf(fileID,'v_s%g_n2_%g %g %g dc 0 ac 0 0\n',j,m,net.s(j).n2(1),net.s(j).n2(m+1));
                end
            end
        end
    else
        if ~isempty(net.param.s_n1i)
            for j=1:length(net.param.s_n1i)
                % Create the skin-electrode imepdance
                if  net.param.z_en==0
                    for m=1:size(net.s(j).n2,2)
                        % s_n1
                        fprintf(fileID,'r1_s%g_n1_%g_se %g %g_se %g\n',j,net.s(j).n1(m),net.s(j).n1(m),net.s(j).n1(m),net.param.r1_se); % % Skin-Electrode R1
                        % s_n2
                        fprintf(fileID,'r1_s%g_n2_%g_se %g %g_se %g\n',j,net.s(j).n2(m),net.s(j).n2(m),net.s(j).n2(m),net.param.r1_se); % % Skin-Electrode R1
                        
                    end
                else
                    for m=1:size(net.s(j).n2,2)
                        % s_n1
                        fprintf(fileID,'r1_s%g_n1_%g_se %g %g_se_i %g\n',j,net.s(j).n1(m),net.s(j).n1(m),net.s(j).n1(m),net.param.r1_se); % % Skin-Electrode R1
                        fprintf(fileID,'c1_s%g_n1_%g_se %g %g_se_i %g\n',j,net.s(j).n1(m),net.s(j).n1(m),net.s(j).n1(m),net.param.c1_se); % % Skin-Electrode C1
                        fprintf(fileID,'r2_s%g_n1_%g_se %g_se %g_se_i %g\n',j,net.s(j).n1(m),net.s(j).n1(m),net.s(j).n1(m),net.param.r2_se); % % Skin-Electrode R2
                        fprintf(fileID,'c2_s%g_n1_%g_se %g_se %g_se_i %g\n',j,net.s(j).n1(m),net.s(j).n1(m),net.s(j).n1(m),net.param.c2_se); % % Skin-Electrode C2
                        % s_n2
                        fprintf(fileID,'r1_s%g_n2_%g_se %g %g_se_i %g\n',j,net.s(j).n2(m),net.s(j).n2(m),net.s(j).n2(m),net.param.r1_se); % % Skin-Electrode R1
                        fprintf(fileID,'c1_s%g_n2_%g_se %g %g_se_i %g\n',j,net.s(j).n2(m),net.s(j).n2(m),net.s(j).n2(m),net.param.c1_se); % % Skin-Electrode C1
                        fprintf(fileID,'r2_s%g_n2_%g_se %g_se %g_se_i %g\n',j,net.s(j).n2(m),net.s(j).n2(m),net.s(j).n2(m),net.param.r2_se); % % Skin-Electrode R2
                        fprintf(fileID,'c2_s%g_n2_%g_se %g_se %g_se_i %g\n',j,net.s(j).n2(m),net.s(j).n2(m),net.s(j).n2(m),net.param.c2_se); % % Skin-Electrode C2
                    end
                end
                % Short circuit voltage sensing electrodes
                for m=1:size(net.s(j).n2,2)-1
                    fprintf(fileID,'v_s%g_n1_%g %g_se %g_se dc 0 ac 0 0\n',j,m,net.s(j).n1(1),net.s(j).n1(m+1));
                    fprintf(fileID,'v_s%g_n2_%g %g_se %g_se dc 0 ac 0 0\n',j,m,net.s(j).n2(1),net.s(j).n2(m+1));
                end
            end
        end
    end
    
    % Current source
    for i=1:length(net.curr)
        %if net.curr(i).iac>0
            curr_n1=net.curr(i).n1;
            curr_n2=net.curr(i).n2;
%         else
%             curr_n1=net.curr(i).n2;
%             curr_n2=net.curr(i).n1;
%         end
        if net.param.se_en==0
            % Short Ciruit current electrodes
            for j=1:size(net.curr(i).n2,2)-1
                fprintf(fileID,'v_i%g_n1_%g %g %g dc 0 ac 0 0\n',i,j,curr_n1(1),curr_n1(j+1));
                fprintf(fileID,'v_i%g_n2_%g %g %g dc 0 ac 0 0\n',i,j,curr_n2(1),curr_n2(j+1));
            end
            % Current source
            fprintf(fileID,'i%g %g %g dc %g ac %g 0\n',i,curr_n2(1),curr_n1(1),(net.curr(i).idc),(net.curr(i).iac));
            
        else % Skin-Electrode
            % Create the skin-electrode imepdance
            if  net.param.z_en==0
                for j=1:size(net.curr(i).n1,2)
                    % curr_n1
                    fprintf(fileID,'r1_i%g_n1_%g_se %g %g_se %g\n',i,curr_n1(j),curr_n1(j),curr_n1(j),net.param.r1_se); % % Skin-Electrode R1
                    % curr_n2
                    fprintf(fileID,'r1_i%g_n2_%g_se %g %g_se %g\n',i,curr_n2(j),curr_n2(j),curr_n2(j),net.param.r1_se); % % Skin-Electrode R1                    
                end
            else
                for j=1:size(net.curr(i).n1,2)
                    % curr_n1
                    fprintf(fileID,'r1_i%g_n1_%g_se %g %g_se_i %g\n',i,curr_n1(j),curr_n1(j),curr_n1(j),net.param.r1_se); % % Skin-Electrode R1
                    fprintf(fileID,'c1_i%g_n1_%g_se %g %g_se_i %g\n',i,curr_n1(j),curr_n1(j),curr_n1(j),net.param.c1_se); % % Skin-Electrode C1
                    fprintf(fileID,'r2_i%g_n1_%g_se %g_se %g_se_i %g\n',i,curr_n1(j),curr_n1(j),curr_n1(j),net.param.r2_se); % % Skin-Electrode R2
                    fprintf(fileID,'c2_i%g_n1_%g_se %g_se %g_se_i %g\n',i,curr_n1(j),curr_n1(j),curr_n1(j),net.param.c2_se); % % Skin-Electrode C2
                    % curr_n2
                    fprintf(fileID,'r1_i%g_n2_%g_se %g %g_se_i %g\n',i,curr_n2(j),curr_n2(j),curr_n2(j),net.param.r1_se); % % Skin-Electrode R1
                    fprintf(fileID,'c1_i%g_n2_%g_se %g %g_se_i %g\n',i,curr_n2(j),curr_n2(j),curr_n2(j),net.param.c1_se); % % Skin-Electrode C1
                    fprintf(fileID,'r2_i%g_n2_%g_se %g_se %g_se_i %g\n',i,curr_n2(j),curr_n2(j),curr_n2(j),net.param.r2_se); % % Skin-Electrode R2
                    fprintf(fileID,'c2_i%g_n2_%g_se %g_se %g_se_i %g\n',i,curr_n2(j),curr_n2(j),curr_n2(j),net.param.c2_se); % % Skin-Electrode C2
                end
            end
            
            % Short Ciruit current electrodes
            l=0;
            for j=1:size(net.curr(i).n2,2)-1
                l=l+1;
                fprintf(fileID,'v_i%g_n1_%g_se %g_se %g_se dc 0 ac 0 0\n',i,l,curr_n1(1),curr_n1(j+1));
                l=l+1;
                fprintf(fileID,'v_i%g_n2_%g_se %g_se %g_se dc 0 ac 0 0\n',i,l,curr_n2(1),curr_n2(j+1));
            end
            % Current source
            fprintf(fileID,'i%g %g_se %g_se dc %g ac %g 0\n',i,curr_n1(1),curr_n2(1),(net.curr(i).idc),(net.curr(i).iac));
        end
    end
    
    % GND short circuit
    if net.param.se_en==0
        fprintf(fileID,'vgnd 0 %g dc 0 ac 0 0\n',curr_n2(1));
    else
        fprintf(fileID,'vgnd 0 %g_se dc 0 ac 0 0\n',curr_n2(1));
    end
    if net.param.mode==0
        % Step Analysis
        if net.param.re_var_max>0
            fprintf(fileID,'.tran %g %g\n',tstep,tend);
        else
            fprintf(fileID,'.dc i1 %g %g 1\n',(net.curr(1).idc),(net.curr(1).idc));
        end
    else
        ac_freq_type_str={'lin', 'oct', 'dec'};
        fprintf(fileID,'.ac %s %g %g %g\n',ac_freq_type_str{net.param.ac_freq(4)+1},net.param.ac_freq(1),net.param.ac_freq(2),net.param.ac_freq(3));  %n_pts, fstart,fend
    end
    fprintf(fileID,'.options reltol=%g abstol=%g Noopiter GminSteps=0\n',net.param.reltol,net.param.abstol);
    fprintf(fileID,'.backanno\n');
    fprintf(fileID,'.end\n');
    fclose(fileID);
    fprintf('... Netlist Generated')
else
    %     if net.param.plot_eit==2
    %         figure;
    %         sgtitle(strcat('Figure '))
    %         subplot(1,2,1)
    %         %heatmap([0:net.Nj],[0:net.Ni],squeeze(net.node_arr(:,:,1)),'Colormap',white,'MissingDataLabel','Elec.');
    %         mesh_plot_fn(net.img_arr(:,:,1),1)
    %         Vs_arrow_plot_fn(net)
    %         xlabel('X');ylabel('Y')
    %         title('Impedance Image Top View, I(Blue), V(Green), Z=0')
    %         if net.Nk>0
    %             subplot(1,2,2)
    %             %heatmap([0:net.Nj],[0:net.Nk],squeeze(net.node_arr(1,:,:))','Colormap',white,'MissingDataLabel','Artery');
    %             mesh_plot_fn(net.img_arr(1,:,:))
    %             xlabel('X');ylabel('Z')
    %             title('Impedance Image Cross Section, I(Blue), V(Green), Art.(red)')
    %         end
    %         set(gcf,'position',[0.0010    0.0410    1.2800    0.6073]*1e3);
    %         %set(gcf,'position',[0.1026    0.3058    1.3656    0.4624]*1e3)
    %     end
end
