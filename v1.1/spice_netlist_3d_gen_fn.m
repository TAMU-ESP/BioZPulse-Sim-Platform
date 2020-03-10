function net=spice_netlist_3d_gen_fn(net)

%   Bassem Ibrahim
%   2019.08.22
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%   
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation 
%   Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling, 
%   IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan.
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
% 8     net.param.r_art_dc		Artery's extra-cellular resistance (RE)
% 9     net.param.ri_art_dc		Artery's intra-cellular resistance (RI)
% 10	net.param.ci_art_dc		Artery's cell membrane capacitance (Cm)
% 11	net.param.r_var_max		Artery's variable extra-cellular resistance (?RE)
% 12	net.param.ri_var_max	Artery's variable intra-cellular resistance (?RI)
% 13	net.param.ci_var_max	Artery's variable cell membrane capacitance (?Cm)
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
% 39	net.curr.iac			AC current amplitude
% 40	net.curr.freq			AC current frequency
% 30	net.curr.dia_x			Size of current electrode in X-axis
% 31	net.curr.dia_y			Size of current electrode in Y-axis
% 41	net.param.r_var_time	Simulation time
% 42	net.param.tstep			Simulation time step




file_name=net.param.circuit_file;
if exist(file_name,'file')
    delete(file_name);
end

r_dc=net.param.r_dc;
if length(r_dc)==3
    r_dc_x=r_dc(1,1);
    r_dc_y=r_dc(1,2);
    r_dc_z=r_dc(1,3);
elseif length(r_dc)==1
    r_dc_x=r_dc;
    r_dc_y=r_dc;
    r_dc_z=r_dc;
    % Edge resistors
    r_dc_x_e=r_dc;
    r_dc_y_e=r_dc;
    r_dc_z_e=r_dc;
    if net.param.small_edges==1
        r_dc_x_e=r_dc/2;
        r_dc_y_e=r_dc/2;
        r_dc_z_e=r_dc/2;
    end
end

Ni=net.Ni;
Nj=net.Nj;
Nk=net.Nk;

freq=net.param.r_var_freq;
tend=net.param.r_var_time;
td=net.param.r_var_delay;
tstep=net.param.tstep;

% Create Node Array
for i=1:Ni+1
    for j=1:Nj+1
        for k=1:Nk+1
            node_arr(i,j,k)=xyz2node_fn(i,j,k,Ni,Nj,Nk);
        end
    end
end
% Variable Resistor nodes
net.r_var.delayi=[];
if net.param.r_var_max(1)>0
    for s=1:length(net.param.r_var_n1j)
        delay_line=linspace(0,td,net.param.r_var_seg(s));
        S=floor((net.param.r_var_num(s))/net.param.r_var_seg(s));
        net.r_var.delayi(s,:)=reshape(repmat(delay_line,S,1),1,[]);
        net.r_var.delayi(s,:)=[net.r_var.delayi(s,:) repmat(net.r_var.delayi(s,end),1,net.Ni-S*net.param.r_var_seg(s))];
        if net.param.r_var_dim==1
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
        elseif net.param.r_var_dim==2
            %             n1i=net.param.r_var_n1i+1+[0 0 1 0];
            %             n1j=net.param.r_var_n1j+1+[0 1 0 0];
            %             n2i=net.param.r_var_n1i+1+[0 1 1 1];
            %             n2j=net.param.r_var_n1j+1+[1 1 1 0];
            %             n1k=net.param.r_var_n1k+1+[0 0 0 0];
            %             n2k=net.param.r_var_n1k+1+[0 0 0 0];
            %             delay=[0,0,0,0];
            n1i=[];n1j=[];n1k=[];n2i=[];n2j=[];n2k=[];delay=[];
            for j=1:net.param.r_var_dia
                n1ibase=net.param.r_var_n1i(s)+1;
                n1jbase=net.param.r_var_n1j(s)+1+j-1;
                n1kbase=net.param.r_var_n1k(s)+1;

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
        elseif net.param.r_var_dim==3
            n1i=[];n1j=[];n1k=[];n2i=[];n2j=[];n2k=[];delay=[];
            for j=1:net.param.r_var_dia
                for k=1:net.param.r_var_dia
                    n1ibase=net.param.r_var_n1i(s)+1;
                    n1jbase=net.param.r_var_n1j(s)+1+j-1;
                    n1kbase=net.param.r_var_n1k(s)+1+k-1;
                    
                    n1i=[n1i n1ibase+[0 0 0 0]];
                    n1j=[n1j n1jbase+[0 1 0 0]];
                    n1k=[n1k n1kbase+[0 0 1 0]];
                    
                    n2i=[n2i n1ibase+[0 0 0 0]];
                    n2j=[n2j n1jbase+[1 1 1 0]];
                    n2k=[n2k n1kbase+[0 1 1 1]];
                    delay=[delay repmat(net.r_var.delayi(s,1),1,4)];
                    for i=1:net.param.r_var_num(s)
                        n1i=[n1i n1ibase+[0 0 0 0 1 1 1 1]];
                        n1j=[n1j n1jbase+[0 1 0 1 0 1 0 0]];
                        n1k=[n1k n1kbase+[0 0 1 1 0 0 1 0]];
                        
                        n2i=[n2i n1ibase+[1 1 1 1 1 1 1 1]];
                        n2j=[n2j n1jbase+[0 1 0 1 1 1 1 0]];
                        n2k=[n2k n1kbase+[0 0 1 1 0 1 1 1]];
                        
                        delay=[delay repmat(net.r_var.delayi(s,i),1,8)];
                        
                        n1ibase=n1ibase+1;
                    end
                end
            end
        end
        net.r_var.n1(s,:)=xyz2node_fn(n1i,n1j,n1k,Ni,Nj,Nk);
        net.r_var.n2(s,:)=xyz2node_fn(n2i,n2j,n2k,Ni,Nj,Nk);
        net.r_var.n1i(s,:)=n1i;
        net.r_var.n1j(s,:)=n1j;
        net.r_var.n1k(s,:)=n1k;
        net.r_var.n2i(s,:)=n2i;
        net.r_var.n2j(s,:)=n2j;
        net.r_var.n2k(s,:)=n2k;
        net.r_var.delay(s,:)=delay;
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
    net.curr(1).n1i_1=net.param.i_src_n1i+1;
    net.curr(1).n1j_1=net.param.i_src_n1j+1;
    net.curr(1).n1k=0+1;
    
    net.curr(1).n2i_1=net.param.i_src_n2i+1;
    net.curr(1).n2j_1=net.param.i_src_n2j+1;
    net.curr(1).n2k=0+1;
end


for j=1:length(net.curr)
    net.curr(j).n1i=[];
    net.curr(j).n1j=[];
    net.curr(j).n2i=[];
    net.curr(j).n2j=[];
    for i=1:net.curr(j).dia_y+1
        for m=1:net.curr(j).dia_x+1
            net.curr(j).n1i=[net.curr(j).n1i net.curr(j).n1i_1+i-1];
            net.curr(j).n1j=[net.curr(j).n1j net.curr(j).n1j_1+m-1];
            net.curr(j).n2i=[net.curr(j).n2i net.curr(j).n2i_1+i-1];
            net.curr(j).n2j=[net.curr(j).n2j net.curr(j).n2j_1+m-1];
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
                net.s(j).n1i=[net.s(j).n1i net.param.s_n1i(j)+1+i-1];
                net.s(j).n1j=[net.s(j).n1j net.param.s_n1j(j)+1+m-1];
                net.s(j).n2i=[net.s(j).n2i net.param.s_n2i(j)+1+i-1];
                net.s(j).n2j=[net.s(j).n2j net.param.s_n2j(j)+1+m-1];
                net.s(j).n1k=[net.s(j).n1k 0+1];
                net.s(j).n2k=[net.s(j).n2k 0+1];
            end
        end
        net.s(j).n1=xyz2node_fn(net.s(j).n1i,net.s(j).n1j,net.s(j).n1k,Ni,Nj,Nk);
        net.s(j).n2=xyz2node_fn(net.s(j).n2i,net.s(j).n2j,net.s(j).n2k,Ni,Nj,Nk);
    end
end

% Node Array NaN
for i=1:length(net.curr)
    node_arr(net.curr(i).n1i,net.curr(i).n1j,net.curr(i).n1k)=NaN;
    node_arr(net.curr(i).n2i,net.curr(i).n2j,net.curr(i).n2k)=NaN;
end
if net.param.r_var_max>0
    node_arr(net.r_var.n1i,net.r_var.n1j,net.r_var.n1k)=NaN;
    node_arr(n2i,n2j,n2k)=NaN;
end
if ~isempty(net.param.s_n1i)
    for j=1:length(net.s)
        node_arr(net.s(j).n1i,net.s(j).n1j,net.s(j).n1k)=NaN;
        node_arr(net.s(j).n2i,net.s(j).n2j,net.s(j).n2k)=NaN;
    end
end
net.node_arr=node_arr;
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
    r_art_all_arr=[];
    ri_art_all_arr=[];
    ci_art_all_arr=[];
    % Horizontal Resistors (x)
    nodes_horz=[];
    for i=1:Ni+1
        for j=1:Nj
            for k=1:Nk+1
                n1=xyz2node_fn(i,j,k,Ni,Nj,Nk);
                n2=n1+1;
                nodes_horz=[nodes_horz;[n1, n2]];
                found=0;
                for s=1:size(net.r_var.n1,1)
                    for l=1:size(net.r_var.n1,2)
                        if n1==net.r_var.n1(s,l) && n2==net.r_var.n2(s,l) && found==0
                            found=1;
                            if net.param.r_var_max>0
                                if net.param.z_en==0
                                    if net.param.mode==0
                                        fprintf(fileID,'r_%g_%g %g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,net.param.r_art_dc,net.param.r_var_max(s),freq(s),net.r_var.delay(s,l));
                                    else
                                        r_art_all=net.param.r_art_dc-net.param.r_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var.delay(s,l))-pi/2);
                                        r_art_all_arr=[r_art_all_arr; r_art_all];
                                        fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,r_art_all);
                                    end
                                else
                                    if net.param.mode==0
                                        fprintf(fileID,'r_%g_%g_e %g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,net.param.r_art_dc,net.param.r_var_max(s),freq(s),net.r_var.delay(s,l));    %delta Re
                                        fprintf(fileID,'r_%g_%g_i %g_%g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,n2,net.param.ri_art_dc,net.param.ri_var_max(s),freq(s),net.r_var.delay(s,l)); %delta Ri
                                        fprintf(fileID,'c_%g_%g_i %g %g_%g Q=(%g-%g*sin(2*pi*%g*(time+%g)-pi/2))*x\n',n1,n2,n1,n1,n2,net.param.ci_art_dc,net.param.ci_var_max(s),freq(s),net.r_var.delay(s,l)); %delta Ci
                                    else
                                        r_art_all=net.param.r_art_dc-net.param.r_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var.delay(s,l))-pi/2);
                                        ri_art_all=net.param.ri_art_dc-net.param.ri_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var.delay(s,l))-pi/2);
                                        ci_art_all=net.param.ci_art_dc-net.param.ci_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var.delay(s,l))-pi/2);
                                        r_art_all_arr=[r_art_all_arr; r_art_all];
                                        ri_art_all_arr=[ri_art_all_arr; ri_art_all];
                                        ci_art_all_arr=[ci_art_all_arr; ci_art_all];
                                        fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,r_art_all);    %delta Re
                                        fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,ri_art_all); %delta Ri
                                        fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,ci_art_all); %delta Ci
                                    end
                                end
                            else
                                if net.param.z_en==0
                                    fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,net.r_art_dc);
                                else
                                    fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,net.r_art_dc);  % Re
                                    fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,net.ri_art_dc); % Ri
                                    fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,net.ci_art_dc); % Ci
                                end
                            end
                        end
                    end
                end
                if found==0
                    if net.param.z_en==0
                        fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,net.r_dc_arr(i,j,k).x);
                    else
                        fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,net.r_dc_arr(i,j,k).x);  % Re
                        fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,net.ri_dc_arr(i,j,k).x); % Ri
                        fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,net.ci_dc_arr(i,j,k).x); % Ci
                    end
                end
            end
        end
    end
    % Vertical Resistors (y)
    nodes_vert=[];
    for i=1:Ni
        for j=1:Nj+1
            for k=1:Nk+1
                n1=xyz2node_fn(i,j,k,Ni,Nj,Nk);
                n2=n1+Nj+1;
                nodes_vert=[nodes_vert;[n1, n2]];
                found=0;
                for s=1:size(net.r_var.n1,1)
                    for l=1:size(net.r_var.n1,2)
                        if n1==net.r_var.n1(s,l) && n2==net.r_var.n2(s,l) && found==0
                            found=1;
                            if net.param.r_var_max>0
                                if net.param.z_en==0
                                    if net.param.mode==0
                                        fprintf(fileID,'r_%g_%g %g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,net.param.r_art_dc,net.param.r_var_max(s),freq(s),net.r_var.delay(s,l));
                                    else
                                        r_art_all=net.param.r_art_dc-net.param.r_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var.delay(s,l))-pi/2);
                                        r_art_all_arr=[r_art_all_arr; r_art_all];
                                        fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,r_art_all);
                                    end
                                else
                                    if net.param.mode==0
                                        fprintf(fileID,'r_%g_%g_e %g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,net.param.r_art_dc,net.param.r_var_max(s),freq(s),net.r_var.delay(s,l));    %delta Re
                                        fprintf(fileID,'r_%g_%g_i %g_%g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,n2,net.param.ri_art_dc,net.param.ri_var_max(s),freq(s),net.r_var.delay(s,l)); %delta Ri
                                        fprintf(fileID,'c_%g_%g_i %g %g_%g Q=(%g-%g*sin(2*pi*%g*(time+%g)-pi/2))*x\n',n1,n2,n1,n1,n2,net.param.ci_art_dc,net.param.ci_var_max(s),freq(s),net.r_var.delay(s,l)); %delta Ci
                                    else
                                        r_art_all=net.param.r_art_dc-net.param.r_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var.delay(s,l))-pi/2);
                                        ri_art_all=net.param.ri_art_dc-net.param.ri_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var.delay(s,l))-pi/2);
                                        ci_art_all=net.param.ci_art_dc-net.param.ci_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var.delay(s,l))-pi/2);
                                        r_art_all_arr=[r_art_all_arr; r_art_all];
                                        ri_art_all_arr=[ri_art_all_arr; ri_art_all];
                                        ci_art_all_arr=[ci_art_all_arr; ci_art_all];
                                        fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,r_art_all);    %delta Re
                                        fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,ri_art_all); %delta Ri
                                        fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,ci_art_all); %delta Ci
                                    end
                                end
                            else
                                if net.param.z_en==0
                                    fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,net.r_art_dc);
                                else
                                    fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,net.r_art_dc);  % Re
                                    fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,net.ri_art_dc); % Ri
                                    fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,net.ci_art_dc); % Ci
                                end
                            end
                        end
                    end
                end
                if found==0
                    if net.param.z_en==0
                        fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,net.r_dc_arr(i,j,k).y);
                    else
                        fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,net.r_dc_arr(i,j,k).y);  % Re
                        fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,net.ri_dc_arr(i,j,k).y); % Ri
                        fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,net.ci_dc_arr(i,j,k).y); % Ci
                    end
                end
            end
        end
    end
    % Depth Resistors (z)
    nodes_depth=[];
    for i=1:Ni+1
        for j=1:Nj+1
            for k=1:Nk
                n1=xyz2node_fn(i,j,k,Ni,Nj,Nk);
                n2=n1+(Nj+1)*(Ni+1);
                nodes_depth=[nodes_depth;[n1, n2]];
                found=0;
                for s=1:size(net.r_var.n1,1)
                    for l=1:size(net.r_var.n1,2)
                        if n1==net.r_var.n1(s,l) && n2==net.r_var.n2(s,l) && found==0
                            found=1;
                            if net.param.r_var_max>0
                                if net.param.z_en==0
                                    if net.param.mode==0
                                        fprintf(fileID,'r_%g_%g %g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,net.param.r_art_dc,net.param.r_var_max(s),freq(s),net.r_var.delay(s,l));
                                    else
                                        r_art_all=net.param.r_art_dc-net.param.r_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var.delay(s,l))-pi/2);
                                        r_art_all_arr=[r_art_all_arr; r_art_all];
                                        fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,r_art_all);
                                    end
                                else
                                    if net.param.mode==0
                                        fprintf(fileID,'r_%g_%g_e %g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,net.param.r_art_dc,net.param.r_var_max(s),freq(s),net.r_var.delay(s,l));    %delta Re
                                        fprintf(fileID,'r_%g_%g_i %g_%g %g R=%g-%g*sin(2*pi*%g*(time+%g)-pi/2)\n',n1,n2,n1,n2,n2,net.param.ri_art_dc,net.param.ri_var_max(s),freq(s),net.r_var.delay(s,l)); %delta Ri
                                        fprintf(fileID,'c_%g_%g_i %g %g_%g Q=(%g-%g*sin(2*pi*%g*(time+%g)-pi/2))*x\n',n1,n2,n1,n1,n2,net.param.ci_art_dc,net.param.ci_var_max(s),freq(s),net.r_var.delay(s,l)); %delta Ci
                                    else
                                        r_art_all=net.param.r_art_dc-net.param.r_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var.delay(s,l))-pi/2);
                                        ri_art_all=net.param.ri_art_dc-net.param.ri_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var.delay(s,l))-pi/2);
                                        ci_art_all=net.param.ci_art_dc-net.param.ci_var_max(s)*sin(2*pi*freq(s)*(net.param.time+net.r_var.delay(s,l))-pi/2);
                                        r_art_all_arr=[r_art_all_arr; r_art_all];
                                        ri_art_all_arr=[ri_art_all_arr; ri_art_all];
                                        ci_art_all_arr=[ci_art_all_arr; ci_art_all];
                                        fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,r_art_all);    %delta Re
                                        fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,ri_art_all); %delta Ri
                                        fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,ci_art_all); %delta Ci
                                    end
                                end
                            else
                                if net.param.z_en==0
                                    fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,net.r_art_dc);
                                else
                                    fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,net.r_art_dc);  % Re
                                    fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,net.ri_art_dc); % Ri
                                    fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,net.ci_art_dc); % Ci
                                end
                            end
                        end
                    end
                end
                if found==0
                    if net.param.z_en==0
                        fprintf(fileID,'r_%g_%g %g %g %g\n',n1,n2,n1,n2,net.r_dc_arr(i,j,k).z);
                    else
                        fprintf(fileID,'r_%g_%g_e %g %g %g\n',n1,n2,n1,n2,net.r_dc_arr(i,j,k).z);  % Re
                        fprintf(fileID,'r_%g_%g_i %g_%g %g %g\n',n1,n2,n1,n2,n2,net.ri_dc_arr(i,j,k).z); % Ri
                        fprintf(fileID,'c_%g_%g_i %g %g_%g %g\n',n1,n2,n1,n1,n2,net.ci_dc_arr(i,j,k).z); % Ci
                    end
                end
            end
        end
    end
    net.r_var.r_art_all_arr=r_art_all_arr;
    net.r_var.ri_art_all_arr=ri_art_all_arr;
    net.r_var.ci_art_all_arr=ci_art_all_arr;
    % Current source
    for i=1:length(net.curr)
        if net.curr(i).ival>0
            curr_n1=net.curr(i).n1;
            curr_n2=net.curr(i).n2;
        else
            curr_n1=net.curr(i).n2;
            curr_n2=net.curr(i).n1;
        end
        fprintf(fileID,'i%g %g %g dc %g ac -%g 0\n',i,curr_n2(1),curr_n1(1),abs(net.curr(i).ival),(net.curr(i).iac));
        
        l=0;
        for j=1:size(net.curr(i).n2,2)-1
            l=l+1;
            fprintf(fileID,'v%g_%g %g %g dc 0 ac 0 0\n',i,l,curr_n2(1),curr_n2(j+1));
            l=l+1;
            fprintf(fileID,'v%g_%g %g %g dc 0 ac 0 0\n',i,l,curr_n1(1),curr_n1(j+1));
        end
    end
    % Sensing short circuit
    if ~isempty(net.param.s_n1i)
        for j=1:length(net.param.s_n1i)
            for m=1:size(net.s(j).n2,2)-1
                l=l+1;
                fprintf(fileID,'v%g_%g %g %g dc 0 ac 0 0\n',i,l,net.s(j).n2(1),net.s(j).n2(m+1));
                l=l+1;
                fprintf(fileID,'v%g_%g %g %g dc 0 ac 0 0\n',i,l,net.s(j).n1(1),net.s(j).n1(m+1));
            end
        end
    end
    
    % GND short circuit
    fprintf(fileID,'vgnd 0 %g dc 0 ac 0 0\n',curr_n2(1));
    if net.param.mode==0
        % Step Analysis
        if net.param.r_var_max>0
            fprintf(fileID,'.tran %g %g 0 %g\n',tstep,tend,tstep);            
        else
            fprintf(fileID,'.dc i1 %g %g 1\n',abs(net.curr(1).ival),abs(net.curr(1).ival));
        end
    else
        fprintf(fileID,'.ac lin %g %g %g\n',net.param.ac_freq(1),net.param.ac_freq(2),net.param.ac_freq(3));  %n_pts, fstart,fend
    end
    fprintf(fileID,'.backanno\n');
    fprintf(fileID,'.end\n');
    fclose(fileID);
    fprintf('... Netlist Generated')
else
    figure;
    sgtitle(strcat('Figure '))
    subplot(1,2,1)
    heatmap([0:net.Nj],[0:net.Ni],squeeze(net.node_arr(:,:,1)),'Colormap',white,'MissingDataLabel','Elec.');
    xlabel('X');ylabel('Y')
    title('Electrodes Location (Black square), Z=0')
    if net.Nk>0
        subplot(1,2,2)
        heatmap([0:net.Nj],[0:net.Nk],squeeze(net.node_arr(1,:,:))','Colormap',white,'MissingDataLabel','Artery');
        xlabel('X');ylabel('Z')
        title('Artery Location (Black square), Z=0')
    end
    set(gcf,'position',[0.1026    0.3058    1.3656    0.4624]*1e3)
end

