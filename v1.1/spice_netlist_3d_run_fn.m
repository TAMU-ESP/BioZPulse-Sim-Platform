function net=spice_netlist_3d_run_fn(net)

%   Bassem Ibrahim
%   2019.08.22
%   BioZPulse: Bio-impedance Simulation Platform for Arterial Pulse Wave Modeling
%
%   Bassem Ibrahim, Drew A. Hall, Roozbeh Jafari, Bio-impedance Simulation
%   Platform using 3D Time-Varying Impedance Grid for Arterial Pulse Wave Modeling,
%   IEEE Biomedical Circuits and Systems Conference (BioCAS), October 17-19, 2019, Nara, Japan.
%
%   function net=spice_netlist_3d_run_fn(net)
%   Run LTSPICE Simulator

net.param.ltspice_path='D:\Program Files\LTC\LTspiceXVII\XVIIx64.exe';

fprintf('L=%g , Dimensions=(%2g,%2g,%2g)',net.param.L,net.Ni+1,net.Nj+1,net.Nk+1);
fprintf(' ,I1 loc.=(%2g,%2g) I2 loc.=(%2g,%2g), ZA=(%2g,%2g)\n',net.param.i_src_n1i,net.param.i_src_n1j,net.param.i_src_n2i,net.param.i_src_n2j,net.param.r_var_n1k,net.param.r_var_dia);

if ~isempty(net.param.s_n1i)
    fprintf('Elec. size=(%g,%g)\n',net.param.s_dia_x,net.param.s_dia_y);
    for i=1:length(net.param.s_n1i)
        fprintf('Sense Elec. 1 loc.=(%g,%g), Sense Elec. 2 loc.=(%g,%g)\n',net.param.s_n1i(i),net.param.s_n1j(i),net.param.s_n2i(i),net.param.s_n2j(i));
    end
end

res_folder=strcat('./results/',net.param.date,'_',net.param.write_netlist,'/');

if net.param.write_figures==1 || isempty(strmatch(net.param.write_netlist,'0'))
    net.param.path=res_folder;
    if ~exist(net.param.path,'dir')
        mkdir(net.param.path);
    end
else
    net.param.path='./';
end



net.param.circuit_file=strcat(net.param.path,'netlist.cir');
net.param.output_file=strcat(net.param.path,'netlist.raw');
net.param.output2_file=strcat(net.param.path,'netlist.op.raw');


if exist(net.param.output_file,'file')
    delete(net.param.output_file);
end
if exist(net.param.output2_file,'file')
    delete(net.param.output2_file);
end
if exist(net.param.circuit_file,'file')
    delete(net.param.circuit_file);
end


% Generate the SPICE Netlist
net=spice_netlist_3d_gen_fn(net);

if net.param.debug==0
    while ~exist(net.param.circuit_file,'file')
        pause(0.1);
    end
    fprintf('... Netlist found')
    status=1;
    while status~=0
        fprintf('... Spice Running')
        cmd=strcat({'"D:\Program Files\LTC\LTspiceXVII\XVIIx64.exe" -b -ascii '},net.param.circuit_file);
        status=system(cmd{1});
        pause(0.1);
    end
    if status==0
        fprintf('... Spice Finished')
    else
        fprintf('... Spice Failed')
    end
    d=0;
    repeat=0;
    while ~exist(net.param.output_file,'file')
        pause(0.1);
        d=d+1;
        if d==50
            repeat=1;
            break;
        end
    end
    
    if exist(net.param.output_file,'file')
        fprintf('... Output File Found')
    else
        fprintf('... Output File Not Found')
    end
    if repeat==1
        fprintf('... Repeat Spice')
    end
    if repeat==1
        while ~exist(net.param.circuit_file,'file')
            pause(0.1);
        end
        while ~exist(net.param.output_file,'file')
            fprintf('... Spice Running')
            cmd=strcat({'"D:\Program Files\LTC\LTspiceXVII\XVIIx64.exe" -b -ascii '},net.param.circuit_file);
            status=system(cmd{1});
            pause(0.1);
        end
    end
    
    if exist(net.param.output_file,'file')
        fprintf('... Output File Found')
    else
        fprintf('... Output File Not Found')
    end
end