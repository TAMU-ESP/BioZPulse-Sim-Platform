clc
clear
close all
fclose('all');
set(0,'DefaultFigureWindowStyle','docked')

% General settings
write_netlist='Example';               % Name of results folder
write_figures=0;                       % Save Figures into results folder
write_mat=1;                           % Save Matlab Worksapce into .mat file results folder
complete_model=0;                      % 1=> Complete Wrist Model, 0=> Section of Wrist Model for faster simulation
debug=0;                               % Run/Skip SPICE Simulation 0=>Run Simulation, 1=> Skip simulation and plot only conductivity images for debugging purposes

% Electrode Configuration
ac_freq=[1 10e3 10e3 1];      % Frequency Analysis [number of freq. pts, fstart, fend, type] type:0=>linear, 1=>octave, 2=>decade
z_en=1;                       % Enable RC model (0=> R model, 1=> RC model)
se_en=1;                      % Enable Skin-electrode impedance RC model (0=> S.C., 1=> RC model)
icalc_en=1;                   % Enable current calculations

date_str=datestr(now,'yyyy-mm-dd_HH_MM');
for t1=1
    clear net;
    % Model parameters. Description in the paper
    if complete_model==1
        Lx=85;Ly=68;Lz=31;  
    else            
        Lx=30;Ly=68;Lz=11;   
    end
    Ex=14;Ey=7;
    L=2;    
    SVy=[10 10];
    SIy=SVy(1)*5;
    PIx=22;PIy=Ly*0.5-SVy(1)/2-(SIy-SVy(1))/2;
    PVx=[PIx PIx];PVy=PIy+[SVy(1),SVy(1)*3];
     
    % Whole Wrist Model (Define tissue type, center point (X,Y,Z) and the size for each dimension of the tissue cube)
    tissue_arr={'Fat','Muscle','Fat','Bone Cortical','Bone Cortical','Bone Cancellous','Bone Cancellous','Blood','Blood'};
    center_z_arr=[19,16,2,20,20,20,20,6,6]./Lz;
    center_x_arr=[43,43,43,26,68,26,68,22,66]./Lx;
    center_y_arr=[0.5,0.5,0.25,0.5,0.5,0.5,0.5,0.5,0.5];
    len_z_arr=[38,24,0,18,18,14,14,3,3]./Lz;
    len_x_arr=[85,77,0,38,25,34,21,3,2]./Lx;
    len_y_arr=[1,1,0,1,1,1,1,1,1];
    if complete_model==1
        art_en=[0,0,0,0,0,0,0,1,1];                  % Enable Radial and Ulnar arteries' model
        tilt_en=[0,0,0,0,0,0,0,1,1];                 % Enable tilt of Radial and Ulnar arteries
    else
        art_en=[0,0,0,0,0,0,0,1,0];                     % Enable Radial artery's model
        tilt_en=[0,0,0,0,0,0,0,1,0];                    % Enable tilt of Radial artery
    end    
    r_art_ratio=[0,0,0,0,0,0,0,1,3/4]*1.14e-2;        % Artery's impedance change for blood flow
    r_ratio=280;                                      % Scaling factor of tissue impedance at L=2 mm    
    se_ratio=18.86;                                   % Scaling factor of Skin impedance at L=2 mm
    
    net.param.write_netlist=write_netlist;
    net.param.write_figures=write_figures;
    net.param.write_mat=write_mat;
    net.param.debug=debug;
    net.param.date=date_str;
    net.param.path=strcat('./netlist/',net.param.write_netlist,'/');
    net.param.ac_freq=ac_freq;
    net.param.z_en=z_en;
    net.param.se_en=se_en;
    net.param.icalc_en=icalc_en;
    net.param.Lx=Lx;net.param.Ly=Ly;net.param.Lz=Lz;
    net.param.PIy=PIy;net.param.PIx=PIx;net.param.SIy=SIy;
    net.param.PVy=PVy;net.param.PVx=PVx;net.param.SVy=SVy;
    net.param.Ex=Ex;net.param.Ey=Ey;
    net.param.L=L;
    net.param.se_ratio=se_ratio;
    for i=1:length(center_x_arr)
        net.param.obj(i).tissue=tissue_arr(i);
        net.param.obj(i).r_ratio=r_ratio;
        net.param.obj(i).center.x=center_x_arr(i);
        net.param.obj(i).center.y=center_y_arr(i);
        net.param.obj(i).center.z=center_z_arr(i);
        net.param.obj(i).len.x=len_x_arr(i);
        net.param.obj(i).len.y=len_y_arr(i);
        net.param.obj(i).len.z=len_z_arr(i);
        net.param.obj(i).art_en=art_en(i);
        net.param.obj(i).tilt_en=tilt_en(i);
        net.param.obj(i).r_art_ratio=r_art_ratio(i);
    end
    
    % Reference
    %Create Grid
    net=img_gen_fn(net);
    net=grid_gen_fn(net);
    %Create electrodes
    [net]=pulse_elec_array_single_fn(net);
    % Simulate SPICE Netlist
    net=run_sim_fn(net);
    net_arr{t1}=net;   
    
end

% Save Data
if net(1).param.write_mat==1 && net(1).param.debug==0
    if ~exist(net(1).param.path,'dir')
        mkdir(net(1).param.path);
    end
    save(strcat(net(1).param.path,'netlist.mat'));
end


%%
set(0,'DefaultFigureWindowStyle','normal')
plot_img_arr_fn(net_arr,3,[1 55]*r_ratio,3); %0 -> impedance, 1-> conductivity, 2->conductivity+electrodes, 3->impedance+electrodes,colormap

if net(1).param.debug==0
    plot_spice_pulse_arr_fn(net_arr(1),[1],[1],2,1);
end

if net(1).param.debug==0
    for t1=1:length(net_arr)
        for i=1:length(net_arr{t1})
            net_arr{t1}(i)=spice_netlist_3d_pp_elec_fn(net_arr{t1}(i));
        end
    end
    plot_pulse_time_fn(net_arr{1},dsearchn(net_arr{1}(1).f_var(:,1),10e3)); % freq index
end




