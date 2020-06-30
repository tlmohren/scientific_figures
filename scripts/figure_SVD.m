% -----------------------------------
% tlmohren 2020-06-03
% COMSOL with Euler-Lagrange data

%-------------------------------------
clc;clear all;close all
 
mydir  = pwd;
idcs   = strfind( mydir, filesep);
base_dir =  mydir(1:idcs(end)-1);
 
config_dir = fullfile( base_dir,'data', 'sim_standard_config.mat');
load_config = load( config_dir );

filenames = { fullfile( base_dir,'data', 'COMSOLwing_Yaw0_E3.0e+08.csv'), ...
     fullfile( base_dir,'data', 'ELwing_Yaw0_E3.0e+08.csv')}; 
  
Pars = load_config.Pars; 
Pars.plot_colors = [
    228,26,28
    55,126,184
    77,175,74]/255;  % insert colors from e.g. colorbrewer

%% load data 
SVD_data = struct;
header_lines = 10;

for i = 1:length(filenames)
    SVD_data(i).filename = filenames{i};
    
    X = csvread(SVD_data(i).filename, header_lines,0);  
    X_end = X(:,end-119:end);
     
    [SVD_data(i).U, SVD_data(i).S , SVD_data(i).V] = svd(X_end,'econ'); 
    SVD_data(i).t = Pars.dt*(0:length(SVD_data(i).V)-1);
    SVD_data(i).sigma = diag(SVD_data(i).S);
end

%% data plotting 

fig = figure();

width = 6.6;
height = 6; 

set(fig,'InvertHardcopy','on');
set(fig,'PaperUnits', 'inches');
set(fig,'units','inch');
papersize = get(fig, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];

set(fig, 'Position',  myfiguresize, 'PaperPosition', myfiguresize); 
 
Pars.n_modes = 5;
dy = 1/Pars.n_modes*0.95;

for i = 1:Pars.n_modes
    y_base = 0.07+dy*(i-1); 
    n_mode = Pars.n_modes -i+1;
    Pars.n_mode = n_mode;
    
    pos1 = [0.05, y_base,  0.13, dy-0.05]; 
    ax0 = subplot('Position',pos1); 
    hold on
    plot_singularvalues(ax0, SVD_data, Pars) 
    
    
    wing_height = (dy/2-0.025) -0.01;
    col2_width= 0.2;
    col2_start = 0.22;
    totaly = y_base+dy-0.05;
    
    pos1 = [col2_start, y_base,  col2_width, wing_height]; 
    ax0 = subplot('Position',pos1); 
    heatmap_wing(ax0, SVD_data(2), Pars)  
    
    pos1 = [col2_start, totaly-wing_height ,  col2_width, wing_height]; 
    ax0 = subplot('Position',pos1);
    heatmap_wing(ax0, SVD_data(1), Pars)  
    

    pos1 = [0.52, y_base,  0.22, dy-0.05]; 
    ax0 = subplot('Position',pos1); 
    hold on
    plot_time(ax0, SVD_data, Pars)  


    pos1 = [0.8, y_base,  0.18, dy-0.05];  
    ax0 = subplot('Position',pos1); 
    hold on
    plot_fft(ax0, SVD_data, Pars) 

end
 
%% Saving image
figname = fullfile( base_dir,'figs', 'figure_SVD');
print(fig, figname, '-dpng', '-r300');
print(fig, figname, '-dsvg', '-r300');
print(fig, figname, '-dpdf', '-r300');
 
%% Panel functions 

function plot_singularvalues(ax, SVD_data, Pars) 
    n = 6;
    for i = 1:length(SVD_data)
       singular_vals = SVD_data(i).sigma;
       pl = plot(ax, singular_vals,...
           '-', 'color',Pars.plot_colors(i,:) );
       sc = scatter(ax, 1:n, singular_vals(1:n),...
            10, Pars.plot_colors(i,:));
       sc.MarkerFaceAlpha = 0.5;
       sc.MarkerEdgeAlpha = 0.5;
       pl.Color(4) = 0.5;
       
       scatter(ax, Pars.n_mode, singular_vals(Pars.n_mode),...
            30, Pars.plot_colors(i,:),'filled')
    end 
    xlim([0,n+1])
    ylim([1e-7, 1])
    yticks([1e-6, 1e-3, 1])
    ax.YGrid = 'on';
    ax.XGrid = 'on';
    xticks([1:n])
    set(ax, 'YScale', 'log') 
    set(ax, 'color', 'none')
    if (Pars.n_mode == Pars.n_modes)
        tick_labels = repelem(" ",[6]);
        tick_labels(1) = "1";
        tick_labels(3) = "3"; 
        tick_labels(5) = "5"; 
        xticklabels(tick_labels);
        xlabel('Mode #')
    elseif (Pars.n_mode ==1)
       title('\Sigma')
    end
    if (Pars.n_mode ~= Pars.n_modes) 
        xticklabels(repelem(" ",[6])) 
    end
end


function heatmap_wing(ax, SVD_data, Pars)
    Un = reshape(SVD_data.U(:,Pars.n_mode) ,size(Pars.x));  
  
%     ph = pcolor(ax, Pars.x,Pars.y, Un );
%     set(ph, 'EdgeColor', 'none'); 

    imagesc(ax, Pars.x(:,1) ,Pars.y(1,:), Un );
 
    ax.Box= 'off';
   
    set(ax,'xticklabel',[])
    set(ax,'yticklabel',[]) 
    if Pars.n_mode==1
        if contains(SVD_data.filename,'COMSOL')
           th = title('U (COMSOL)');
            th.Color = Pars.plot_colors(1,:);
        else
           th = title('U (Euler)');
            th.Color = Pars.plot_colors(2,:);
        end
        titlePos = get( th , 'position');
        titlePos(2) = 0; 
        set( th , 'position' , titlePos);
    end
    axis(ax,'tight') 
    ax.YTick = [12]; 
    ax.XColor = 'none';
    ax.YAxis.LineWidth = 1.5;
end


function plot_time(ax, SVD_data, Pars)
    
    for i = 1:length(SVD_data)
          
        plot(ax, SVD_data(i).t, SVD_data(i).V(:,Pars.n_mode) ,...
            '-' , ...
            'color',Pars.plot_colors(i,:))
        set(ax, 'XLimSpec', 'Tight');
    end   
    set(ax, 'color', 'none')
    xticks(0:0.04:0.12)
    xlim([0,0.1201])
    ax.XGrid = 'on';
    
    if (Pars.n_mode == Pars.n_modes)
        xlabel('Time (s)')
    elseif (Pars.n_mode ==1)
       title('V')
    end
    if (Pars.n_mode ~= Pars.n_modes)
        xticklabels([])
    end
end


function [f, P1] = do_fft(time, signal)

    Y = fft(signal);
    L = length(signal);
    T = time(2)-time(1);
    Fs = 1/T;

    P2 = abs(Y/L);
    P1= P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L; 
    
end


function plot_fft(ax, SVD_data, Pars)
    
    for i = 1:length(SVD_data)
        
        time = SVD_data(i).t;
        signal = SVD_data(i).V(:,Pars.n_mode);
        [freq, Amp] = do_fft(time, signal);
        plot(ax, freq, Amp ,'-', 'color',Pars.plot_colors(i,:) )
        set(ax, 'XLimSpec', 'Tight');
        xlim([0, 110])
        ylim([0, 0.14])
        yticks([0,0.1])
        grid('on')
        xticks(0:25:100)
        if (Pars.n_mode == Pars.n_modes)
            xlabel('Frequency (Hz)')
        elseif (Pars.n_mode ==1)
           title('V (fft)')
        end
        if (Pars.n_mode ~= Pars.n_modes) 
            xticklabels([])
        end
    end   
    set(ax, 'color', 'none')
     
end