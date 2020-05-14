function FD_seis_sim(Vp,Vs,rho,X,Y,Z,topo,X_topo,Y_topo,Q0,f0,tsl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vp (required): 3D matrix with P velocity values
% Vs (required): 3D matrix with S velocity values
% rho (required): 3D matrix with density values
% X, Y, Z (required): 1D vector with velocity node locations
% topo (optional): 2D matrix with elevation values
% X_topo, Y_topo (necessary only with topo): 1D vector with topography node locations

% For the above input data units must be equivalent (e.g. km/s for Vp, Vs,
% Km for X, Y, Z, etc...). The Vp, Vs and rho matrices must have the same
% size. The size of X must be: length(X)=length(Vp(1,:,1)). Similarly, 
% length(Y)=length(Vp(:,1,1)) and length(Z)=length(Vp(1,1,:)). Topography
% node locations can be different from velocity, as they will be later
% interpolated into the velocity node locations. 
% similarly, length(X_topo)=length(topo(1,:)) and length(Y_topo)=length(Vp(:,1))

% Q0 (optional): 3D matrix with Q0 or TP values
% f0 (optional): 3D matrix with f0 or TS values
% tsl (optional): 3D matrix with tsl values
 
% !!! In order to use 3D attenuation, topography must also be used. If no
% topography data is available, dummy values for topo, X_topo, Y_topo can
% be used but they must be singular values, so they will be ignored in the
% calculations. Q0, f0 and/or tsl matrices must be the same size as Vp.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT DATA CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load model data and other parameters
% Check if source is properly defined
tic
if nargin<6
    msgbox('Incorrect number of input variables','input');
    return
end

%% Check if velocity model is properly defined
if isequal(ndims(Vp),3)==0 || isequal(ndims(Vs),3)==0 || isequal(ndims(rho),3)==0
    msgbox('Vp, Vs and/or density are not 3D matrices','input');
    return
elseif isequal(size(Vp),size(Vs),size(rho))==0
    msgbox('Vp, Vs and/or density matrices have different dimensions','input');
    return
elseif isequal(length(Vp(1,:,1)),length(X))==0
    msgbox('X node indices do not match velocity model nodes','input');
    return
elseif isequal(length(Vp(:,1,1)),length(Y))==0
    msgbox('Y node indices do not match velocity model nodes','input');
    return
elseif isequal(length(Vp(1,1,:)),length(Z))==0
    msgbox('Z node indices do not match velocity model nodes','input');
    return
else
    [XX_out,YY_out,ZZ_out]=meshgrid(X,Y,Z);
end

%% Check if topographic grid is given and fix model based on topography
topoload=0;
if nargin>6
    if isequal(length(topo(1,:)),length(X_topo))==0
        msgbox('X node indices do not match topography nodes','input');
        return
    elseif isequal(length(topo(:,1)),length(Y_topo))==0
        msgbox('Y node indices do not match topography nodes','input');
        return
    end
    topoload=1;
end

%% Check if 3D attenuation model is given (Q0/TP and/or f0/TS and/or tsl)
att_load=0;
if nargin==10
    if isequal(ndims(Q0),3)==0 
        msgbox('Q0/TP is not a 3D matrix','input');
        return
    elseif isequal(size(Vp),size(Q0))==0
        msgbox('Q0/TP matrix has different dimensions fom Vp model','input');
        return
    end   
att_load=1;    
end

if nargin==11
    if isequal(ndims(f0),3)==0 
        msgbox('f0/TS is not a 3D matrix','input');
        return
    elseif isequal(size(Vp),size(f0))==0
        msgbox('f0/TS matrix has different dimensions fom Vp model','input');
        return
    end   
att_load=2;    
end

if nargin==12
    if isequal(ndims(tsl),3)==0 
        msgbox('tsl is not a 3D matrix','input');
        return
    elseif isequal(size(Vp),size(tsl))==0
        msgbox('tsl matrix has different dimensions fom Vp model','input');
        return
    end   
att_load=3;    
end


%% Load the GUI
if topoload==0 && att_load==0
    FD_sim_gui(Vp,Vs,rho,XX_out,YY_out,ZZ_out,topoload,att_load)
elseif topoload==1 && att_load==0
    FD_sim_gui(Vp,Vs,rho,XX_out,YY_out,ZZ_out,topoload,att_load,topo,X_topo,Y_topo)
elseif topoload==1 && att_load==1
    FD_sim_gui(Vp,Vs,rho,XX_out,YY_out,ZZ_out,topoload,att_load,topo,X_topo,Y_topo,Q0)
elseif topoload==1 && att_load==2
    FD_sim_gui(Vp,Vs,rho,XX_out,YY_out,ZZ_out,topoload,att_load,topo,X_topo,Y_topo,Q0,f0)
elseif topoload==1 && att_load==3
    FD_sim_gui(Vp,Vs,rho,XX_out,YY_out,ZZ_out,topoload,att_load,topo,X_topo,Y_topo,Q0,f0,tsl)
elseif topoload==0 && att_load==1
    FD_sim_gui(Vp,Vs,rho,XX_out,YY_out,ZZ_out,topoload,att_load,Q0,f0)
elseif topoload==0 && att_load==2
    FD_sim_gui(Vp,Vs,rho,XX_out,YY_out,ZZ_out,topoload,att_load,Q0,f0)
elseif topoload==0 && att_load==3
    FD_sim_gui(Vp,Vs,rho,XX_out,YY_out,ZZ_out,topoload,att_load,Q0,f0,tsl)
end

