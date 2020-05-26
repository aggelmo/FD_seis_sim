function w3=FD_sim_point_4_Q2(Vp,Vs,rho,X,Y,Z,source,pulse,TP,TS,tsl,Dt,maxt,posx,posy,posz,gpumode,N,b,topoload,topo,X_topo,Y_topo,savemode,savetype,skipt,outpath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%INPUT DATA%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vp: 3D P-velocity grid
% Vp: 3D S-velocity grid
% rho: 3D Pdensity grid
% X,Y,Z: x,y,z node locations (e.g. 0 100 200,... meters)
% source: vector of 5 elements: [strike dip rake mantissa exponent]
% pulse: N x 2 matrix. pulse(:,1)=time, pulse(:,2)= amplitude of source
% TP: P-waves relaxation times (Bohlen, 2002)
% TS: S-waves relaxation times (Bohlen, 2002)
% tsl: Stress relaxation times (Bohlen, 2002)
% Dt: Simulation sampliong period
% maxt: maximum simulation steps
% posx, posy, posz: x,y,z node index of source
% gpumode: 0-1 disable or enable gpu processing
% N: Number of nodes used as transition zone for ABC
% b:constant with values 0.3-0.5
% topoload: 0-1 flag stating that topography data is present
% topo: 2D elevation grid
% topo_X, topo_Y: x,y node locations of topographic grid
% savemode: 1=save whole 3D grids, 2=save only surface values
% savetype: 1=save as mat, else save as binary
% skipt: Save output per "skipt" timesteps
% outpath: Path where results will be stored
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

orpath=pwd; % Get the existing path
cd(outpath) % Change directory to save path
tic

%% Calculate node spacings
Dhx=single(zeros(size(Vp(1,:,1))));
Dhx(1:end-1)=abs(X(1,2:end,1)-X(1,1:end-1,1));
Dhx(end)=Dhx(end-1);
Dhy=single(zeros(size(Vp(:,1,1))));
Dhy(1:end-1)=abs(Y(2:end,1,1)-Y(1:end-1,1,1));
Dhy(end)=Dhy(end-1);
Dhz=single(zeros(size(Vp(1,1,:))));
Dhz(1:end-1)=abs(Z(1,1,2:end)-Z(1,1,1:end-1));
Dhz(end)=Dhz(end-1);

%% Prepare topography mask
if topoload==1   
    % Add nodes over the model if maximum topographic altitude is higher than maximum model altitude
    max_topo=max(topo(:));             % Get maximum topographic altitude
    max_model=max(Z(:));               % Get upper model boundary
    
    if max_topo>max_model              % Fix model if topographic altitude higher than model grid altitude.
        
        new_inc=Z(1,1,1)-Z(1,1,2);
        z_new=(max_model+new_inc:new_inc:max_topo+new_inc);     % Create new Z vector including the new nodes
        Z2=Z(1,1,:);
        Z_out=sort([z_new';Z2(:)],'descend');                % Sort the new vector
        [X,Y,Z]=meshgrid(X(1,:,1),Y(:,1,1),Z_out);           % Create new mesh

        % Initialize new Vp and Vs and density matrices
        Vp_new=zeros(size(Z));Vs_new=zeros(size(Z));rho_new=zeros(size(Z));         
        % Populate lower part of matrices
        Vp_new(:,:,length(z_new)+1:end)=Vp;Vs_new(:,:,length(z_new)+1:end)=Vs;rho_new(:,:,length(z_new)+1:end)=rho; 
        % Populate new nodes on the upper matrix part
        for i=1:length(z_new)
            Vp_new(:,:,i)=Vp(:,:,1);Vs_new(:,:,i)=Vs(:,:,1);rho_new(:,:,i)=rho(:,:,1);
        end    
        clear Vp Vs rho
        Vp=single(Vp_new);Vs=single(Vs_new);rho=single(rho_new);
        
        clear Vp_new Vs_new rho_new z_new new_inc Z2
    end
    
    % Reevaluate Z node spacings
    Dhz=single(zeros(size(Vp(1,1,:))));
    Dhz(1:end-1)=abs(Z(1,1,2:end)-Z(1,1,1:end-1));
    Dhz(end)=Dhz(end-1);
        
    disp('Gridding the topography into the velocity nodes')
    % Grid the topography into the normal stress nodes of the 3D model 
    
    w3=griddata(double(X_topo),double(Y_topo),double(topo),double(X(:,:,1)),...
        double(Y(:,:,1)),'linear');
    
    if max(isnan(w3(:)))==1
        msgbox('topography grid contains nan values after interpolation',...
            'set larger area for the topography in order to cover all nodes of the input model','input');
        return
    end
    w3=round(w3/Dhz(1))*Dhz(1);          % round towards dhz values
    % Now remove values from velocity model that fall over the topography.
    
    g_mask=single(zeros(size(Vp)));
    for i=1:length(Z(1,1,:))
        temp1=w3-(Z(:,:,i)+single(Dhz(1)/2)); 	% Create grid mask
        temp1(temp1<0)=0;                       % (free surface coincides with cell surface
        temp1(temp1>0)=1;                       % but normal stresses are always below)                    
        g_mask(:,:,i)=temp1;
    end
    
    % Find indices of normal stress nodes that lie over the topographic
    % surface
    g_mask_tzzi=single(find(g_mask==0));
    
    % Following part is to determine which nodes lie on vertical and horizontal surfaces of the 
    % free surface for Ohminato & Chouet surface BC
    % Interpolate grid on txy stress nodes   
    g_mask_txy=interp3(single(X),single(Y),...
        single(Z),single(g_mask),single(X)+single(Dhx)./2,...
        single(Y)+single(Dhy)./2,single(Z),'linear');
    g_mask_txy(abs(g_mask_txy)<1)=0;
    g_mask_txy(1,:,:)=g_mask_txy(2,:,:);g_mask_txy(:,1,:)=g_mask_txy(:,2,:);
    g_mask_txy(end,:,:)=g_mask_txy(end-1,:,:);g_mask_txy(:,end,:)=g_mask_txy(:,end-1,:);
    g_mask_txy(:,:,end)=g_mask_txy(:,:,end-1);
    g_mask_txyi=single(find(g_mask_txy==0));
    clear g_mask_txy
    
    % Interpolate grid on txz stress nodes   
    g_mask_txz=interp3(single(X),single(Y),...
        single(Z),single(g_mask),single(X)+single(Dhx)./2,single(Y),...
        single(Z)-single(Dhz)./2,'linear');
    g_mask_txz(g_mask_txz<1)=0;
    g_mask_txz(1,:,:)=g_mask_txz(2,:,:);g_mask_txz(:,1,:)=g_mask_txz(:,2,:);
    g_mask_txz(end,:,:)=g_mask_txz(end-1,:,:);g_mask_txz(:,end,:)=g_mask_txz(:,end-1,:);
    g_mask_txz(:,:,end)=g_mask_txz(:,:,end-1);
    g_mask_txzi=single(find(g_mask_txz==0));
    clear g_mask_txz
    
    % Interpolate grid on tyz stress nodes   
    g_mask_tyz=interp3(single(X),single(Y),...
        single(Z),single(g_mask),single(X),single(Y)+single(Dhy)./2,...
        single(Z)-single(Dhz)./2,'linear');
    g_mask_tyz(abs(g_mask_tyz)<1)=0;
    g_mask_tyz(1,:,:)=g_mask_tyz(2,:,:);g_mask_tyz(:,1,:)=g_mask_tyz(:,2,:);
    g_mask_tyz(end,:,:)=g_mask_tyz(end-1,:,:);g_mask_tyz(:,end,:)=g_mask_tyz(:,end-1,:);
    g_mask_tyz(:,:,end)=g_mask_tyz(:,:,end-1);
    g_mask_tyzi=single(find(g_mask_tyz==0));
    clear g_mask_tyz g_mask

    disp('Finished adjusting topography')
end

% remove unecessary matrices
clear X Y max_topo max_model X_topo Y_topo w3 temp1 X_out Y_out Z_out
%% Prepare source time function
a=interp1(pulse(:,1),pulse(:,2),0:Dt:(maxt*Dt)-Dt);  
a(isnan(a))=0;
a(abs(a)<0.0000001)=0;
a=a./max(abs(a)); % normalize time series

% Convert strike/dip/rake to 6 moment tensor components
[Gtxx,Gtxy,Gtxz,Gtyy,Gtyz,Gtzz]=sdr2mt(source(1)-90,source(2),source(3),source(4));

clear pulse
%% Find indices of surface nodes. Will be used only if output of surface values is requested
[ny,nx,nz]=size(Vp);    
if savemode==2   
    if topoload==1
        kk=1;
        topo_index=single(ones(nx*ny,1));
        for i=1:ny
            for j=1:nx
                [~,n1]=min(abs(Z(i,j,:)-topo(i,j)));
                topo_index(kk)=sub2ind(size(Vp),i,j,n1);
                kk=kk+1;
            end
        end
    end
end

clear n1 n2 Z kk i j topo
%% Matrix preallocation
% If GPU mode is selected, check if harware supports CUDA for GPU processing
if gpumode==1
try
    gpuDevice;
catch 
    disp('No GPU with CUDA found. Performing calculations on CPU')
    gpumode=0;
end
end

% Obtain Lame parameters from Vp, Vs and density
lamda=(rho.*Vp.^2)-(2.*rho.*Vs.^2);
mu=rho.*(Vs.^2); 

% remove unecessary matrices to free up memory
clear Vp Vs

disp('Preallocating matrices')    

if gpumode==1
    % Preallocate velocity component matrices
    Vx_new=gpuArray(single(zeros(ny,nx,nz)));   
    Vy_new=gpuArray(single(zeros(ny,nx,nz)));
    Vz_new=gpuArray(single(zeros(ny,nx,nz)));
    
    % Preallocate matrices for the application of boundary conditions
    Vz_old_y1=gpuArray(Vz_new(1:N+1,:,:));Vx_old_y1=gpuArray(Vx_new(1:N+1,:,:));Vy_old_y1=gpuArray(Vy_new(1:N+1,:,:)); %y=1 boundary
    Vz_old_x1=gpuArray(Vz_new(:,1:N+1,:));Vy_old_x1=gpuArray(Vy_new(:,1:N+1,:));Vx_old_x1=gpuArray(Vx_new(:,1:N+1,:)); %x=1 boundart
    Vz_old_yend=gpuArray(Vz_new(end-N:end,:,:));Vx_old_yend=gpuArray(Vx_new(end-N:end,:,:));Vy_old_yend=gpuArray(Vy_new(end-N:end,:,:)); %y=end boundary
    Vz_old_xend=gpuArray(Vz_new(:,end-N:end,:));Vx_old_xend=gpuArray(Vx_new(:,end-N:end,:));Vy_old_xend=gpuArray(Vy_new(:,end-N:end,:)); %x=end boundary
    Vz_old_zend=gpuArray(Vz_new(:,:,end-N:end));Vx_old_zend=gpuArray(Vx_new(:,:,end-N:end));Vy_old_zend=gpuArray(Vy_new(:,:,end-N:end)); %z=end boundary
    
    % Preallocate stress matrices
    tyy=gpuArray(single(zeros(ny,nx,nz)));txy=gpuArray(single(zeros(ny,nx,nz)));tyz=gpuArray(single(zeros(ny,nx,nz)));
    txx=gpuArray(single(zeros(ny,nx,nz)));tzz=gpuArray(single(zeros(ny,nx,nz)));txz=gpuArray(single(zeros(ny,nx,nz)));
    
    % Preallocate relaxation matrices
    ryy=gpuArray(single(zeros(ny,nx,nz)));rxy=gpuArray(single(zeros(ny,nx,nz)));ryz=gpuArray(single(zeros(ny,nx,nz)));
    rxx=gpuArray(single(zeros(ny,nx,nz)));rzz=gpuArray(single(zeros(ny,nx,nz)));rxz=gpuArray(single(zeros(ny,nx,nz)));

    % Convert velocity values to mechanical properties
    rho=gpuArray(rho);  %temporary!!!
    lamda=gpuArray(lamda);
    mu=gpuArray(mu);

    % Preallocate initial pulse for each mt component. Check scaling for
    % correct output units !!
    F_txx=gpuArray((10^source(5)).*Gtxx.*a);
    F_tyy=gpuArray((10^source(5)).*Gtyy.*a);
    F_tzz=gpuArray((10^source(5)).*Gtzz.*a);
    F_txy=gpuArray((10^source(5)).*Gtxy.*a);
    F_txz=gpuArray((10^source(5)).*Gtxz.*a);
    F_tyz=gpuArray((10^source(5)).*Gtyz.*a);

else % Do the same without GPU processing
    % Preallocate velocity component matrices
    Vx_new=single(zeros(ny,nx,nz));   
    Vy_new=single(zeros(ny,nx,nz));
    Vz_new=single(zeros(ny,nx,nz));

    % Preallocate matrices for the application of boundary conditions
    Vz_old_y1=Vz_new(1:N+1,:,:);Vx_old_y1=Vx_new(1:N+1,:,:);Vy_old_y1=Vy_new(1:N+1,:,:); %y=1 boundary
    Vz_old_x1=Vz_new(:,1:N+1,:);Vy_old_x1=Vy_new(:,1:N+1,:);Vx_old_x1=Vx_new(:,1:N+1,:); %x=1 boundart
    Vz_old_yend=Vz_new(end-N:end,:,:);Vx_old_yend=Vx_new(end-N:end,:,:);Vy_old_yend=Vy_new(end-N:end,:,:); %y=end boundary
    Vz_old_xend=Vz_new(:,end-N:end,:);Vx_old_xend=Vx_new(:,end-N:end,:);Vy_old_xend=Vy_new(:,end-N:end,:); %x=end boundary
    Vz_old_zend=Vz_new(:,:,end-N:end);Vx_old_zend=Vx_new(:,:,end-N:end);Vy_old_zend=Vy_new(:,:,end-N:end); %z=end boundary

    % Preallocate stress matrices
    tyy=single(zeros(ny,nx,nz));txy=single(zeros(ny,nx,nz));tyz=single(zeros(ny,nx,nz));
    txx=single(zeros(ny,nx,nz));tzz=single(zeros(ny,nx,nz));txz=single(zeros(ny,nx,nz));
   
    % Preallocate relaxation matrices
    ryy=(single(zeros(ny,nx,nz)));rxy=(single(zeros(ny,nx,nz)));ryz=(single(zeros(ny,nx,nz)));
    rxx=(single(zeros(ny,nx,nz)));rzz=(single(zeros(ny,nx,nz)));rxz=(single(zeros(ny,nx,nz)));

    % Preallocate initial pulse for each mt component. Check scaling for
    % correct output units !!
    F_txx=((10^source(5)).*Gtxx.*a);
    F_tyy=((10^source(5)).*Gtyy.*a);
    F_tzz=((10^source(5)).*Gtzz.*a);
    F_txy=((10^source(5)).*Gtxy.*a);
    F_txz=((10^source(5)).*Gtxz.*a);
    F_tyz=((10^source(5)).*Gtyz.*a);
end

% Check if 3D attenuation is given and set appropriate index of nodes for
% its application on the relaxation functions and all stress components. 
% This step is neccessary in order to avoid using a 3D matrix (or a separate script) if
% only singular TP/TS/tsl values are given.
if length(TP)==1
    indtii_TP1=1;indtii_TP2=1;indtii_TP3=1;
else
    indtii_TP1=3:ny-2;indtii_TP2=3:nx-2;indtii_TP3=3:nz-2;
end

if length(TS)==1
    indtii_TS1=1;indtii_TS2=1;indtii_TS3=1;
    indtxy_TS1=1;indtxy_TS2=1;indtxy_TS3=1;
    indtxz_TS1=1;indtxz_TS2=1;indtxz_TS3=1;
    indtyz_TS1=1;indtyz_TS2=1;indtyz_TS3=1;
else
    indtii_TS1=3:ny-2;indtii_TS2=3:nx-2;indtii_TS3=3:nz-2;
    indtxy_TS1=3:ny-1;indtxy_TS2=3:nx-1;indtxy_TS3=2:nz-2;
    indtxz_TS1=2:ny-2;indtxz_TS2=3:nx-1;indtxz_TS3=3:nz-1;
    indtyz_TS1=3:ny-1;indtyz_TS2=2:nx-2;indtyz_TS3=3:nz-1;
end

if length(tsl)==1
    indtii_tsl1=1;indtii_tsl2=1;indtii_tsl3=1;
    indtxy_tsl1=1;indtxy_tsl2=1;indtxy_tsl3=1;
    indtxz_tsl1=1;indtxz_tsl2=1;indtxz_tsl3=1;
    indtyz_tsl1=1;indtyz_tsl2=1;indtyz_tsl3=1;
else
    indtii_tsl1=3:ny-2;indtii_tsl2=3:nx-2;indtii_tsl3=3:nz-2;
    indtxy_tsl1=3:ny-1;indtxy_tsl2=3:nx-1;indtxy_tsl3=2:nz-2;
    indtxz_tsl1=2:ny-2;indtxz_tsl2=3:nx-1;indtxz_tsl3=3:nz-1;
    indtyz_tsl1=3:ny-1;indtyz_tsl2=2:nx-2;indtyz_tsl3=3:nz-1;
end
    

dura=toc;
%% Main calculation section
tic
kk=1;tt=1;
a=9/8;b2=-1/24; % 4th order coefficients

% remove unecessary matrices
clear XX_out YY_out ZZ_out Vp Vs topo X_topo Y_topo

% Loop over time
for t=1:maxt
    
    disp([num2str(100*t/maxt) '%']) % Show completion percentage
    
    % Calculate X velocity component for whole grid 
    Vx_new(2:end-2,2:end-2,2:end-2)=Vx_new(2:end-2,2:end-2,2:end-2)+...
        ((Dt.*(1./(rho(2:end-2,2:end-2,2:end-2)))).*...
        ((a.*(txx(2:end-2,3:end-1,2:end-2)-txx(2:end-2,2:end-2,2:end-2))+...
        b2.*(txx(2:end-2,4:end,2:end-2)-txx(2:end-2,1:end-3,2:end-2)))./Dhx(2:end-2)+...
        (a.*(txy(3:end-1,3:end-1,2:end-2)-txy(2:end-2,3:end-1,2:end-2))+...
        b2.*(txy(4:end,3:end-1,2:end-2)-txy(1:end-3,3:end-1,2:end-2)))./Dhy(2:end-2)+...
        (a.*(txz(2:end-2,3:end-1,3:end-1)-txz(2:end-2,3:end-1,2:end-2))+...
        b2.*(txz(2:end-2,3:end-1,4:end)-txz(2:end-2,3:end-1,1:end-3)))./Dhz(2:end-2)));


    % Insert source at Vx
    Vx_new(posy,posx-1,posz)=Vx_new(posy,posx-1,posz)+(Dt./(rho(posy,posx-1,posz))).*(F_txx(t)./((Dhx(posx-1).^2.*Dhy(posy).*Dhz(posz))));%
    Vx_new(posy,posx,posz)=Vx_new(posy,posx,posz)+(Dt./(rho(posy,posx,posz))).*(F_txx(t)./((Dhx(posx).^2.*Dhy(posy).*Dhz(posz))));%
    Vx_new(posy+1,posx-1,posz)=Vx_new(posy+1,posx-1,posz)+(Dt./(rho(posy+1,posx-1,posz))).*(F_txy(t)./((Dhx(posx-1).*4.*Dhy(posy+1).^2.*Dhz(posz))));% 
    Vx_new(posy-1,posx-1,posz)=Vx_new(posy-1,posx-1,posz)-(Dt./(rho(posy-1,posx-1,posz))).*(F_txy(t)./((Dhx(posx-1).*4.*Dhy(posy-1).^2.*Dhz(posz))));% 
    Vx_new(posy+1,posx,posz)=Vx_new(posy+1,posx,posz)+(Dt./(rho(posy+1,posx,posz))).*(F_txy(t)./((Dhx(posx).*4.*Dhy(posy+1).^2.*Dhz(posz))));% 
    Vx_new(posy-1,posx,posz)=Vx_new(posy-1,posx,posz)+(Dt./(rho(posy-1,posx,posz))).*(F_txy(t)./((Dhx(posx).*4.*Dhy(posy-1).^2.*Dhz(posz))));% 
    Vx_new(posy,posx-1,posz+1)=Vx_new(posy,posx-1,posz+1)+(Dt./(rho(posy,posx-1,posz+1))).*(F_txz(t)./((Dhx(posx-1).*Dhy(posy).*4.*Dhz(posz+1).^2)));% 
    Vx_new(posy,posx-1,posz-1)=Vx_new(posy,posx-1,posz-1)+(Dt./(rho(posy,posx-1,posz-1))).*(F_txz(t)./((Dhx(posx-1).*Dhy(posy).*4.*Dhz(posz-1).^2)));% 
    Vx_new(posy,posx,posz+1)=Vx_new(posy,posx,posz+1)+(Dt./(rho(posy,posx,posz+1))).*(F_txz(t)./((Dhx(posx).*Dhy(posy).*4.*Dhz(posz+1).^2)));% 
    Vx_new(posy,posx,posz-1)=Vx_new(posy,posx,posz-1)+(Dt./(rho(posy,posx,posz-1))).*(F_txz(t)./((Dhx(posx).*Dhy(posy).*4.*Dhz(posz-1).^2)));% 

    % Calculate Y velocity component for whole grid  
    Vy_new(2:end-2,2:end-2,2:end-2)=Vy_new(2:end-2,2:end-2,2:end-2)+...
        ((Dt.*(1./(rho(2:end-2,2:end-2,2:end-2)))).*...
        ((a.*(tyy(3:end-1,2:end-2,2:end-2)-tyy(2:end-2,2:end-2,2:end-2))+...
        b2.*(tyy(4:end,2:end-2,2:end-2)-tyy(1:end-3,2:end-2,2:end-2)))./Dhy(2:end-2)+...
        (a.*(txy(3:end-1,3:end-1,2:end-2)-txy(3:end-1,2:end-2,2:end-2))+...
        b2.*(txy(3:end-1,4:end,2:end-2)-txy(3:end-1,1:end-3,2:end-2)))./Dhx(2:end-2)+...
        (a.*(tyz(3:end-1,2:end-2,3:end-1)-tyz(3:end-1,2:end-2,2:end-2))+...
        b2.*(tyz(3:end-1,2:end-2,4:end)-tyz(3:end-1,2:end-2,1:end-3)))./Dhz(2:end-2)));

    % Insert source at Vx
    Vy_new(posy-1,posx,posz)=Vy_new(posy-1,posx,posz)+(Dt./(rho(posy-1,posx,posz))).*(F_tyy(t)./((Dhx(posx).*Dhy(posy).^2.*Dhz(posz))));% 
    Vy_new(posy,posx,posz)=Vy_new(posy,posx,posz)+(Dt./(rho(posy,posx,posz))).*(F_tyy(t)./((Dhx(posx).*Dhy(posy).^2.*Dhz(posz))));% 
    Vy_new(posy-1,posx+1,posz)=Vy_new(posy-1,posx+1,posz)+(Dt./(rho(posy-1,posx+1,posz))).*(F_txy(t)./((4.*Dhx(posx+1).^2.*Dhy(posy-1).*Dhz(posz))));% 
    Vy_new(posy-1,posx-1,posz)=Vy_new(posy-1,posx-1,posz)+(Dt./(rho(posy-1,posx-1,posz))).*(F_txy(t)./((4.*Dhx(posx-1).^2.*Dhy(posy-1).*Dhz(posz))));% 
    Vy_new(posy,posx+1,posz)=Vy_new(posy,posx+1,posz)+(Dt./(rho(posy,posx+1,posz))).*(F_txy(t)./((4.*Dhx(posx+1).^2.*Dhy(posy).*Dhz(posz))));% 
    Vy_new(posy,posx-1,posz)=Vy_new(posy,posx-1,posz)+(Dt./(rho(posy,posx-1,posz))).*(F_txy(t)./((4.*Dhx(posx-1).^2.*Dhy(posy).*Dhz(posz))));% 
    Vy_new(posy-1,posx,posz+1)=Vy_new(posy-1,posx,posz+1)+(Dt./(rho(posy-1,posx,posz+1))).*(F_tyz(t)./((Dhx(posx).*Dhy(posy-1).*4.*Dhz(posz+1).^2)));% 
    Vy_new(posy-1,posx,posz-1)=Vy_new(posy-1,posx,posz-1)+(Dt./(rho(posy-1,posx,posz-1))).*(F_tyz(t)./((Dhx(posx).*Dhy(posy-1).*4.*Dhz(posz-1).^2)));% 
    Vy_new(posy,posx,posz+1)=Vy_new(posy,posx,posz+1)+(Dt./(rho(posy,posx,posz+1))).*(F_tyz(t)./((Dhx(posx).*Dhy(posy).*4.*Dhz(posz+1).^2)));% 
    Vy_new(posy,posx,posz-1)=Vy_new(posy,posx,posz-1)+(Dt./(rho(posy,posx,posz-1))).*(F_tyz(t)./((Dhx(posx).*Dhy(posy).*4.*Dhz(posz-1).^2)));% 
    
    % Calculate Z velocity component for whole grid  
    Vz_new(2:end-2,2:end-2,2:end-2)=Vz_new(2:end-2,2:end-2,2:end-2)+...
        ((Dt.*(1./(rho(2:end-2,2:end-2,2:end-2)))).*...
        ((a.*(tzz(2:end-2,2:end-2,3:end-1)-tzz(2:end-2,2:end-2,2:end-2))+...
        b2.*(tzz(2:end-2,2:end-2,4:end)-tzz(2:end-2,2:end-2,3:end-1)))./Dhz(2:end-2)+...
        (a.*(txz(2:end-2,3:end-1,3:end-1)-txz(2:end-2,2:end-2,3:end-1))+...
        b2.*(txz(2:end-2,4:end,3:end-1)-txz(2:end-2,1:end-3,3:end-1)))./Dhx(2:end-2)+...
        (a.*(tyz(3:end-1,2:end-2,3:end-1)-tyz(2:end-2,2:end-2,3:end-1))+...
        b2.*(tyz(4:end,2:end-2,3:end-1)-tyz(3:end-1,2:end-2,3:end-1)))./Dhy(2:end-2)));


    % Insert source at Vz
    Vz_new(posy,posx,posz-1)=Vz_new(posy,posx,posz-1)+(Dt./(rho(posy,posx,posz-1))).*(F_tzz(t)./((Dhx(posx).*Dhy(posy).*Dhz(posz).^2)));% 
    Vz_new(posy,posx,posz)=Vz_new(posy,posx,posz)+(Dt./(rho(posy,posx,posz))).*(F_tzz(t)./((Dhx(posx).*Dhy(posy).*Dhz(posz).^2)));% 
    Vz_new(posy,posx+1,posz-1)=Vz_new(posy,posx+1,posz-1)+(Dt./(rho(posy,posx+1,posz-1))).*(F_txz(t)./((Dhx(posx+1).^2.*4.*Dhy(posy).*Dhz(posz-1))));% 
    Vz_new(posy,posx-1,posz-1)=Vz_new(posy,posx-1,posz-1)+(Dt./(rho(posy,posx-1,posz-1))).*(F_txz(t)./((Dhx(posx-1).^2.*4.*Dhy(posy).*Dhz(posz-1))));% 
    Vz_new(posy,posx+1,posz)=Vz_new(posy,posx+1,posz)+(Dt./(rho(posy,posx+1,posz))).*(F_txz(t)./((Dhx(posx+1).^2.*4.*Dhy(posy).*Dhz(posz))));% 
    Vz_new(posy,posx-1,posz)=Vz_new(posy,posx-1,posz)+(Dt./(rho(posy,posx-1,posz))).*(F_txz(t)./((Dhx(posx-1).^2.*4.*Dhy(posy).*Dhz(posz))));% 
    Vz_new(posy+1,posx,posz-1)=Vz_new(posy+1,posx,posz-1)+(Dt./(rho(posy+1,posx,posz-1))).*(F_tyz(t)./((4.*Dhx(posx).*Dhy(posy+1).^2.*Dhz(posz-1))));% 
    Vz_new(posy-1,posx,posz-1)=Vz_new(posy-1,posx,posz-1)+(Dt./(rho(posy-1,posx,posz-1))).*(F_tyz(t)./((4.*Dhx(posx).*Dhy(posy-1).^2.*Dhz(posz-1))));% 
    Vz_new(posy+1,posx,posz)=Vz_new(posy+1,posx,posz)+(Dt./(rho(posy+1,posx,posz))).*(F_tyz(t)./((4.*Dhx(posx).*Dhy(posy+1).^2.*Dhz(posz))));% 
    Vz_new(posy-1,posx,posz)=Vz_new(posy-1,posx,posz)+(Dt./(rho(posy-1,posx,posz))).*(F_tyz(t)./((4.*Dhx(posx).*Dhy(posy-1).^2.*Dhz(posz))));% 
  
    
    % Higdon ABC. Vp and Vs values of normal stress nodes are used!! (no apparent difference if otherwise...)
    wb=((1:N)-1)/N;
    % Boundary y=1
    beta=(1+(sqrt((lamda(1:N,:,:)+2.*mu(1:N,:,:))./rho(1:N,:,:))./sqrt(mu(1:N,:,:)./rho(1:N,:,:))))./2;
    r=(sqrt(mu(1:N,:,:)./rho(1:N,:,:)).*Dt)./Dhy(1:N);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vz_new(1:N,:,:)=(wb'.*Vz_new(1:N,:,:))+((1-wb').*(-qx.*Vz_new(2:N+1,:,:)-qt.*Vz_old_y1(1:N,:,:)-qxt.*Vz_old_y1(2:N+1,:,:)));     
    Vx_new(1:N,:,:)=(wb'.*Vx_new(1:N,:,:))+((1-wb').*(-qx.*Vx_new(2:N+1,:,:)-qt.*Vx_old_y1(1:N,:,:)-qxt.*Vx_old_y1(2:N+1,:,:))); 
    r=(sqrt((lamda(1:N,:,:)+2.*mu(1:N,:,:))./rho(1:N,:,:)).*Dt)./Dhy(1:N);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vy_new(1:N,:,:)=(wb'.*Vy_new(1:N,:,:))+((1-wb').*(-qx.*Vy_new(2:N+1,:,:)-qt.*Vy_old_y1(1:N,:,:)-qxt.*Vy_old_y1(2:N+1,:,:)));     
    clear beta r qx qt
    % Boundary x=1
    beta=(1+(sqrt((lamda(:,1:N,:)+2.*mu(:,1:N,:))./rho(:,1:N,:))./sqrt(mu(:,1:N,:)./rho(:,1:N,:))))./2;
    r=(sqrt(mu(:,1:N,:)./rho(:,1:N,:)).*Dt)./Dhx(1:N);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vz_new(:,1:N,:)=(wb.*Vz_new(:,1:N,:))+((1-wb).*(-qx.*Vz_new(:,2:N+1,:)-qt.*Vz_old_x1(:,1:N,:)-qxt.*Vz_old_x1(:,2:N+1,:)));   
    Vy_new(:,1:N,:)=(wb.*Vy_new(:,1:N,:))+((1-wb).*(-qx.*Vy_new(:,2:N+1,:)-qt.*Vy_old_x1(:,1:N,:)-qxt.*Vy_old_x1(:,2:N+1,:)));    
    r=(sqrt((lamda(:,1:N,:)+2.*mu(:,1:N,:))./rho(:,1:N,:)).*Dt)./Dhx(1:N);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vx_new(:,1:N,:)=(wb.*Vx_new(:,1:N,:))+((1-wb).*(-qx.*Vx_new(:,2:N+1,:)-qt.*Vx_old_x1(:,1:N,:)-qxt.*Vx_old_x1(:,2:N+1,:)));    
    clear beta r qx qt
    % Boundary y=end
    beta=((1+(sqrt((lamda(end-(N-1):end,:,:)+2.*mu(end-(N-1):end,:,:))./rho(end-(N-1):end,:,:))./sqrt(mu(end-(N-1):end,:,:)./rho(end-(N-1):end,:,:))))./2);
    r=(sqrt(mu(end-(N-1):end,:,:)./rho(end-(N-1):end,:,:)).*Dt)./Dhy(end-(N-1):end);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vz_new(end-(N-1):end,:,:)=((1-wb').*Vz_new(end-(N-1):end,:,:))+(wb'.*(-qx.*Vz_new(end-N:end-1,:,:)-qt.*Vz_old_yend(end-(N-1):end,:,:)-qxt.*Vz_old_yend(end-N:end-1,:,:)));     
    Vx_new(end-(N-1):end,:,:)=((1-wb').*Vx_new(end-(N-1):end,:,:))+(wb'.*(-qx.*Vx_new(end-N:end-1,:,:)-qt.*Vx_old_yend(end-(N-1):end,:,:)-qxt.*Vx_old_yend(end-N:end-1,:,:)));       
    r=(sqrt((lamda(end-(N-1):end,:,:)+2.*mu(end-(N-1):end,:,:))./rho(end-(N-1):end,:,:)).*Dt)./Dhy(end-(N-1):end);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vy_new(end-(N-1):end,:,:)=((1-wb').*Vy_new(end-(N-1):end,:,:))+(wb'.*(-qx.*Vy_new(end-N:end-1,:,:)-qt.*Vy_old_yend(end-(N-1):end,:,:)-qxt.*Vy_old_yend(end-N:end-1,:,:)));  
    clear beta r qx qt
    % Boundary x=end
    beta=(1+(sqrt((lamda(:,end-(N-1):end,:)+2.*mu(:,end-(N-1):end,:))./rho(:,end-(N-1):end,:))./sqrt(mu(:,end-(N-1):end,:)./rho(:,end-(N-1):end,:))))./2;
    r=(sqrt(mu(:,end-(N-1):end,:)./rho(:,end-(N-1):end,:)).*Dt)./Dhx(end-(N-1):end);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vz_new(:,end-(N-1):end,:)=((1-wb).*Vz_new(:,end-(N-1):end,:))+(wb.*(-qx.*Vz_new(:,end-N:end-1,:)-qt.*Vz_old_xend(:,end-(N-1):end,:)-qxt.*Vz_old_xend(:,end-N:end-1,:)));    
    Vy_new(:,end-(N-1):end,:)=((1-wb).*Vy_new(:,end-(N-1):end,:))+(wb.*(-qx.*Vy_new(:,end-N:end-1,:)-qt.*Vy_old_xend(:,end-(N-1):end,:)-qxt.*Vy_old_xend(:,end-N:end-1,:)));   
    r=(sqrt((lamda(:,end-(N-1):end,:)+2.*mu(:,end-(N-1):end,:))./rho(:,end-(N-1):end,:)).*Dt)./Dhx(end-(N-1):end);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vx_new(:,end-(N-1):end,:)=((1-wb).*Vx_new(:,end-(N-1):end,:))+(wb.*(-qx.*Vx_new(:,end-N:end-1,:)-qt.*Vx_old_xend(:,end-(N-1):end,:)-qxt.*Vx_old_xend(:,end-N:end-1,:)));       
    clear beta r qx qt
    % Boundary z=end
    wb2(1,1,:)=wb;
    beta=(1+(sqrt((lamda(:,:,end-(N-1):end)+2.*mu(:,:,end-(N-1):end))./rho(:,:,end-(N-1):end))./sqrt(mu(:,:,end-(N-1):end)./rho(:,:,end-(N-1):end))))./2;
    r=(sqrt((lamda(:,:,end-(N-1):end)+2.*mu(:,:,end-(N-1):end))./rho(:,:,end-(N-1):end)).*Dt)./Dhz(end-(N-1):end);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vz_new(:,:,end-(N-1):end)=((1-wb2).*Vz_new(:,:,end-(N-1):end))+(wb2.*(-qx.*Vz_new(:,:,end-N:end-1)-qt.*Vz_old_zend(:,:,end-(N-1):end)-qxt.*Vz_old_zend(:,:,end-N:end-1)));       
    r=(sqrt(mu(:,:,end-(N-1):end)./rho(:,:,end-(N-1):end)).*Dt)./Dhz(end-(N-1):end);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vx_new(:,:,end-(N-1):end)=((1-wb2).*Vx_new(:,:,end-(N-1):end))+(wb2.*(-qx.*Vx_new(:,:,end-N:end-1)-qt.*Vx_old_zend(:,:,end-(N-1):end)-qxt.*Vx_old_zend(:,:,end-N:end-1)));     
    Vy_new(:,:,end-(N-1):end)=((1-wb2).*Vy_new(:,:,end-(N-1):end))+(wb2.*(-qx.*Vy_new(:,:,end-N:end-1)-qt.*Vy_old_zend(:,:,end-(N-1):end)-qxt.*Vy_old_zend(:,:,end-N:end-1)));
    clear beta r qx qt
    
    % Keep old velocity component values for the calculation of ABC on the
    % next time step
    Vz_old_y1=Vz_new(1:N+1,:,:);Vx_old_y1=Vx_new(1:N+1,:,:);Vy_old_y1=Vy_new(1:N+1,:,:); %y=1 boundary
    Vz_old_x1=Vz_new(:,1:N+1,:);Vy_old_x1=Vy_new(:,1:N+1,:);Vx_old_x1=Vx_new(:,1:N+1,:); %x=1 boundart
    Vz_old_yend=Vz_new(end-N:end,:,:);Vx_old_yend=Vx_new(end-N:end,:,:);Vy_old_yend=Vy_new(end-N:end,:,:); %y=end boundary
    Vz_old_xend=Vz_new(:,end-N:end,:);Vx_old_xend=Vx_new(:,end-N:end,:);Vy_old_xend=Vy_new(:,end-N:end,:); %x=end boundary
    Vz_old_zend=Vz_new(:,:,end-N:end);Vx_old_zend=Vx_new(:,:,end-N:end);Vy_old_zend=Vy_new(:,:,end-N:end); %z=end boundary
  

    % Create ouput based on user choice
    if kk==skipt
        if savemode==1
            if savetype==1
                % Save whole 3D as .mat files
                wfield=single(gather(Vx_new));
                save(['Vx_' num2str(tt)],'wfield')
                wfield=single(gather(Vy_new));
                save(['Vy_' num2str(tt)],'wfield')
                wfield=single(gather(Vz_new));
                save(['Vz_' num2str(tt)],'wfield')
                clear wfield
                    
            else
                % Save whole 3D as binary files
                fileID1 = fopen(['Vx_' num2str(tt) '.bin'],'w');
                fileID2 = fopen(['Vy_' num2str(tt) '.bin'],'w');
                fileID3 = fopen(['Vz_' num2str(tt) '.bin'],'w');

                fwrite(fileID1,single(gather(Vx_new)),'float');
                fwrite(fileID2,single(gather(Vy_new)),'float');
                fwrite(fileID3,single(gather(Vz_new)),'float');
                fclose(fileID1);
                fclose(fileID2);
                fclose(fileID3);
            end
            tt=tt+1;
        elseif savemode==2 && topoload==1
            if savetype==1
                % Save only surface values as .mat files
                wfield=gather(Vx_new(topo_index));wfield=single(reshape(wfield,nx,ny))';
                save(['Vx_' num2str(tt)],'wfield')
                wfield=gather(Vy_new(topo_index));wfield=single(reshape(wfield,nx,ny))';
                save(['Vy_' num2str(tt)],'wfield')
                wfield=gather(Vz_new(topo_index));wfield=single(reshape(wfield,nx,ny))';
                save(['Vz_' num2str(tt)],'wfield')
                clear wfield
            else
                % Save only surface values as binary files
                fileID1 = fopen(['Vx_' num2str(tt) '.bin'],'w');
                fileID2 = fopen(['Vy_' num2str(tt) '.bin'],'w');
                fileID3 = fopen(['Vz_' num2str(tt) '.bin'],'w');

                wfield=Vx_new(topo_index);wfield=reshape(wfield,nx,ny);
                fwrite(fileID1,single(gather(wfield))','float');
                wfield=Vy_new(topo_index);wfield=reshape(wfield,nx,ny);
                fwrite(fileID2,single(gather(wfield))','float');
                wfield=Vz_new(topo_index);wfield=reshape(wfield,nx,ny);
                fwrite(fileID3,single(gather(wfield))','float');
                fclose(fileID1);fclose(fileID2);fclose(fileID3);
                clear wfield
            end
            tt=tt+1;
        else
            if savetype==1
                % Save only surface values as .mat files (for flat
                % topographic surface)
                wfield=single(gather(Vx_new(:,:,2)));
                save(['Vx_' num2str(tt)],'wfield')
                wfield=single(gather(Vy_new(:,:,2)));
                save(['Vy_' num2str(tt)],'wfield')
                wfield=single(gather(Vz_new(:,:,2)));
                save(['Vz_' num2str(tt)],'wfield')
                clear wfield
                    
            else
                % Save only surface values as binary files (for flat
                % topographic surface)
                fileID1 = fopen(['Vx_' num2str(tt) '.bin'],'w');
                fileID2 = fopen(['Vy_' num2str(tt) '.bin'],'w');
                fileID3 = fopen(['Vz_' num2str(tt) '.bin'],'w');

                fwrite(fileID1,single(gather(Vx_new(:,:,2))),'float');
                fwrite(fileID2,single(gather(Vy_new(:,:,2))),'float');
                fwrite(fileID3,single(gather(Vz_new(:,:,2))),'float');
                fclose(fileID1);
                fclose(fileID2);
                fclose(fileID3);
            end
            tt=tt+1;
        end
    end

    % Calculate rxx relaation funciton
    rxx(3:end-2,3:end-2,3:end-2)=((1+(Dt./2.*tsl(indtii_tsl1,indtii_tsl2,indtii_tsl3))).^-1).*...
        (((1-(Dt./2.*tsl(indtii_tsl1,indtii_tsl2,indtii_tsl3))).*rxx(3:end-2,3:end-2,3:end-2))-...
        ((lamda(3:end-2,3:end-2,3:end-2).*Dt./tsl(indtii_tsl1,indtii_tsl2,indtii_tsl3)).*...
        TP(indtii_TP1,indtii_TP2,indtii_TP3).*...
        (((a.*(Vx_new(3:end-2,3:end-2,3:end-2)-Vx_new(3:end-2,2:end-3,3:end-2))+...
        b2.*(Vx_new(3:end-2,4:end-1,3:end-2)-Vx_new(3:end-2,1:end-4,3:end-2)))./Dhx(3:end-2))+...
        ((a.*(Vy_new(3:end-2,3:end-2,3:end-2)-Vy_new(2:end-3,3:end-2,3:end-2))+...
        b2.*(Vy_new(4:end-1,3:end-2,3:end-2)-Vy_new(1:end-4,3:end-2,3:end-2)))./Dhy(3:end-2))+...
        ((a.*(Vz_new(3:end-2,3:end-2,3:end-2)-Vz_new(3:end-2,3:end-2,2:end-3))+...
        b2.*(Vz_new(3:end-2,3:end-2,4:end-1)-Vz_new(3:end-2,3:end-2,1:end-4)))./Dhz(3:end-2))))+...
        (2.*Dt.*mu(3:end-2,3:end-2,3:end-2)./tsl(indtii_tsl1,indtii_tsl2,indtii_tsl3)).*...
        TS(indtii_TS1,indtii_TS2,indtii_TS3).*...
        (((a.*(Vy_new(3:end-2,3:end-2,3:end-2)-Vy_new(2:end-3,3:end-2,3:end-2))+...
        b2.*(Vy_new(4:end-1,3:end-2,3:end-2)-Vy_new(1:end-4,3:end-2,3:end-2)))./Dhy(3:end-2))+...
        ((a.*(Vz_new(3:end-2,3:end-2,3:end-2)-Vz_new(3:end-2,3:end-2,2:end-3))+...
        b2.*(Vz_new(3:end-2,3:end-2,4:end-1)-Vz_new(3:end-2,3:end-2,1:end-4)))./Dhz(3:end-2))));

    % Calculate ryy relaation funciton
    ryy(3:end-2,3:end-2,3:end-2)=((1+(Dt./2.*tsl(indtii_tsl1,indtii_tsl2,indtii_tsl3))).^-1).*...
        (((1-(Dt./2.*tsl(indtii_tsl1,indtii_tsl2,indtii_tsl3))).*ryy(3:end-2,3:end-2,3:end-2))-...
        ((lamda(3:end-2,3:end-2,3:end-2).*Dt./tsl(indtii_tsl1,indtii_tsl2,indtii_tsl3)).*...
        TP(indtii_TP1,indtii_TP2,indtii_TP3).*...
        (((a.*(Vx_new(3:end-2,3:end-2,3:end-2)-Vx_new(3:end-2,2:end-3,3:end-2))+...
        b2.*(Vx_new(3:end-2,4:end-1,3:end-2)-Vx_new(3:end-2,1:end-4,3:end-2)))./Dhx(3:end-2))+...
        ((a.*(Vy_new(3:end-2,3:end-2,3:end-2)-Vy_new(2:end-3,3:end-2,3:end-2))+...
        b2.*(Vy_new(4:end-1,3:end-2,3:end-2)-Vy_new(1:end-4,3:end-2,3:end-2)))./Dhy(3:end-2))+...
        ((a.*(Vz_new(3:end-2,3:end-2,3:end-2)-Vz_new(3:end-2,3:end-2,2:end-3))+...
        b2.*(Vz_new(3:end-2,3:end-2,4:end-1)-Vz_new(3:end-2,3:end-2,1:end-4)))./Dhz(3:end-2))))+...
        (2.*Dt.*mu(3:end-2,3:end-2,3:end-2)./tsl(indtii_tsl1,indtii_tsl2,indtii_tsl3)).*...
        TS(indtii_TS1,indtii_TS2,indtii_TS3).*...
        (((a.*(Vx_new(3:end-2,3:end-2,3:end-2)-Vx_new(3:end-2,2:end-3,3:end-2))+...
        b2.*(Vx_new(3:end-2,4:end-1,3:end-2)-Vx_new(3:end-2,1:end-4,3:end-2)))./Dhx(3:end-2))+...
        ((a.*(Vz_new(3:end-2,3:end-2,3:end-2)-Vz_new(3:end-2,3:end-2,2:end-3))+...
        b2.*(Vz_new(3:end-2,3:end-2,4:end-1)-Vz_new(3:end-2,3:end-2,1:end-4)))./Dhz(3:end-2))));

    % Calculate rzz relaation funciton
    rzz(3:end-2,3:end-2,3:end-2)=((1+(Dt./2.*tsl(indtii_tsl1,indtii_tsl2,indtii_tsl3))).^-1).*...
        (((1-(Dt./2.*tsl(indtii_tsl1,indtii_tsl2,indtii_tsl3))).*rzz(3:end-2,3:end-2,3:end-2))-...
        ((lamda(3:end-2,3:end-2,3:end-2).*Dt./tsl(indtii_tsl1,indtii_tsl2,indtii_tsl3)).*...
        TP(indtii_TP1,indtii_TP2,indtii_TP3).*...
        (((a.*(Vx_new(3:end-2,3:end-2,3:end-2)-Vx_new(3:end-2,2:end-3,3:end-2))+...
        b2.*(Vx_new(3:end-2,4:end-1,3:end-2)-Vx_new(3:end-2,1:end-4,3:end-2)))./Dhx(3:end-2))+...
        ((a.*(Vy_new(3:end-2,3:end-2,3:end-2)-Vy_new(2:end-3,3:end-2,3:end-2))+...
        b2.*(Vy_new(4:end-1,3:end-2,3:end-2)-Vy_new(1:end-4,3:end-2,3:end-2)))./Dhy(3:end-2))+...
        ((a.*(Vz_new(3:end-2,3:end-2,3:end-2)-Vz_new(3:end-2,3:end-2,2:end-3))+...
        b2.*(Vz_new(3:end-2,3:end-2,4:end-1)-Vz_new(3:end-2,3:end-2,1:end-4)))./Dhz(3:end-2))))+...
        (2.*Dt.*mu(3:end-2,3:end-2,3:end-2)./tsl(indtii_tsl1,indtii_tsl2,indtii_tsl3)).*...
        TS(indtii_TS1,indtii_TS2,indtii_TS3).*...
        (((a.*(Vx_new(3:end-2,3:end-2,3:end-2)-Vx_new(3:end-2,2:end-3,3:end-2))+...
        b2.*(Vx_new(3:end-2,4:end-1,3:end-2)-Vx_new(3:end-2,1:end-4,3:end-2)))./Dhx(3:end-2))+...
        ((a.*(Vy_new(3:end-2,3:end-2,3:end-2)-Vy_new(2:end-3,3:end-2,3:end-2))+...
        b2.*(Vy_new(4:end-1,3:end-2,3:end-2)-Vy_new(1:end-4,3:end-2,3:end-2)))./Dhy(3:end-2))));

    % Calculate rxy relaation funciton
    rxy(3:end-1,3:end-1,2:end-2)=((1+(Dt./2.*tsl(indtxy_tsl1,indtxy_tsl2,indtxy_tsl3))).^-1).*...
        (((1-(Dt./2.*tsl(indtxy_tsl1,indtxy_tsl2,indtxy_tsl3))).*rxy(3:end-1,3:end-1,2:end-2))-...
        ((mu(3:end-1,3:end-1,2:end-2).*Dt./tsl(indtxy_tsl1,indtxy_tsl2,indtxy_tsl3)).*...
        TS(indtxy_TS1,indtxy_TS2,indtxy_TS3).*...
        (((a.*(Vx_new(3:end-1,2:end-2,2:end-2)-Vx_new(2:end-2,2:end-2,2:end-2))+...
        b2.*(Vx_new(4:end,2:end-2,2:end-2)-Vx_new(1:end-3,2:end-2,2:end-2)))./Dhy(2:end-2))+...
        ((a.*(Vy_new(2:end-2,3:end-1,2:end-2)-Vy_new(2:end-2,2:end-2,2:end-2))+...
        b2.*(Vy_new(2:end-2,4:end,2:end-2)-Vy_new(2:end-2,1:end-3,2:end-2)))./Dhx(2:end-2)))));

    % Calculate rxz relaation funciton
    rxz(2:end-2,3:end-1,3:end-1)=((1+(Dt./2.*tsl(indtxz_tsl1,indtxz_tsl2,indtxz_tsl3))).^-1).*...
        (((1-(Dt./2.*tsl(indtxz_tsl1,indtxz_tsl2,indtxz_tsl3))).*rxz(2:end-2,3:end-1,3:end-1))-...
        ((mu(2:end-2,3:end-1,3:end-1).*Dt./tsl(indtxz_tsl1,indtxz_tsl2,indtxz_tsl3)).*...
        TS(indtxz_TS1,indtxz_TS2,indtxz_TS3).*...
        (((a.*(Vx_new(2:end-2,2:end-2,3:end-1)-Vx_new(2:end-2,2:end-2,2:end-2))+...
        b2.*(Vx_new(2:end-2,2:end-2,4:end)-Vx_new(2:end-2,2:end-2,1:end-3)))./Dhz(2:end-2))+...
        ((a.*(Vz_new(2:end-2,3:end-1,2:end-2)-Vz_new(2:end-2,2:end-2,2:end-2))+...
        b2.*(Vz_new(2:end-2,4:end,2:end-2)-Vz_new(2:end-2,1:end-3,2:end-2)))./Dhx(2:end-2)))));

    % Calculate ryz relaation funciton
    ryz(3:end-1,2:end-2,3:end-1)=((1+(Dt./2.*tsl(indtyz_tsl1,indtyz_tsl2,indtyz_tsl3))).^-1).*...
        (((1-(Dt./2.*tsl(indtyz_tsl1,indtyz_tsl2,indtyz_tsl3))).*ryz(3:end-1,2:end-2,3:end-1))-...
        ((mu(3:end-1,2:end-2,3:end-1).*Dt./tsl(indtyz_tsl1,indtyz_tsl2,indtyz_tsl3)).*...
        TS(indtyz_TS1,indtyz_TS2,indtyz_TS3).*...
        (((a.*(Vz_new(3:end-1,2:end-2,2:end-2)-Vz_new(2:end-2,2:end-2,2:end-2))+...
        b2.*(Vz_new(4:end,2:end-2,2:end-2)-Vz_new(1:end-3,2:end-2,2:end-2)))./Dhy(2:end-2))+...
        ((a.*(Vy_new(2:end-2,2:end-2,3:end-1)-Vy_new(2:end-2,2:end-2,2:end-2))+...
        b2.*(Vy_new(2:end-2,2:end-2,4:end)-Vy_new(2:end-2,2:end-2,3:end-1)))./Dhz(2:end-2)))));

    % Stress field calculation
    % Calculate txx stress component
    txx(3:end-2,3:end-2,3:end-2)=txx(3:end-2,3:end-2,3:end-2)+(Dt.*lamda(3:end-2,3:end-2,3:end-2)).*...
        ((1+TP(indtii_TP1,indtii_TP2,indtii_TP3)).*...
        (((a.*(Vx_new(3:end-2,3:end-2,3:end-2)-Vx_new(3:end-2,2:end-3,3:end-2))+...
        b2.*(Vx_new(3:end-2,4:end-1,3:end-2)-Vx_new(3:end-2,1:end-4,3:end-2)))./Dhx(3:end-2))+...
        ((a.*(Vy_new(3:end-2,3:end-2,3:end-2)-Vy_new(2:end-3,3:end-2,3:end-2))+...
        b2.*(Vy_new(4:end-1,3:end-2,3:end-2)-Vy_new(1:end-4,3:end-2,3:end-2)))./Dhy(3:end-2))+...
        ((a.*(Vz_new(3:end-2,3:end-2,3:end-2)-Vz_new(3:end-2,3:end-2,2:end-3))+...
        b2.*(Vz_new(3:end-2,3:end-2,4:end-1)-Vz_new(3:end-2,3:end-2,1:end-4)))./Dhz(3:end-2))))-...
        (Dt.*2.*mu(3:end-2,3:end-2,3:end-2)).*((1+TS(indtii_TS1,indtii_TS2,indtii_TS3)).*...
        (((a.*(Vy_new(3:end-2,3:end-2,3:end-2)-Vy_new(2:end-3,3:end-2,3:end-2))+...
        b2.*(Vy_new(4:end-1,3:end-2,3:end-2)-Vy_new(1:end-4,3:end-2,3:end-2)))./Dhy(3:end-2))+...
        ((a.*(Vz_new(3:end-2,3:end-2,3:end-2)-Vz_new(3:end-2,3:end-2,2:end-3))+...
        b2.*(Vz_new(3:end-2,3:end-2,4:end-1)-Vz_new(3:end-2,3:end-2,1:end-4)))./Dhz(3:end-2))))+...
        (Dt/2).*(rxx(3:end-2,3:end-2,3:end-2));

    % Calculate tyy stress component
    tyy(3:end-2,3:end-2,3:end-2)=tyy(3:end-2,3:end-2,3:end-2)+(Dt.*lamda(3:end-2,3:end-2,3:end-2)).*...
        ((1+TP(indtii_TP1,indtii_TP2,indtii_TP3)).*...
        (((a.*(Vx_new(3:end-2,3:end-2,3:end-2)-Vx_new(3:end-2,2:end-3,3:end-2))+...
        b2.*(Vx_new(3:end-2,4:end-1,3:end-2)-Vx_new(3:end-2,1:end-4,3:end-2)))./Dhx(3:end-2))+...
        ((a.*(Vy_new(3:end-2,3:end-2,3:end-2)-Vy_new(2:end-3,3:end-2,3:end-2))+...
        b2.*(Vy_new(4:end-1,3:end-2,3:end-2)-Vy_new(1:end-4,3:end-2,3:end-2)))./Dhy(3:end-2))+...
        ((a.*(Vz_new(3:end-2,3:end-2,3:end-2)-Vz_new(3:end-2,3:end-2,2:end-3))+...
        b2.*(Vz_new(3:end-2,3:end-2,4:end-1)-Vz_new(3:end-2,3:end-2,1:end-4)))./Dhz(3:end-2))))-...
        (Dt.*2.*mu(3:end-2,3:end-2,3:end-2)).*((1+TS(indtii_TS1,indtii_TS2,indtii_TS3)).*...
        (((a.*(Vx_new(3:end-2,3:end-2,3:end-2)-Vx_new(3:end-2,2:end-3,3:end-2))+...
        b2.*(Vx_new(3:end-2,4:end-1,3:end-2)-Vx_new(3:end-2,1:end-4,3:end-2)))./Dhx(3:end-2))+...
        ((a.*(Vz_new(3:end-2,3:end-2,3:end-2)-Vz_new(3:end-2,3:end-2,2:end-3))+...
        b2.*(Vz_new(3:end-2,3:end-2,4:end-1)-Vz_new(3:end-2,3:end-2,1:end-4)))./Dhz(3:end-2))))+...
        (Dt/2).*(ryy(3:end-2,3:end-2,3:end-2));

    % Calculate tzz stress component
    tzz(3:end-2,3:end-2,3:end-2)=tzz(3:end-2,3:end-2,3:end-2)+(Dt.*lamda(3:end-2,3:end-2,3:end-2)).*...
        ((1+TP(indtii_TP1,indtii_TP2,indtii_TP3)).*...
        (((a.*(Vx_new(3:end-2,3:end-2,3:end-2)-Vx_new(3:end-2,2:end-3,3:end-2))+...
        b2.*(Vx_new(3:end-2,4:end-1,3:end-2)-Vx_new(3:end-2,1:end-4,3:end-2)))./Dhx(3:end-2))+...
        ((a.*(Vy_new(3:end-2,3:end-2,3:end-2)-Vy_new(2:end-3,3:end-2,3:end-2))+...
        b2.*(Vy_new(4:end-1,3:end-2,3:end-2)-Vy_new(1:end-4,3:end-2,3:end-2)))./Dhy(3:end-2))+...
        ((a.*(Vz_new(3:end-2,3:end-2,3:end-2)-Vz_new(3:end-2,3:end-2,2:end-3))+...
        b2.*(Vz_new(3:end-2,3:end-2,4:end-1)-Vz_new(3:end-2,3:end-2,1:end-4)))./Dhz(3:end-2))))-...
        (Dt.*2.*mu(3:end-2,3:end-2,3:end-2)).*((1+TS(indtii_TS1,indtii_TS2,indtii_TS3)).*...
        (((a.*(Vx_new(3:end-2,3:end-2,3:end-2)-Vx_new(3:end-2,2:end-3,3:end-2))+...
        b2.*(Vx_new(3:end-2,4:end-1,3:end-2)-Vx_new(3:end-2,1:end-4,3:end-2)))./Dhx(3:end-2))+...
        ((a.*(Vy_new(3:end-2,3:end-2,3:end-2)-Vy_new(2:end-3,3:end-2,3:end-2))+...
        b2.*(Vy_new(4:end-1,3:end-2,3:end-2)-Vy_new(1:end-4,3:end-2,3:end-2)))./Dhy(3:end-2))))+...
        (Dt/2).*(rzz(3:end-2,3:end-2,3:end-2));

    % Calculate txy stress component
    txy(3:end-1,3:end-1,2:end-2)=txy(3:end-1,3:end-1,2:end-2)+...
        (mu(3:end-1,3:end-1,2:end-2).*Dt).*((1+TS(indtxy_TS1,indtxy_TS2,indtxy_TS3)).*...
        (((a.*(Vx_new(3:end-1,2:end-2,2:end-2)-Vx_new(2:end-2,2:end-2,2:end-2))+...
        b2.*(Vx_new(4:end,2:end-2,2:end-2)-Vx_new(1:end-3,2:end-2,2:end-2)))./Dhy(2:end-2))+...
        ((a.*(Vy_new(2:end-2,3:end-1,2:end-2)-Vy_new(2:end-2,2:end-2,2:end-2))+...
        b2.*(Vy_new(2:end-2,4:end,2:end-2)-Vy_new(2:end-2,1:end-3,2:end-2)))./Dhx(2:end-2))))+...
        (Dt/2).*(rxy(3:end-1,3:end-1,2:end-2));

    % Calculate txz stress component
    txz(2:end-2,3:end-1,3:end-1)=txz(2:end-2,3:end-1,3:end-1)+...
        (mu(2:end-2,3:end-1,3:end-1).*Dt).*((1+TS(indtxz_TS1,indtxz_TS2,indtxz_TS3)).*...
        (((a.*(Vx_new(2:end-2,2:end-2,3:end-1)-Vx_new(2:end-2,2:end-2,2:end-2))+...
        b2.*(Vx_new(2:end-2,2:end-2,4:end)-Vx_new(2:end-2,2:end-2,1:end-3)))./Dhz(2:end-2))+...
        ((a.*(Vz_new(2:end-2,3:end-1,2:end-2)-Vz_new(2:end-2,2:end-2,2:end-2))+...
        b2.*(Vz_new(2:end-2,4:end,2:end-2)-Vz_new(2:end-2,1:end-3,2:end-2)))./Dhx(2:end-2))))+...
        (Dt/2).*(rxz(2:end-2,3:end-1,3:end-1));

    % Calculate tyz stress component
    tyz(3:end-1,2:end-2,3:end-1)=tyz(3:end-1,2:end-2,3:end-1)+...
        (mu(3:end-1,2:end-2,3:end-1).*Dt).*((1+TS(indtyz_TS1,indtyz_TS2,indtyz_TS3)).*...
        (((a.*(Vz_new(3:end-1,2:end-2,2:end-2)-Vz_new(2:end-2,2:end-2,2:end-2))+...
        b2.*(Vz_new(4:end,2:end-2,2:end-2)-Vz_new(1:end-3,2:end-2,2:end-2)))./Dhy(2:end-2))+...
        ((a.*(Vy_new(2:end-2,2:end-2,3:end-1)-Vy_new(2:end-2,2:end-2,2:end-2))+...
        b2.*(Vy_new(2:end-2,2:end-2,4:end)-Vy_new(2:end-2,2:end-2,3:end-1)))./Dhz(2:end-2))))+...
        (Dt/2).*(ryz(3:end-1,2:end-2,3:end-1));


    % Apply free surface boundary conditions
    if topoload ==1
        % If arbitary surface is used
        txx(g_mask_tzzi)=0;
        tyy(g_mask_tzzi)=0;
        tzz(g_mask_tzzi)=0;
        txy(g_mask_txyi)=0;
        txz(g_mask_txzi)=0;
        tyz(g_mask_tyzi)=0;
    else
        % If flat surface is used
        tzz(:,:,1)=0;tyy(:,:,1)=0;tzz(:,:,1)=0;txy(:,:,1)=0;
        Vx_new(:,:,1)=0;Vy_new(:,:,1)=0;
        txz(:,:,1)=0;tyz(:,:,1)=0;
    end


    % Update save counter
    if kk==skipt
            kk=1;
    else
        kk=kk+1;
    end
    
    % Stop script execution if simulation is unstable
    if max(abs(Vx_new(:)))==inf
        msgbox('Reached inf amplitude value. Aborting simulation','input');
        return
    end
end
durb=toc;
disp(['Preprocessing duration: ' num2str(dura) ' sec'])
disp(['Simulation duration: ' num2str(durb) ' sec'])
cd(orpath)
end
