function w3=FD_sim_fault_2(Vp,Vs,rho,X,Y,Z,source,pulse,Q,f0,Dt,maxt,posx,posy,posz,fault,xrup,zrup,gpumode,N,b,topoload,topo,X_topo,Y_topo,savemode,savetype,skipt,outpath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%INPUT DATA%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vp: 3D P-velocity grid
% Vp: 3D S-velocity grid
% rho: 3D Pdensity grid
% X,Y,Z: x,y,z node locations (e.g. 0 100 200,... meters)
% source: vector of 5 elements: [strike dip rake mantissa exponent]
% pulse: N x 2 matrix. pulse(:,1)=time, pulse(:,2)= amplitude of source
% Q: Single value or 3D matrix with Q0
% f0: Single value or 3D matrix with central frequency of attenuation
% Dt: Simulation sampliong period
% maxt: maximum simulation steps
% posx, posy, posz: x,y,z node index of source (or fault upper left edge)
% fault: Fault length, width and rupture velocity ([i, j, k])
% xrup: X location of rupture initiation within fault (xrup.*length)
% zrup: Y (Z), Location of rupture initiation within fault (zrup.*width)
% gpumode: 0-1 disable or enable gpu processing
% N: Number of nodes used as transition zone for ABC
% b
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
clear max_topo max_model X_topo Y_topo w3 temp1 X_out Y_out Z_out
%% Prepare source 
disp('Adjusting fault plane to calculation grid');
X=single(X);Y=single(Y);Z=single(Z);
% Adjust fault surface to calculation grid
[faultpoints, ttime, nfact]=fin_fault(posx,posy,posz,X,Y,Z,fault(1),fault(2),source(1),-source(2),xrup,zrup,fault(3),Dt,min(abs(Dhx(:))));
  
% Obtain fault nodes indices for all velocity and stress nodes
[faultpoints_x0_y0_z0,faultpoints_x1_y0_z0, faultpoints_x1_y1_z0, ...
faultpoints_x1_y11_z0,faultpoints_x0_y1_z0,faultpoints_x0_y11_z0,...
faultpoints_x0_y0_z1,faultpoints_x0_y0_z11,faultpoints_x1_y0_z11,...
faultpoints_x11_y0_z0,faultpoints_x11_y1_z0, faultpoints_x0_y1_z11,...
faultpoints_x0_y1_z1,faultpoints_x1_y0_z1,faultpoints_x11_y0_z1,...
faultpoints_x0_y11_z1]=fault_shift(faultpoints);

% Interpolate and adjust source time function to simulation sampling period
a=interp1(pulse(:,1),pulse(:,2),0:Dt:(maxt*Dt)-Dt);  
a(isnan(a))=0;
a(abs(a)<0.0000001)=0;
a=a./max(abs(a));
% Expand STF duration to cover all sub-source nodes
pulse3(1:max(ttime(:)))=0;
pulse3(end+1:end+length(a))=a;
pulse3(end+1:end+10000)=0;
pulse3=pulse3./(max(abs(pulse3)));
clear a
a=pulse3;
a=(a./max(a));

% Convert strike/dip/rake to 6 moment tensor components
[Gtxx,Gtxy,Gtxz,Gtyy,Gtyz,Gtzz]=sdr2mt(source(1)-90,source(2),source(3),source(4));

disp('Finished preparing finite source');
disp(['Total number of sub-sources: ' num2str(round(nfact))])
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
clear Vp Vs X Y

if gpumode==1
    % Preallocate velocity component matrices
    Vx_new=gpuArray(single(zeros(ny,nx,nz)));   
    Vy_new=gpuArray(single(zeros(ny,nx,nz)));
    Vz_new=gpuArray(single(zeros(ny,nx,nz)));
    
    % Preallocate matrices for the application of boundary conditions
    Vz_old_y1=gpuArray(Vz_new(1:N+1,:,:));Vx_old_y1=gpuArray(Vx_new(1:N+1,:,:));Vy_old_y1=gpuArray(Vy_new(1:N+1,:,:)); %y=1 boundary
    Vz_old_x1=gpuArray(Vz_new(:,1:N+1,:));Vy_old_x1=gpuArray(Vy_new(:,1:N+1,:));Vx_old_x1=gpuArray(Vx_new(:,1:N+1,:)); %x=1 boundary
    Vz_old_yend=gpuArray(Vz_new(end-N:end,:,:));Vx_old_yend=gpuArray(Vx_new(end-N:end,:,:));Vy_old_yend=gpuArray(Vy_new(end-N:end,:,:)); %y=end boundary
    Vz_old_xend=gpuArray(Vz_new(:,end-N:end,:));Vx_old_xend=gpuArray(Vx_new(:,end-N:end,:));Vy_old_xend=gpuArray(Vy_new(:,end-N:end,:)); %x=end boundary
    Vz_old_zend=gpuArray(Vz_new(:,:,end-N:end));Vx_old_zend=gpuArray(Vx_new(:,:,end-N:end));Vy_old_zend=gpuArray(Vy_new(:,:,end-N:end)); %z=end boundary

    % Preallocate stress matrices
    tyy=gpuArray(single(zeros(ny,nx,nz)));txy=gpuArray(single(zeros(ny,nx,nz)));tyz=gpuArray(single(zeros(ny,nx,nz)));
    txx=gpuArray(single(zeros(ny,nx,nz)));tzz=gpuArray(single(zeros(ny,nx,nz)));txz=gpuArray(single(zeros(ny,nx,nz)));

    rho=gpuArray(rho);
    lamda=gpuArray(lamda);
    mu=gpuArray(mu);
    
    % Preallocate constant attenuation matrix (or singular value if no 3D matrix was given as input) 
    Att=gpuArray(single(exp((-pi.*f0.*Dt)./Q)));
        
    faultpoints_x0_y0_z0=gpuArray(faultpoints_x0_y0_z0);
	faultpoints_x1_y0_z0=gpuArray(faultpoints_x1_y0_z0);
	faultpoints_x1_y1_z0=gpuArray(faultpoints_x1_y1_z0);
    faultpoints_x1_y11_z0=gpuArray(faultpoints_x0_y11_z0);
	faultpoints_x0_y1_z0=gpuArray(faultpoints_x0_y1_z0);
	faultpoints_x0_y11_z0=gpuArray(faultpoints_x0_y11_z0);
    faultpoints_x0_y0_z1=gpuArray(faultpoints_x0_y0_z1);
	faultpoints_x0_y0_z11=gpuArray(faultpoints_x0_y0_z11);
	faultpoints_x1_y0_z11=gpuArray(faultpoints_x1_y0_z11);
    faultpoints_x11_y0_z0=gpuArray(faultpoints_x11_y0_z0);
	faultpoints_x11_y1_z0=gpuArray(faultpoints_x11_y1_z0);
	faultpoints_x0_y1_z11=gpuArray(faultpoints_x0_y1_z11);
    faultpoints_x0_y1_z1=gpuArray(faultpoints_x0_y1_z1);
	faultpoints_x1_y0_z1=gpuArray(faultpoints_x1_y0_z1);
	faultpoints_x11_y0_z1=gpuArray(faultpoints_x11_y0_z1);
    faultpoints_x0_y11_z1=gpuArray(faultpoints_x0_y11_z1);
    

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

    % Preallocate constant attenuation matrix (or singular value if no 3D matrix was given as input)
    Att=single(exp((-pi.*f0.*Dt)./Q));
  
    F_txx=((10^source(5)).*Gtxx.*a);
    F_tyy=((10^source(5)).*Gtyy.*a);
    F_tzz=((10^source(5)).*Gtzz.*a);
    F_txy=((10^source(5)).*Gtxy.*a);
    F_txz=((10^source(5)).*Gtxz.*a);
    F_tyz=((10^source(5)).*Gtyz.*a);
end

% Check if 3D attenuation is given and set appropriate index of nodes for
% its application on all velocity and stress components. This step is
% neccessary in order to avoid using a 3D matrix (or a separate script) if
% only a singular attenuation value is given.
if length(Att)==1
    indV_Q1=1;indV_Q2=1;indV_Q3=1;
    indtii_Q1=1;indtii_Q2=1;indtii_Q3=1;
    indtxy_Q1=1;indtxy_Q2=1;indtxy_Q3=1;
    indtxz_Q1=1;indtxz_Q2=1;indtxz_Q3=1;
    indtyz_Q1=1;indtyz_Q2=1;indtyz_Q3=1;
else
    indV_Q1=1:ny-1;indV_Q2=1:nx-1;indV_Q3=1:nz-1;
    indtii_Q1=2:ny-1;indtii_Q2=2:nx-1;indtii_Q3=2:nz-1;
    indtxy_Q1=2:ny;indtxy_Q2=2:nx;indtxy_Q3=1:nz-1;
    indtxz_Q1=1:ny-1;indtxz_Q2=2:nx;indtxz_Q3=2:nz;
    indtyz_Q1=2:ny;indtyz_Q2=1:nx-1;indtyz_Q3=2:nz;
end

dura=toc;
%% Main calculation section
tic
kk=1;tt=1;

for t=1:maxt
    
    disp([num2str(100*t/maxt) '%'])
    
    % Calculate X velocity component for whole grid    
    Vx_new(1:end-1,1:end-1,1:end-1)=Att(indV_Q1,indV_Q2,indV_Q3).*(Vx_new(1:end-1,1:end-1,1:end-1)+...
        ((Dt.*(1./(rho(1:end-1,1:end-1,1:end-1)))).*...
        (((txx(1:end-1,2:end,1:end-1)-txx(1:end-1,1:end-1,1:end-1))./Dhx(1:end-1))+...
        ((txy(2:end,2:end,1:end-1)-txy(1:end-1,2:end,1:end-1))./Dhy(1:end-1))+...
        ((txz(1:end-1,2:end,2:end)-txz(1:end-1,2:end,1:end-1))./Dhz(1:end-1)))));

    % Insert source at Vx
    Vx_new(faultpoints_x0_y0_z0)=Vx_new(faultpoints_x0_y0_z0)+(Dt./(rho(faultpoints_x0_y0_z0))).*(F_txx(t+max(ttime(:))-ttime(faultpoints_x0_y0_z0))'./((Dhx(posx).^2).*Dhy(posy).*Dhz(posz)));%
    Vx_new(faultpoints_x1_y0_z0)=Vx_new(faultpoints_x1_y0_z0)+(Dt./(rho(faultpoints_x1_y0_z0))).*(F_txx(t+max(ttime(:))-ttime(faultpoints_x1_y0_z0))'./((Dhx(posx).^2).*Dhy(posy).*Dhz(posz)));% 
    Vx_new(faultpoints_x1_y11_z0)=Vx_new(faultpoints_x1_y11_z0)+(Dt./(rho(faultpoints_x1_y11_z0))).*(F_txy(t+max(ttime(:))-ttime(faultpoints_x1_y11_z0))'./((Dhx(posx).^2).*4.*Dhy(posy).*Dhz(posz)));% 
    Vx_new(faultpoints_x1_y1_z0)=Vx_new(faultpoints_x1_y1_z0)+(Dt./(rho(faultpoints_x1_y1_z0))).*(F_txy(t+max(ttime(:))-ttime(faultpoints_x1_y1_z0))'./((Dhx(posx).^2).*4.*Dhy(posy).*Dhz(posz)));% 
    Vx_new(faultpoints_x0_y1_z0)=Vx_new(faultpoints_x0_y1_z0)+(Dt./(rho(faultpoints_x0_y1_z0))).*(F_txy(t+max(ttime(:))-ttime(faultpoints_x0_y1_z0))'./((Dhx(posx).^2).*4.*Dhy(posy).*Dhz(posz)));% 
    Vx_new(faultpoints_x0_y11_z0)=Vx_new(faultpoints_x0_y11_z0)+(Dt./(rho(faultpoints_x0_y11_z0))).*(F_txy(t+max(ttime(:))-ttime(faultpoints_x0_y11_z0))'./((Dhx(posx).^2).*4.*Dhy(posy).*Dhz(posz)));% 
    Vx_new(faultpoints_x0_y0_z11)=Vx_new(faultpoints_x0_y0_z11)+(Dt./(rho(faultpoints_x0_y0_z11))).*(F_txz(t+max(ttime(:))-ttime(faultpoints_x0_y0_z11))'./((Dhx(posx).^2).*4.*Dhy(posy).*Dhz(posz)));% 
    Vx_new(faultpoints_x1_y0_z11)=Vx_new(faultpoints_x1_y0_z11)+(Dt./(rho(faultpoints_x1_y0_z11))).*(F_txz(t+max(ttime(:))-ttime(faultpoints_x1_y0_z11))'./((Dhx(posx).^2).*4.*Dhy(posy).*Dhz(posz)));% 
    Vx_new(faultpoints_x1_y0_z1)=Vx_new(faultpoints_x1_y0_z1)+(Dt./(rho(faultpoints_x1_y0_z1))).*(F_txz(t+max(ttime(:))-ttime(faultpoints_x1_y0_z1))'./((Dhx(posx).^2).*4.*Dhy(posy).*Dhz(posz)));% 
    Vx_new(faultpoints_x0_y0_z1)=Vx_new(faultpoints_x0_y0_z1)+(Dt./(rho(faultpoints_x0_y0_z1))).*(F_txz(t+max(ttime(:))-ttime(faultpoints_x0_y0_z1))'./((Dhx(posx).^2).*4.*Dhy(posy).*Dhz(posz)));% 

    % Calculate Y velocity component for whole grid
    Vy_new(1:end-1,1:end-1,1:end-1)=Att(indV_Q1,indV_Q2,indV_Q3).*(Vy_new(1:end-1,1:end-1,1:end-1)+...
        ((Dt.*(1./(rho(1:end-1,1:end-1,1:end-1)))).*...
        (((tyy(2:end,1:end-1,1:end-1)-tyy(1:end-1,1:end-1,1:end-1))./Dhy(1:end-1))+...
        ((txy(2:end,2:end,1:end-1)-txy(2:end,1:end-1,1:end-1))./Dhx(1:end-1))+...
        ((tyz(2:end,1:end-1,2:end)-tyz(2:end,1:end-1,1:end-1))./Dhz(1:end-1)))));

    % Insert source at Vy
    Vy_new(faultpoints_x0_y0_z0)=Vy_new(faultpoints_x0_y0_z0)+(Dt./(rho(faultpoints_x0_y0_z0))).*(F_tyy(t+max(ttime(:))-ttime(faultpoints_x0_y0_z0))'./(Dhx(posx).*(Dhy(posy).^2).*Dhz(posz)));% 
    Vy_new(faultpoints_x0_y1_z0)=Vy_new(faultpoints_x0_y1_z0)+(Dt./(rho(faultpoints_x0_y1_z0))).*(F_tyy(t+max(ttime(:))-ttime(faultpoints_x0_y1_z0))'./(Dhx(posx).*(Dhy(posy).^2).*Dhz(posz)));% 
    Vy_new(faultpoints_x11_y1_z0)=Vy_new(faultpoints_x11_y1_z0)+(Dt./(rho(faultpoints_x11_y1_z0))).*(F_txy(t+max(ttime(:))-ttime(faultpoints_x11_y1_z0))'./(4.*Dhx(posx).*(Dhy(posy).^2).*Dhz(posz)));%
    Vy_new(faultpoints_x1_y1_z0)=Vy_new(faultpoints_x1_y1_z0)+(Dt./(rho(faultpoints_x1_y1_z0))).*(F_txy(t+max(ttime(:))-ttime(faultpoints_x1_y1_z0))'./(4.*Dhx(posx).*(Dhy(posy).^2).*Dhz(posz)));%
    Vy_new(faultpoints_x1_y0_z0)=Vy_new(faultpoints_x1_y0_z0)+(Dt./(rho(faultpoints_x1_y0_z0))).*(F_txy(t+max(ttime(:))-ttime(faultpoints_x1_y0_z0))'./(4.*Dhx(posx).*(Dhy(posy).^2).*Dhz(posz)));%
    Vy_new(faultpoints_x11_y0_z0)=Vy_new(faultpoints_x11_y0_z0)+(Dt./(rho(faultpoints_x11_y0_z0))).*(F_txy(t+max(ttime(:))-ttime(faultpoints_x11_y0_z0))'./(4.*Dhx(posx).*(Dhy(posy).^2).*Dhz(posz)));%
    Vy_new(faultpoints_x0_y0_z11)=Vy_new(faultpoints_x0_y0_z11)+(Dt./(rho(faultpoints_x0_y0_z11))).*(F_tyz(t+max(ttime(:))-ttime(faultpoints_x0_y0_z11))'./(4.*Dhx(posx).*(Dhy(posy).^2).*Dhz(posz)));%
    Vy_new(faultpoints_x0_y1_z11)=Vy_new(faultpoints_x0_y1_z11)+(Dt./(rho(faultpoints_x0_y1_z11))).*(F_tyz(t+max(ttime(:))-ttime(faultpoints_x0_y1_z11))'./(4.*Dhx(posx).*(Dhy(posy).^2).*Dhz(posz)));%
    Vy_new(faultpoints_x0_y1_z1)=Vy_new(faultpoints_x0_y1_z1)+(Dt./(rho(faultpoints_x0_y1_z1))).*(F_tyz(t+max(ttime(:))-ttime(faultpoints_x0_y1_z1))'./(4.*Dhx(posx).*(Dhy(posy).^2).*Dhz(posz)));%
    Vy_new(faultpoints_x0_y1_z1)=Vy_new(faultpoints_x0_y1_z1)+(Dt./(rho(faultpoints_x0_y1_z1))).*(F_tyz(t+max(ttime(:))-ttime(faultpoints_x0_y1_z1))'./(4.*Dhx(posx).*(Dhy(posy).^2).*Dhz(posz)));%

    % Calculate Z velocity component for whole grid
    Vz_new(1:end-1,1:end-1,1:end-1)=Att(indV_Q1,indV_Q2,indV_Q3).*(Vz_new(1:end-1,1:end-1,1:end-1)+...
        ((Dt.*(1./(rho(1:end-1,1:end-1,1:end-1)))).*...
        (((tzz(1:end-1,1:end-1,2:end)-tzz(1:end-1,1:end-1,1:end-1))./Dhz(1:end-1))+...
        ((txz(1:end-1,2:end,2:end)-txz(1:end-1,1:end-1,2:end))./Dhx(1:end-1))+...
        ((tyz(2:end,1:end-1,2:end)-tyz(1:end-1,1:end-1,2:end))./Dhy(1:end-1)))));
    
    % Insert source at Vz
    Vz_new(faultpoints_x0_y0_z0)=Vz_new(faultpoints_x0_y0_z0)+(Dt./(rho(faultpoints_x0_y0_z0))).*(F_tzz(t+max(ttime(:))-ttime(faultpoints_x0_y0_z0))'./(Dhx(posx).*Dhy(posy).*(Dhz(posz).^2)));% 
    Vz_new(faultpoints_x0_y0_z1)=Vz_new(faultpoints_x0_y0_z1)+(Dt./(rho(faultpoints_x0_y0_z1))).*(F_tzz(t+max(ttime(:))-ttime(faultpoints_x0_y0_z1))'./(Dhx(posx).*Dhy(posy).*(Dhz(posz).^2)));% 
    Vz_new(faultpoints_x11_y0_z1)=Vz_new(faultpoints_x11_y0_z1)+(Dt./(rho(faultpoints_x11_y0_z1))).*(F_txz(t+max(ttime(:))-ttime(faultpoints_x11_y0_z1))'./(4.*Dhx(posx).*Dhy(posy).*(Dhz(posz).^2)));% 
    Vz_new(faultpoints_x1_y0_z1)=Vz_new(faultpoints_x1_y0_z1)+(Dt./(rho(faultpoints_x1_y0_z1))).*(F_txz(t+max(ttime(:))-ttime(faultpoints_x1_y0_z1))'./(4.*Dhx(posx).*Dhy(posy).*(Dhz(posz).^2)));%  
    Vz_new(faultpoints_x1_y0_z0)=Vz_new(faultpoints_x1_y0_z0)+(Dt./(rho(faultpoints_x1_y0_z0))).*(F_txz(t+max(ttime(:))-ttime(faultpoints_x1_y0_z0))'./(4.*Dhx(posx).*Dhy(posy).*(Dhz(posz).^2)));% 
    Vz_new(faultpoints_x11_y0_z0)=Vz_new(faultpoints_x11_y0_z0)+(Dt./(rho(faultpoints_x11_y0_z0))).*(F_txz(t+max(ttime(:))-ttime(faultpoints_x11_y0_z0))'./(4.*Dhx(posx).*Dhy(posy).*(Dhz(posz).^2)));% 
    Vz_new(faultpoints_x0_y1_z0)=Vz_new(faultpoints_x0_y1_z0)+(Dt./(rho(faultpoints_x0_y1_z0))).*(F_tyz(t+max(ttime(:))-ttime(faultpoints_x0_y1_z0))'./(4.*Dhx(posx).*Dhy(posy).*(Dhz(posz).^2)));% 
    Vz_new(faultpoints_x0_y11_z0)=Vz_new(faultpoints_x0_y11_z0)+(Dt./(rho(faultpoints_x0_y11_z0))).*(F_tyz(t+max(ttime(:))-ttime(faultpoints_x0_y11_z0))'./(4.*Dhx(posx).*Dhy(posy).*(Dhz(posz).^2)));%
    Vz_new(faultpoints_x0_y11_z1)=Vz_new(faultpoints_x0_y11_z1)+(Dt./(rho(faultpoints_x0_y11_z1))).*(F_tyz(t+max(ttime(:))-ttime(faultpoints_x0_y11_z1))'./(4.*Dhx(posx).*Dhy(posy).*(Dhz(posz).^2)));%
    Vz_new(faultpoints_x0_y1_z1)=Vz_new(faultpoints_x0_y1_z1)+(Dt./(rho(faultpoints_x0_y1_z1))).*(F_tyz(t+max(ttime(:))-ttime(faultpoints_x0_y1_z1))'./(4.*Dhx(posx).*Dhy(posy).*(Dhz(posz).^2)));%



    % Higdon ABC
    wb=((1:N)-1)/N;
    %Boundary y=1
    beta=(1+(sqrt((lamda(1:N,:,:)+2.*mu(1:N,:,:))./rho(1:N,:,:))./sqrt(mu(1:N,:,:)./rho(1:N,:,:))))./2;
    r=(sqrt(mu(1:N,:,:)./rho(1:N,:,:)).*Dt)./Dhy(1:N);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vz_new(1:N,:,:)=(wb'.*Vz_new(1:N,:,:))+((1-wb').*(-qx.*Vz_new(2:N+1,:,:)-qt.*Vz_old_y1(1:N,:,:)-qxt.*Vz_old_y1(2:N+1,:,:)));     
    Vx_new(1:N,:,:)=(wb'.*Vx_new(1:N,:,:))+((1-wb').*(-qx.*Vx_new(2:N+1,:,:)-qt.*Vx_old_y1(1:N,:,:)-qxt.*Vx_old_y1(2:N+1,:,:))); 
    r=(sqrt((lamda(1:N,:,:)+2.*mu(1:N,:,:))./rho(1:N,:,:)).*Dt)./Dhy(1:N);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vy_new(1:N,:,:)=(wb'.*Vy_new(1:N,:,:))+((1-wb').*(-qx.*Vy_new(2:N+1,:,:)-qt.*Vy_old_y1(1:N,:,:)-qxt.*Vy_old_y1(2:N+1,:,:)));     
    clear beta r qx qt
    %Boundarty  x=1
    beta=(1+(sqrt((lamda(:,1:N,:)+2.*mu(:,1:N,:))./rho(:,1:N,:))./sqrt(mu(:,1:N,:)./rho(:,1:N,:))))./2;
    r=(sqrt(mu(:,1:N,:)./rho(:,1:N,:)).*Dt)./Dhx(1:N);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vz_new(:,1:N,:)=(wb.*Vz_new(:,1:N,:))+((1-wb).*(-qx.*Vz_new(:,2:N+1,:)-qt.*Vz_old_x1(:,1:N,:)-qxt.*Vz_old_x1(:,2:N+1,:)));   
    Vy_new(:,1:N,:)=(wb.*Vy_new(:,1:N,:))+((1-wb).*(-qx.*Vy_new(:,2:N+1,:)-qt.*Vy_old_x1(:,1:N,:)-qxt.*Vy_old_x1(:,2:N+1,:)));    
    r=(sqrt((lamda(:,1:N,:)+2.*mu(:,1:N,:))./rho(:,1:N,:)).*Dt)./Dhx(1:N);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vx_new(:,1:N,:)=(wb.*Vx_new(:,1:N,:))+((1-wb).*(-qx.*Vx_new(:,2:N+1,:)-qt.*Vx_old_x1(:,1:N,:)-qxt.*Vx_old_x1(:,2:N+1,:)));    
    clear beta r qx qt
    %Boundarty y=end
    beta=((1+(sqrt((lamda(end-(N-1):end,:,:)+2.*mu(end-(N-1):end,:,:))./rho(end-(N-1):end,:,:))./sqrt(mu(end-(N-1):end,:,:)./rho(end-(N-1):end,:,:))))./2);
    r=(sqrt(mu(end-(N-1):end,:,:)./rho(end-(N-1):end,:,:)).*Dt)./Dhy(end-(N-1):end);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vz_new(end-(N-1):end,:,:)=((1-wb').*Vz_new(end-(N-1):end,:,:))+(wb'.*(-qx.*Vz_new(end-N:end-1,:,:)-qt.*Vz_old_yend(end-(N-1):end,:,:)-qxt.*Vz_old_yend(end-N:end-1,:,:)));     
    Vx_new(end-(N-1):end,:,:)=((1-wb').*Vx_new(end-(N-1):end,:,:))+(wb'.*(-qx.*Vx_new(end-N:end-1,:,:)-qt.*Vx_old_yend(end-(N-1):end,:,:)-qxt.*Vx_old_yend(end-N:end-1,:,:)));       
    r=(sqrt((lamda(end-(N-1):end,:,:)+2.*mu(end-(N-1):end,:,:))./rho(end-(N-1):end,:,:)).*Dt)./Dhy(end-(N-1):end);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vy_new(end-(N-1):end,:,:)=((1-wb').*Vy_new(end-(N-1):end,:,:))+(wb'.*(-qx.*Vy_new(end-N:end-1,:,:)-qt.*Vy_old_yend(end-(N-1):end,:,:)-qxt.*Vy_old_yend(end-N:end-1,:,:)));  
    clear beta r qx qt
    %Boundarty x=end
    beta=(1+(sqrt((lamda(:,end-(N-1):end,:)+2.*mu(:,end-(N-1):end,:))./rho(:,end-(N-1):end,:))./sqrt(mu(:,end-(N-1):end,:)./rho(:,end-(N-1):end,:))))./2;
    r=(sqrt(mu(:,end-(N-1):end,:)./rho(:,end-(N-1):end,:)).*Dt)./Dhx(end-(N-1):end);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vz_new(:,end-(N-1):end,:)=((1-wb).*Vz_new(:,end-(N-1):end,:))+(wb.*(-qx.*Vz_new(:,end-N:end-1,:)-qt.*Vz_old_xend(:,end-(N-1):end,:)-qxt.*Vz_old_xend(:,end-N:end-1,:)));    
    Vy_new(:,end-(N-1):end,:)=((1-wb).*Vy_new(:,end-(N-1):end,:))+(wb.*(-qx.*Vy_new(:,end-N:end-1,:)-qt.*Vy_old_xend(:,end-(N-1):end,:)-qxt.*Vy_old_xend(:,end-N:end-1,:)));   
    r=(sqrt((lamda(:,end-(N-1):end,:)+2.*mu(:,end-(N-1):end,:))./rho(:,end-(N-1):end,:)).*Dt)./Dhx(end-(N-1):end);
    qx=((b.*(beta+r)-r))./((beta+r).*(1-b));qt=((b.*(beta+r)-beta))./((beta+r).*(1-b));qxt=(b./(b-1));
    Vx_new(:,end-(N-1):end,:)=((1-wb).*Vx_new(:,end-(N-1):end,:))+(wb.*(-qx.*Vx_new(:,end-N:end-1,:)-qt.*Vx_old_xend(:,end-(N-1):end,:)-qxt.*Vx_old_xend(:,end-N:end-1,:)));       
    clear beta r qx qt
    %Boundarty z=end
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

    % Calculate txx stress component
    txx(2:end-1,2:end-1,2:end-1)=Att(indtii_Q1,indtii_Q2,indtii_Q3).*(txx(2:end-1,2:end-1,2:end-1)+Dt.*(((lamda(2:end-1,2:end-1,2:end-1)+2.*mu(2:end-1,2:end-1,2:end-1)).*...
        ((Vx_new(2:end-1,2:end-1,2:end-1)-Vx_new(2:end-1,1:end-2,2:end-1))./Dhx(2:end-1)))+...
        (lamda(2:end-1,2:end-1,2:end-1)).*(((Vy_new(2:end-1,2:end-1,2:end-1)-Vy_new(1:end-2,2:end-1,2:end-1))./Dhy(2:end-1))+...
        ((Vz_new(2:end-1,2:end-1,2:end-1)-Vz_new(2:end-1,2:end-1,1:end-2))./Dhz(2:end-1)))));

    % Calculate tyy stress component
    tyy(2:end-1,2:end-1,2:end-1)=Att(indtii_Q1,indtii_Q2,indtii_Q3).*(tyy(2:end-1,2:end-1,2:end-1)+Dt.*(((lamda(2:end-1,2:end-1,2:end-1)+2.*mu(2:end-1,2:end-1,2:end-1)).*...
        ((Vy_new(2:end-1,2:end-1,2:end-1)-Vy_new(1:end-2,2:end-1,2:end-1))./Dhx(2:end-1)))+...
        (lamda(2:end-1,2:end-1,2:end-1)).*(((Vx_new(2:end-1,2:end-1,2:end-1)-Vx_new(2:end-1,1:end-2,2:end-1))./Dhy(2:end-1))+...
        ((Vz_new(2:end-1,2:end-1,2:end-1)-Vz_new(2:end-1,2:end-1,1:end-2))./Dhz(2:end-1)))));

    % Calculate tzz stress component
    tzz(2:end-1,2:end-1,2:end-1)=Att(indtii_Q1,indtii_Q2,indtii_Q3).*(tzz(2:end-1,2:end-1,2:end-1)+Dt.*(((lamda(2:end-1,2:end-1,2:end-1)+2.*mu(2:end-1,2:end-1,2:end-1)).*...
        ((Vz_new(2:end-1,2:end-1,2:end-1)-Vz_new(2:end-1,2:end-1,1:end-2))./Dhx(2:end-1)))+...    
        (lamda(2:end-1,2:end-1,2:end-1)).*(((Vx_new(2:end-1,2:end-1,2:end-1)-Vx_new(2:end-1,1:end-2,2:end-1))./Dhy(2:end-1))+...
        ((Vy_new(2:end-1,2:end-1,2:end-1)-Vy_new(1:end-2,2:end-1,2:end-1))./Dhz(2:end-1)))));

    % Calculate txy stress component
    txy(2:end,2:end,1:end-1)=Att(indtxy_Q1,indtxy_Q2,indtxy_Q3).*(txy(2:end,2:end,1:end-1)+...
        (mu(2:end,2:end,1:end-1).*Dt).*...
        (((Vx_new(2:end,1:end-1,1:end-1)-Vx_new(1:end-1,1:end-1,1:end-1))./Dhy(1:end-1))+...
        ((Vy_new(1:end-1,2:end,1:end-1)-Vy_new(1:end-1,1:end-1,1:end-1))./Dhx(1:end-1))));

    % Calculate txz stress component
    txz(1:end-1,2:end,2:end)=Att(indtxz_Q1,indtxz_Q2,indtxz_Q3).*(txz(1:end-1,2:end,2:end)+...
        (mu(1:end-1,2:end,2:end).*Dt).*...
        (((Vx_new(1:end-1,1:end-1,2:end)-Vx_new(1:end-1,1:end-1,1:end-1))./Dhz(1:end-1))+...
        ((Vz_new(1:end-1,2:end,1:end-1)-Vz_new(1:end-1,1:end-1,1:end-1))./Dhx(1:end-1))));

    % Calculate tyz stress component
    tyz(2:end,1:end-1,2:end)=Att(indtyz_Q1,indtyz_Q2,indtyz_Q3).*(tyz(2:end,1:end-1,2:end)+...
        (mu(2:end,1:end-1,2:end).*Dt).*...
        (((Vz_new(2:end,1:end-1,1:end-1)-Vz_new(1:end-1,1:end-1,1:end-1))./Dhy(1:end-1))+...
        ((Vy_new(1:end-1,1:end-1,2:end)-Vy_new(1:end-1,1:end-1,1:end-1))./Dhz(1:end-1))));


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