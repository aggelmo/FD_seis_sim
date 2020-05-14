function [faultpoints, ttime, nfact]=fin_fault(posx,posy,posz,X,Y,Z,flength,fwidth,strike,dip,ix,iz,vr,Dt,shiftval)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that obtains the nodes that represent the fault surface within
% the model grid and calculates the rupture travel times from a source node
% towards all other fault nodes.
% posx, posy, posz are the X, Y and Z index of the fault upper left edge.
% X, Y, Z are 3D matrices with the cartesian coordinates of the model grid.
% flength and fwidth are the fault length and width respectively. 
% strike and Dip are the fault strike and dip.
% ix and iz is the location of the initial rupture point within the fault
% (Values: 0-1, e.g. value 0.5 means thet the rupture source is in the
% middle of the fault).
% vr is the rupture velocity.
% Dt is the timestep of the simulations.
% shiftval is the fault "volume" half-width.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check input
if length(posx)>1 || length(posz)>1 || length(flength)>1 || length(fwidth)>1 ...
        || length(strike)>1 || length(dip)>1 || length(ix)>1 || length(iz)>1 ...
        || length(vr)>1 || length(Dt)>1
    disp('Incorrect input format')
    return
end
if ix>1 || ix<0 || iz>1 || iz<0
    disp('ix and iz must  be between 0 and 1')
    return
end

% Initialize a matrix
V=zeros(size(X));

% Form a horizontal surface with the fault dimensions
% Find the X node index of the last surface point in the X direction
[~,indx2]=min(abs(X(1,:,1)-(X(1,posx,1)+flength))); 
% Do the same for the Y node index of the last surface point in the Y
% direction
[~,indy2]=min(abs(Y(:,1,1)-(Y(posy,1,1)+fwidth))); 

% Set the values of this surface to 1
V(posy:indy2,posx:indx2,posz)=1;

% Create vectors with the 3D coordinates and values of matrix V, setting
% the central point as the fault uppermost left edge, defined by posx,posy
% and posz
newX=X(:)-X(1,posx,1);
newY=Y(:)-Y(posy,1,1);
newZ=Z(:)-Z(1,1,posz); %depth must be negative!!!!
fvals=V(:);

% Remove values that do not lie on the initial fault surface
fspace(:,1)=newX(fvals==1);
fspace(:,2)=newY(fvals==1);
fspace(:,3)=newZ(fvals==1);

% Rotate surface relative to fault dip
fspace_new(:,2)=fspace(:,2).*cosd(dip)-fspace(:,3).*sind(dip);
fspace_new(:,3)=fspace(:,2).*sind(dip)+fspace(:,3).*cosd(dip);
fspace_new(:,1)=fspace(:,1);

% Rotate resulting surfrace relative to fault strike
fspace_new2(:,1)=fspace_new(:,1).*cosd(-strike+90)-fspace_new(:,2).*sind(-strike+90);
fspace_new2(:,2)=fspace_new(:,1).*sind(-strike+90)+fspace_new(:,2).*cosd(-strike+90);
fspace_new2(:,3)=fspace_new(:,3);

% Shift resulting matrix back in to the original coordinates
temp(:,1)=fspace_new2(:,1)+X(1,posx,1);temp(:,2)=fspace_new2(:,2)+Y(posy,1,1);temp(:,3)=fspace_new2(:,3)+Z(1,1,posz);temp(:,4)=ones(size(temp(:,1)));

% Create two surfaces shifted relative to the obtained fault surface. These
% surface seperate the fault "volume" from the outer nodes.
d_fault(:,1)=shiftval.*cosd(strike)-shiftval.*sind(strike);
d_fault(:,2)=shiftval.*sind(strike)+shiftval.*cosd(strike);
d_fault(:,3)=0;

temp2(:,1)=temp(:,1)+d_fault(1,1);temp2(:,2)=temp(:,2)+d_fault(1,2);temp2(:,3)=temp(:,3);temp2(:,4)=temp(:,4)-1;
temp3(:,1)=temp(:,1)-d_fault(1,1);temp3(:,2)=temp(:,2)-d_fault(1,2);temp3(:,3)=temp(:,3);temp3(:,4)=temp(:,4)-1;
temp4=[temp;temp2;temp3];

% Obtain values that contain the faultposz nodes in the initial model grid
faultpoints=single(griddata(double(temp4(:,1)),double(temp4(:,2)),double(temp4(:,3)),double(temp4(:,4)),double(X),double(Y),double(Z),'natural'));
faultpoints(isnan(faultpoints))=0;

% Find point of surface contained within the model grid to set as rupture initiation
fpoints_temp(:,1)=X(faultpoints>0.5);fpoints_temp(:,2)=Y(faultpoints>0.5);fpoints_temp(:,3)=Z(faultpoints>0.5);
source_point(1,1)=min(fpoints_temp(:,1))+(max(fpoints_temp(:,1))-min(fpoints_temp(:,1)))./2;
source_point(1,2)=min(fpoints_temp(:,2))+(max(fpoints_temp(:,2))-min(fpoints_temp(:,2))).*ix;
source_point(1,3)=min(fpoints_temp(:,3))+(max(fpoints_temp(:,3))-min(fpoints_temp(:,3))).*iz;

% Calculate rupture travel time from the source node towards all nodes of
% the fault.
% Calculate absolute distance between each node and the source point
dist=sqrt(abs(X-source_point(1,1)).^2+abs(Y-source_point(1,2)).^2+abs(Z-source_point(1,3)).^2);


% Calculate the travel time
ttime=dist./vr;

% Round and convert seconds to time step index
% ttime=(ttime./(Dt));
ttime=round(ttime./(Dt));

faultpoints(faultpoints>0.5)=1;faultpoints(faultpoints<=0.5)=0;

% Zero traveltime values outside the fault surface
% ttime=ttime.*faultpoints;

% Calculate normalization factor relative to total no of subsources
nfact=1/length(temp(:,4));
end