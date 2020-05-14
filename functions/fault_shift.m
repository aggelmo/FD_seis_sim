function [faultpoints_x0_y0_z0,faultpoints_x1_y0_z0, faultpoints_x1_y1_z0, ...
    faultpoints_x1_y11_z0,faultpoints_x0_y1_z0,faultpoints_x0_y11_z0,...
    faultpoints_x0_y0_z1,faultpoints_x0_y0_z11,faultpoints_x1_y0_z11,...
    faultpoints_x11_y0_z0,faultpoints_x11_y1_z0, faultpoints_x0_y1_z11,...
    faultpoints_x0_y1_z1,faultpoints_x1_y0_z1,faultpoints_x11_y0_z1,...
    faultpoints_x0_y11_z1]=fault_shift(faultpoints)
    % Function to obtain the node indices for the fault source relative to
    % each velocity component. 
    %faultpoints for Vx
    faultpoints_x0_y0_z0=find(faultpoints==1);
    faultpoints_x1_y0_z0=find(circshift(faultpoints,1,2)==1);
    temp=circshift(faultpoints,1,2);
    faultpoints_x1_y1_z0=find(circshift(temp,1,1)==1);
    faultpoints_x1_y11_z0=find(circshift(temp,-1,1)==1);
    faultpoints_x0_y1_z0=find(circshift(faultpoints,1,1)==1);
    faultpoints_x0_y11_z0=find(circshift(faultpoints,-1,1)==1);
    faultpoints_x0_y0_z1=find(circshift(faultpoints,1,3)==1);
    faultpoints_x0_y0_z11=find(circshift(faultpoints,-1,3)==1);
    temp=circshift(faultpoints,1,3);
    faultpoints_x1_y0_z1=find(circshift(temp,1,2)==1);
    faultpoints_x1_y0_z11=find(circshift(temp,-1,2)==1);
    
    
    %extra faultpoints for Vy
    faultpoints_x11_y0_z0=find(circshift(faultpoints,-1,2)==1);
    temp=circshift(faultpoints,1,1);
    faultpoints_x11_y1_z0=find(circshift(temp,-1,2)==1);
    faultpoints_x0_y1_z11=find(circshift(temp,-1,3)==1);
    faultpoints_x0_y1_z1=find(circshift(temp,1,3)==1);
    
    %extra faultpoints for Vz
    temp=circshift(faultpoints,-1,1);
    faultpoints_x0_y11_z1=find(circshift(temp,1,3)==1);
    temp=circshift(faultpoints,-1,2);
    faultpoints_x11_y0_z1=find(circshift(temp,1,3)==1);
end


