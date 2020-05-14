%% For surface pgv
load('example2_input.mat')
clear waveformx waveformy waveformz
ny=201;
nx=201;
nt=500;
waveformx=single(zeros(ny,nz,40,nt));waveformy=waveformx;waveformz=waveformx;
for i=1:nt
    i
    load(['Vx_' num2str(i) '.mat'])
    waveformx(:,:,:,i)=single(wfield(:,:,1:40));
    load(['Vy_' num2str(i) '.mat'])
    waveformy(:,:,:,i)=single(wfield(:,:,1:40));
    load(['Vz_' num2str(i) '.mat'])
    waveformz(:,:,:,i)=single(wfield(:,:,1:40));
end
clear test11


test=sqrt(waveformx.^2+waveformy.^2+waveformz.^2);
test2=max(abs(test),[],4);
[~,~,ZZ]=meshgrid(X,Y,Z(1:40)+2.1);
test3=abs(ZZ-elev_data);
[temp,test3]=min(test3,[],3);
[X_topo,Y_topo]=meshgrid(X,Y);

for i=1:nx
    for j=1:ny
        test4(i,j)=test2(i,j,test3(i,j));
    end
end
%%
surf(X_topo,Y_topo,elev_data,test4./max(abs(test4(:))))
caxis([0 0.1])
%% for vertical slice 1
clear waveformx waveformy waveformz XX YY ZZ test test2 test3 test5
ny=201;
nx=201;
nz=122;
nt=449;
waveformx=single(zeros(ny,nz,nt));waveformy=waveformx;waveformz=waveformx;
for i=1:nt
    i
    load(['Vx_' num2str(i) '.mat'])
    waveformx(:,:,i)=single(wfield(100,:,:));
    load(['Vy_' num2str(i) '.mat'])
    waveformy(:,:,i)=single(wfield(100,:,:));
    load(['Vz_' num2str(i) '.mat'])
    waveformz(:,:,i)=single(wfield(100,:,:));
end
waveformx=reshape(waveformx,ny,nz,nt);
waveformy=reshape(waveformx,ny,nz,nt);
waveformz=reshape(waveformx,ny,nz,nt);
clear test11

%%
[XX,ZZ]=meshgrid(X,2.1:-0.1:-10);

test=sqrt(waveformx.^2+waveformy.^2+waveformz.^2);
test2=max(test,[],3);test2(test2==0)=nan;
contourf(XX(1:end-10,11:end-10),ZZ(1:end-10,11:end-10),test2(11:end-10,1:end-10)'./max(abs(test2(:)))',50,'linestyle','none');
caxis([-0 0.3])
xlabel('X (km)')
ylabel('Z (km)')
axis equal
hold on
Vp2=Vp(:,101,:);Vp2=reshape(Vp2,ny,nz-21);
contour(XX(22:end-10,11:end-10),ZZ(22:end-10,11:end-10),Vp2(11:end-10,1:end-10)',5,'black');
hold off
