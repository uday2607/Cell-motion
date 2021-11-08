clear all; clc; close all; 

plot_configuration;
addpath(genpath('../'));
label='f0id1';
fileid='f0id1';

i=1;
sstep=30;

Lx=60;


a=6;
b=4;
lambda=0.4;
% xc=60*(3/4-1/5);
% yc=60/4*sqrt(3);
xc=Lx*3/4-10;
yc=Lx/4*sqrt(3);

% xc2=60*(3/4+1/5);
% yc2=60/4*sqrt(3);
xc=Lx*3/4+10;
yc=Lx/4*sqrt(3);

%%%%%%%%%%%%%----read out system size Lx-------------
fileID=fopen('Lx.txt','r');

sizeA =[1 Inf];

Lx=fscanf(fileID, '%d', sizeA);

fclose(fileID);

%%%-------reading out bonds--------------------

fileID=fopen([label,'bonds.txt'],'r');

sizeA = [2 Inf];

temp=fscanf(fileID,'%d %d',sizeA);

Bondlist=struct();
hold on
for i=1:length(temp)
   G1=temp(1,i);
   G2=temp(2,i);

   Bondlist(i).G1=G1;
   Bondlist(i).G2=G2;
end
% 

fclose(fileID);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID=fopen(['theta',fileid,'.txt'],'r');

sizeA =[1 Inf];

temp=fscanf(fileID,'%f',sizeA);

for k=1:size(temp,2)
   theta_list(k)=temp(1,k);   
end

fclose(fileID);

fileID=fopen(['theta2',fileid,'.txt'],'r');

sizeA =[1 Inf];

temp=fscanf(fileID,'%f',sizeA);

for k=1:size(temp,2)
   theta_list2(k)=temp(1,k);   
end

fileID=fopen(['phase', fileid, '.txt'], 'r');

sizeA =[1 Inf];

temp=fscanf(fileID,'%f',sizeA);

for k=1:size(temp,2)
   phase_list(k)=temp(1,k);   
end

fclose(fileID);

%%--------------------------------------------

fildID=fopen(['phase2', fileid, '.txt'], 'r');

sizeA =[1 Inf];

temp=fscanf(fileID,'%f',sizeA);

for k=1:size(temp,2)
   phase_list2(k)=temp(1,k);   
end

fclose(fileID);

%%--------------------------------------------

fildID=fopen(['period', fileid, '.txt'], 'r');

sizeA =[1 Inf];

temp=fscanf(fileID,'%f',sizeA);

for k=1:size(temp,2)
   period_list(k)=temp(1,k);   
end

fclose(fileID);

%%--------------------------------------------

fildID=fopen(['period2', fileid, '.txt'], 'r');

sizeA =[1 Inf];

temp=fscanf(fileID,'%f',sizeA);

for k=1:size(temp,2)
   period_list2(k)=temp(1,k);   
end

fclose(fileID);

%%--------------------------------------------

fileID=fopen(['xc',fileid,'.txt'],'r');

sizeA =[2 Inf];

temp=fscanf(fileID,'%f',sizeA);

for k=1:size(temp,2)
   xc_list(k,1)=temp(1,k); 
   xc_list(k,2)=temp(2,k);
end

fclose(fileID);

fileID=fopen(['xc2',fileid,'.txt'],'r');

sizeA =[2 Inf];

temp=fscanf(fileID,'%f',sizeA);

for k=1:size(temp,2)
   xc_list2(k,1)=temp(1,k); 
   xc_list2(k,2)=temp(2,k);
end

fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=size(xc_list,1);

for i=1:step
    i
    cla;
    %--------read new node coordinates---------------------
    fileID=fopen(['nodes',num2str(i),label,'.txt'],'r');
    %readout filament infomation from .txt file

    sizeA = [3 Inf];

    temp=fscanf(fileID,'%f %f %d',sizeA);

    Superlist=struct();

    for j=1:length(temp) 
        index=temp(3,j);

        Superlist(index).x1=temp(1,j);
        Superlist(index).y1=temp(2,j);
    end

    fclose(fileID);
    %------------------------------------------------------
    
    for j=1:length(Bondlist)
        G1=Bondlist(j).G1;
        G2=Bondlist(j).G2;
        
        x1=Superlist(G1).x1;
        y1=Superlist(G1).y1;

        x2=Superlist(G2).x1;
        y2=Superlist(G2).y1;

        line([x1,x2],[y1,y2],'LineWidth',2);
        hold on
    end
    
    phase1=phase_list(i);
    phase2=phase_list2(i);
    period1=period_list(i);
    period2=period_list2(i);
    
    
    ellipse(a*exp(-lambda*phase1/period1),b,theta_list(i),xc_list(i,1),xc_list(i,2),'k');
    ellipse(a*exp(-lambda*phase2/period2),b,theta_list2(i),xc_list2(i,1),xc_list2(i,2),'k');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fileID=fopen(['Adhesion_step',num2str(i),'.txt'],'r');
    
    sizeA= [7 Inf];
    temp=fscanf(fileID,'%f %f %d %d %f %f %d',sizeA);

    for j=1:size(temp,2)
        Adhesion(j,1)=temp(1,j);
        Adhesion(j,2)=temp(2,j);
        AdhesionID(j,1)=temp(3,j);
        AdhesionID(j,2)=temp(4,j);
        Adhesion_shift(j,1)=temp(5,j);
        Adhesion_shift(j,2)=temp(6,j);
        ifremove(j)=temp(7,j);
    end
    fclose(fileID);
    
    fileID=fopen(['Adhesion2_step',num2str(i),'.txt'],'r');
    
    sizeA= [7 Inf];
    temp=fscanf(fileID,'%f %f %d %d %f %f %d',sizeA);

    for j=1:size(temp,2)
        Adhesion2(j,1)=temp(1,j);
        Adhesion2(j,2)=temp(2,j);
        AdhesionID2(j,1)=temp(3,j);
        AdhesionID2(j,2)=temp(4,j);
        Adhesion_shift2(j,1)=temp(5,j);
        Adhesion_shift2(j,2)=temp(6,j);
        ifremove2(j)=temp(7,j);
    end
    fclose(fileID);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
    for j=1:size(Adhesion,1)
       if ifremove(j)==1
           continue;
       end
        
       x=Adhesion(j,1)+xc_list(i,1);
       y=Adhesion(j,2)+xc_list(i,2);
       plot(x,y, '.r','MarkerSize',20)
       x0=Superlist(AdhesionID(j,1)).x1+Adhesion_shift(j,1);
       y0=Superlist(AdhesionID(j,1)).y1+Adhesion_shift(j,2);
       
       line([x,x0],[y,y0],'LineWidth',2,'Color','r');
    end
    
    for j=1:size(Adhesion2,1)
       if ifremove2(j)==1
           continue;
       end
        
       x=Adhesion2(j,1)+xc_list2(i,1);
       y=Adhesion2(j,2)+xc_list2(i,2);
       plot(x,y, '.r','MarkerSize',20)
       x0=Superlist(AdhesionID2(j,1)).x1+Adhesion_shift2(j,1);
       y0=Superlist(AdhesionID2(j,1)).y1+Adhesion_shift2(j,2);
       line([x,x0],[y,y0],'LineWidth',2,'Color','r');
    end
    
    xlim([-1,90])
    ylim([-1,55])
    
    
    drawnow
    M(i)=getframe(gcf);
    
    clear Adhesion AdhesionID Adhesion_shift ifremove Adhesion2 AdhesionID2 Adhesion_shift2 ifremove2
end

movie2avi(M,'cell3Movie.avi','fps',6);