%plot the initial config of net
%------------------------------
clear all;
addpath(genpath('../'))
sstep=30;

drawnow
plot_configuration

fileid='f0id1';
a=4;
b=4;
lambda=0.4;
Ly=60;
%%%%%%%%%%%%%----read out system size Lx-------------
fileID=fopen('Lx.txt','r');

sizeA =[1 Inf];

Lx=fscanf(fileID, '%d', sizeA);

fclose(fileID);
%%%-------reading out bonds--------------------

fileID=fopen([fileid,'bonds.txt'],'r');

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

fileID=fopen(['xc',fileid,'.txt'],'r');

sizeA =[2 Inf];

temp=fscanf(fileID,'%f',sizeA);

for k=1:size(temp,2)
   xc_list(k,1)=temp(1,k); 
   xc_list(k,2)=temp(2,k);
end

fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step=size(xc_list,1);

for i=1:step
    i
    clf;
    %--------read new node coordinates---------------------
    fileID=fopen(['nodes',num2str(i),fileid,'.txt'],'r');
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
    t=mod(mod(theta_list(i),6.28)+6.28,6.28);
    text(40,70,['\theta=',num2str(t*180/3.14,'%10.3e')]);
    
    for j=1:length(Bondlist)
        G1=Bondlist(j).G1;
        G2=Bondlist(j).G2;
        
        x1=Superlist(G1).x1;
        y1=Superlist(G1).y1;

        x2=Superlist(G2).x1;
        y2=Superlist(G2).y1;

        
        j1=(G1-1)/Ly;
        j2=(G2-1)/Ly;
        
        if (j1>j2) temp=j2;
        else temp=j1;
        end
        
        if (temp>Ly/2) line([x1,x2],[y1,y2],'LineWidth',2);
        else line([x1,x2],[y1,y2],'LineWidth',2);
        end
        hold on
    end
    
%     if i<=sstep
%         ellipse(a*exp(-lambda*i/sstep),b,theta_list(i),xc_list(i,1),xc_list(i,2),'k');
%     else
%         ellipse(a*exp(-lambda*(i-20)/sstep),b,theta_list(i),xc_list(i,1),xc_list(i,2),'k');
%     end

    ellipse(a*exp(-lambda*mod(i-1,sstep)/sstep),b,theta_list(i),xc_list(i,1),xc_list(i,2),'k');
    
    vector=10*[cos(theta_list(i)),sin(theta_list(i))];
    quiver(xc_list(i,1), xc_list(i,2), vector(1), vector(2),'LineWidth',2,'Color','Cyan');
    
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for j=1:size(Adhesion,1)
       if ifremove(j)==1
           continue;
       end
        
       x=Adhesion(j,1)+xc_list(i,1);
       y=Adhesion(j,2)+xc_list(i,2);
       plot(x,y, '.r','MarkerSize',20)
       x0=Superlist(AdhesionID(j,1)).x1+Adhesion_shift(j,1);
       y0=Superlist(AdhesionID(j,1)).y1+Adhesion_shift(j,2);
       
%        if (y>40||y0>40)
%            fprintf('x=%f, y=%f, x0=%f, y0=%f, j=%d\n',x,y,x0,y0,j)
%        end
       
       
       line([x,x0],[y,y0],'LineWidth',2,'Color','r');
    end
    
%     xlim([xc_list(i,1)-15,xc_list(i,1)+15])
%     ylim([xc_list(i,2)-15,xc_list(i,2)+15])
    
%     xlim([xc_list(i,1)-15,xc_list(i,1)+15])
    axis off
    xlim([20,60])
    ylim([10,44.64])
    plot(linspace(0,80,40),ones(40)*25.98,'k--','LineWidth',2)

    drawnow
    M(i)=getframe(gcf);
   
    
    
    clear Adhesion AdhesionID Adhesion_shift ifremove
end

movie2avi(M,'cell1Movie.avi','fps',6);