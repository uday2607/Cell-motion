function y= cal_displacement(r_seq, width, Superlist0, Superlist, xc, yc)

    y= zeros(size(r_seq));
    n= zeros(size(r_seq));
    
    for i=1: length(Superlist0)
        
        x0=Superlist0(i).x1;
        y0=Superlist0(i).y1;
        
        x1=Superlist(i).x1;
        y1=Superlist(i).y1;
        
        dx= x0-xc;
        dy= y0-yc;

        dl= sqrt(dx*dx+ dy*dy);

        dx= dx/dl;
        dy= dy/dl;

        deltaX= x1- x0;
        deltaY= y1- y0;

        dr= deltaX*dx + deltaY*dy;

        
        for j= 1: length(r_seq)
            
           if (dl> r_seq(j)-width/2) && (dl<= r_seq(j)+width/2)
              
              %%%%%% 
%               fprintf('x0=%d y0=%d x1=%d y1=%d dx=%d dy=%d dl=%d\n',x0, y0, x1, y1, dx, dy, dl) 
              %%%%% 
               
              y(j)= y(j)+dr;
              n(j)= n(j)+1;
              break;
           end
            
            
        end
        
       
        
      
        
    
    end
    
    y= y./n;