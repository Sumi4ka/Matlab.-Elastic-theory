classdef solver
    properties
        r          %Distribution in the radial direction
        z          %Distribution in the axial direction
        N          %Count of nodes in the axial direction
        M          %Count of nodes in the radial direction
        w          %Matrix of Axial displacement
        v          %Matrix of Radial displacement
        z_left_v   %Condition for Radial displacement in the left side
        z_left_w   %Condition for Axial displacement in the left side
        z_right_v  %Condition for Radial displacement in the right side
        z_right_w  %Condition for Axial displacement in the right side
        z_bot_v    %Condition for Radial displacement in the bot side
        z_bot_w    %Condition for Axial displacement in the bot side
        z_top_v    %Condition for Radial displacement in the top side
        z_top_w    %Condition for Axial displacement in the top side
        sigmazr    %Matrix of sigmzr distribution
        sigmazz    %Matrix of sigmzz distribution
        sigmarr    %Matrix of sigmrr distribution
        sigmaff    %Matrix of sigmff distribution
    end    
    methods
        %Constructor
        function obj=solver(N,M,mesh,zone)
            obj.N=N;
            obj.M=M;
            obj.v=zeros(N,M);
            obj.w=zeros(N,M);
            obj.r=mesh.r;
            obj.z=mesh.z;
            obj.z_left_v=zone.z_left_v;
            obj.z_left_w=zone.z_left_w;
            obj.z_right_v=zone.z_right_v;
            obj.z_right_w=zone.z_right_w;
            obj.z_bot_v=zone.z_bot_v;
            obj.z_bot_w=zone.z_bot_w;
            obj.z_top_v=zone.z_top_v;
            obj.z_top_w=zone.z_top_w;
            obj.sigmazr=zeros(N,M);
            obj.sigmazz=zeros(N,M);
            obj.sigmarr=zeros(N,M);
            obj.sigmaff=zeros(N,M);
        end
        %Solver
        function obj=solve(obj,nu,E,eps)
            E=E*10^9;                   %Young's Modulus
            b=2*(1-nu)/(1-2*nu);
            mu=E/2/(1+nu);              %First Lame constant 
            lambda=nu*E/(1+nu)/(1-2*nu);%Second Lame constant 
            c=lambda/(lambda+2*mu);
            t=int16(0);                 %Number of iteration
            eps=10^(eps);               %Accuracy
            con=2*eps;
            tic
            while (con>eps)
                
                t=t+1;
                con=0;
                %Internal zone
                for k=2:obj.N-1
                    for j=2:obj.M-1
                        re=(obj.r(j-1)+obj.r(j))/2;
                        rw=(obj.r(j)+obj.r(j+1))/2;
                        sp1=obj.r(j)^2*(b-1)*(obj.w(k+1,j+1)-obj.w(k-1,j+1)-obj.w(k+1,j-1)+obj.w(k-1,j-1))/(obj.r(j+1)-obj.r(j-1))/(obj.z(k+1)-obj.z(k-1));
                        sp2=obj.r(j)^2*(b-1)*(obj.v(k+1,j+1)-obj.v(k-1,j+1)-obj.v(k+1,j-1)+obj.v(k-1,j-1))/(obj.r(j+1)-obj.r(j-1))/(obj.z(k+1)-obj.z(k-1));
                        sp2=sp2+obj.r(j)*(b-1)*(obj.v(k+1,j)-obj.v(k-1,j))/(obj.z(k+1)-obj.z(k-1));
                        ae=re*obj.r(j)/(obj.r(j)-obj.r(j-1))/(rw-re);
                        aw=rw*obj.r(j)/(obj.r(j+1)-obj.r(j))/(rw-re);
                        at=2*obj.r(j)^2/(obj.z(k)-obj.z(k-1))/(obj.z(k+1)-obj.z(k-1));
                        ab=2*obj.r(j)^2/(obj.z(k+1)-obj.z(k))/(obj.z(k+1)-obj.z(k-1));
                        ap=ae+aw+at+ab;
                        v1=obj.v(k,j);
                        w1=obj.w(k,j);
                        obj.v(k,j)=(b*ae*obj.v(k,j-1)+b*aw*obj.v(k,j+1)+at*obj.v(k-1,j)+ab*obj.v(k+1,j)+sp1)/(b+ap+ae*(b-1)+aw*(b-1));
                        obj.w(k,j)=(ae*obj.w(k,j-1)+aw*obj.w(k,j+1)+b*at*obj.w(k-1,j)+b*ab*obj.w(k+1,j)+sp2)/(ap+at*(b-1)+ab*(b-1));
                        con=con+(obj.v(k,j)-v1)^2+(obj.w(k,j)-w1)^2;
                    end
                end
                %External zone
                %Top_side
                for j=2:obj.M-1
                    v1=obj.v(obj.N,j);
                    w1=obj.w(obj.N,j);
                    if (obj.z_top_w(j)~=0)
                        if (obj.z_top_w(j)==-1)
                            obj.w(obj.N,j)=obj.w(obj.N-1,j)-c*(obj.z(obj.N)-obj.z(obj.N-1))/(obj.r(j+1)-obj.r(j-1))*(obj.v(obj.N,j+1)-obj.v(obj.N,j-1))-c*(obj.z(obj.N)-obj.z(obj.N-1))*obj.v(obj.N,j)/obj.r(j);
                        else
                            obj.w(obj.N,j)=(-1)*obj.z_top_w(j);
                        end
                    end
                    if (obj.z_top_v(j)==-1)
                        obj.v(obj.N,j)=obj.v(obj.N-1,j)-(obj.z(obj.N)-obj.z(obj.N-1))/(obj.r(j+1)-obj.r(j-1))*(obj.w(obj.N,j+1)-obj.w(obj.N,j-1));
                    end
                    con=con+(obj.v(obj.N,j)-v1)^2+(obj.w(obj.N,j)-w1)^2;
                end
                %Bot_side
                for j=2:obj.M-1
                    v1=obj.v(1,j);
                    if (obj.z_bot_v(j)==-1)
                        obj.v(1,j)=obj.v(2,j)+(obj.z(2)-obj.z(1))/(obj.r(j+1)-obj.r(j-1))*(obj.w(1,j+1)-obj.w(1,j-1));
                        con=con+(obj.v(1,j)-v1)^2;
                    end
                end
                %Right_bot_angle_point
                v1=obj.v(1,obj.M);
                if (obj.z_bot_v(obj.M)~=0)
                    if (obj.z_bot_v(obj.M)==-1)
                    obj.v(1,obj.M)=obj.v(1,obj.M-1)-(obj.r(obj.M)-obj.r(obj.M-1))*(obj.v(1,obj.M-1)/obj.r(obj.M-1)+1/c*(obj.w(2,obj.M)-obj.w(1,obj.M))/(obj.z(2)-obj.z(1)));
                    else
                    obj.v(1,obj.M)=(-1)*obj.z_bot_v(obj.M);
                    end
                end
                con=con+(obj.v(1,obj.M)-v1)^2;  
                %Right_side
                for k=2:obj.N-1
                    v1=obj.v(k,obj.M);
                    w1=obj.w(k,obj.M);
                    if (obj.z_right_v(k)~=0)
                        if (obj.z_right_v(k)==-1)
                            obj.v(k,obj.M)=obj.v(k,obj.M-1)-c*(obj.r(obj.M)-obj.r(obj.M-1))*((obj.w(k+1,obj.M)-obj.w(k-1,obj.M))/(obj.z(k+1)-obj.z(k-1))+obj.v(k,obj.M)/obj.r(obj.M));
                        else
                            obj.v(k,obj.M)=(-1)*obj.z_right_v(k);
                        end
                    end
                    if (obj.z_right_w(k)==-1)
                        obj.w(k,obj.M)=obj.w(k,obj.M-1)-(obj.r(obj.M)-obj.r(obj.M-1))*(obj.v(k+1,obj.M)-obj.v(k-1,obj.M))/(obj.z(k+1)-obj.z(k-1));
                    end
                    con=con+(obj.v(k,obj.M)-v1)^2+(obj.w(k,obj.M)-w1)^2;
                end
                %Right angle point
                x=zeros(1,2);
                %A=[1/(obj.r(obj.M)-obj.r(obj.M-1)), 1/(obj.z(obj.N)-obj.z(obj.N-1)); c/(obj.z(obj.N)-obj.z(obj.N-1)), 1/(obj.r(obj.M)-obj.r(obj.M-1))];
                %B=[obj.w(obj.N,obj.M-1)/(obj.r(obj.M)-obj.r(obj.M-1))+obj.v(obj.N-1,obj.M)/(obj.z(obj.N)-obj.z(obj.N-1)); obj.v(obj.N,obj.M-1)/(obj.r(obj.M)-obj.r(obj.M-1))+c*obj.w(obj.N-1,obj.M)/(obj.z(obj.N)-obj.z(obj.N-1))];
                %x=A\B;
                x(1)=obj.w(obj.N,obj.M-1)-(obj.r(obj.M)-obj.r(obj.M-1))*(obj.v(obj.N,obj.M)-obj.v(obj.N-1,obj.M))/(obj.z(obj.N)-obj.z(obj.N-1));
                x(2)=obj.v(obj.N,obj.M-1)-c*(obj.r(obj.M)-obj.r(obj.M-1))*((obj.w(obj.N,obj.M)-obj.w(obj.N-1,obj.M))/(obj.z(obj.N)-obj.z(obj.N-1))+obj.v(obj.N,obj.M)/obj.r(obj.M));
                v1=obj.v(obj.N,obj.M);
                w1=obj.w(obj.N,obj.M);
                if (obj.z_right_w(obj.N)~=0)
                    if (obj.z_right_w(obj.N)==-1)
                        obj.w(obj.N,obj.M)=x(1);
                    else
                        obj.w(obj.N,obj.M)=(-1)*obj.z_right_w(obj.N);
                    end
                end
                if (obj.z_right_v(obj.N)~=0)
                    if (obj.z_right_v(obj.N)==-1)
                        obj.v(obj.N,obj.M)=x(2);
                    else
                        obj.v(obj.N,obj.M)=(-1)*obj.z_right_v(obj.N);
                    end
                end
                con=con+(obj.v(obj.N,obj.M)-v1)^2+(obj.w(obj.N,obj.M)-w1)^2;
                %Left_side
                if (obj.r(1)~=0)
                    for k=2:obj.N
                        obj.v(k,1)=0;
                        obj.w(k,1)=obj.w(k,2);
                        con=con+(obj.v(k,1)-v1)^2+(obj.w(k,1)-w1)^2;
                    end
                else
                    for k=2:obj.N-1
                        v1=obj.v(k,1);
                        w1=obj.w(k,1);
                        if (obj.z_left_v(k)~=0)
                            if (obj.z_left_v(k)==-1)
                                obj.v(k,1)=obj.v(k,2)-c*(obj.r(1)-obj.r(2))*((obj.w(k+1,1)-obj.w(k-1,1))/(obj.z(k+1)-obj.z(k-1))+obj.v(k,1)/obj.r(1));
                            else
                                obj.v(k,1)=(-1)*obj.z_left_v(k);
                            end
                        end
                        if (obj.z_left_w(k)==-1)
                            obj.w(k,1)=obj.w(k,2)-(obj.r(1)-obj.r(2))*(obj.v(k+1,1)-obj.v(k-1,1))/(obj.z(k+1)-obj.z(k-1));
                        end
                        con=con+(obj.v(k,1)-v1)^2+(obj.w(k,1)-w1)^2;
                    end
                    %Left angle point
                    x=zeros(1,2);
                    %A=[1/(obj.r(1)-obj.r(2)), 1/(obj.z(obj.N)-obj.z(obj.N-1)); c/(obj.z(obj.N)-obj.z(obj.N-1)), 1/(obj.r(1)-obj.r(2))+c/obj.r(1)];
                    %B=[obj.w(obj.N,2)/(obj.r(1)-obj.r(2))+obj.v(obj.N-1,1)/(obj.z(obj.N)-obj.z(obj.N-1)); obj.v(obj.N,2)/(obj.r(1)-obj.r(2))+c*obj.w(obj.N-1,1)/(obj.z(obj.N)-obj.z(obj.N-1))];
                    %x=A\B;
                    x(1)=obj.w(obj.N,2)-(obj.r(1)-obj.r(2))*(obj.v(obj.N,1)-obj.v(obj.N-1,1))/(obj.z(obj.N)-obj.z(obj.N-1));
                    x(2)=obj.v(obj.N,2)-c*(obj.r(1)-obj.r(2))*((obj.w(obj.N,1)-obj.w(obj.N-1,1))/(obj.z(obj.N)-obj.z(obj.N-1))+obj.v(obj.N,1)/obj.r(1));
                    v1=obj.v(obj.N,1);
                    w1=obj.w(obj.N,1);
                    if (obj.z_left_w(obj.N)~=0)
                        if (obj.z_left_w(obj.N)==-1)
                            obj.w(obj.N,1)=x(1);
                        else
                            obj.w(obj.N,1)=(-1)*obj.z_left_w(obj.N);
                        end
                    end
                    if (obj.z_left_v(obj.N)~=0)
                        if (obj.z_left_v(obj.N)==-1)
                            obj.v(obj.N,1)=x(2);
                        else
                            obj.v(obj.N,1)=(-1)*obj.z_left_v(obj.N);
                        end
                    end
                    con=con+(obj.v(obj.N,1)-v1)^2+(obj.w(obj.N,1)-w1)^2;
                end
                con=sqrt(con);
                if (~(mod(t,1)))
                    con
                end
                %if (t>obj.N*obj.M)
                %    con=0;
                %end
            end
            %-----------------------------
            %Calculation of stress-strain matrices 
            for k=2:obj.N-1
                for j=2:obj.M-1
                    obj.sigmazr(k,j)=mu*((obj.w(k,j+1)-obj.w(k,j-1))/(obj.r(j+1)-obj.r(j-1))+(obj.v(k+1,j)-obj.v(k-1,j))/(obj.z(k+1)-obj.z(k-1)));
                    obj.sigmazz(k,j)=(lambda+2*mu)*((obj.w(k+1,j)-obj.w(k-1,j))/(obj.z(k+1)-obj.z(k-1)))+lambda*((obj.v(k,j+1)-obj.v(k,j-1))/(obj.r(j+1)-obj.r(j-1))+obj.v(k,j)/obj.r(j));
                    obj.sigmarr(k,j)=(lambda+2*mu)*(obj.v(k,j+1)-obj.v(k,j-1))/(obj.r(j+1)-obj.r(j-1))+lambda*(obj.w(k+1,j)-obj.w(k-1,j))/(obj.z(k+1)-obj.z(k-1));
                    obj.sigmaff(k,j)=(lambda+2*mu)*obj.v(k,j)/obj.r(j)+lambda*((obj.v(k,j+1)-obj.v(k,j-1))/(obj.r(j+1)-obj.r(j-1))+(obj.w(k+1,j)-obj.w(k-1,j))/(obj.z(k+1)-obj.z(k-1)));
                end
            end
            for k=2:obj.N-1
                j=obj.M;
                obj.sigmazr(k,j)=mu*((obj.w(k,j)-obj.w(k,j-1))/(obj.r(j)-obj.r(j-1))+(obj.v(k+1,j)-obj.v(k-1,j))/(obj.z(k+1)-obj.z(k-1)));
                obj.sigmazz(k,j)=(lambda+2*mu)*((obj.w(k+1,j)-obj.w(k-1,j))/(obj.z(k+1)-obj.z(k-1)))+lambda*((obj.v(k,j)-obj.v(k,j-1))/(obj.r(j)-obj.r(j-1))+obj.v(k,j)/obj.r(j));
                obj.sigmarr(k,j)=(lambda+2*mu)*(obj.v(k,j)-obj.v(k,j-1))/(obj.r(j)-obj.r(j-1))+lambda*(obj.w(k+1,j)-obj.w(k-1,j))/(obj.z(k+1)-obj.z(k-1));
            end
            for k=2:obj.N
                obj.sigmazr(k,1)=0;
                obj.sigmazz(k,1)=obj.sigmazz(k,2);
                obj.sigmarr(k,1)=obj.sigmarr(k,2);
            end
            t
            toc
        end
    end
end
