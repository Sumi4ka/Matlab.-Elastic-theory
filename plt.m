classdef plt
    properties
        N          %Count of nodes in the axial direction
        M          %Count of nodes in the radial direction
        h          %Height
        R1         %Inner radius
        R2         %Outter radius
        r1         %Radius distribution
        z1         %Height distribution
        v          %Axial displacement distribution
        w          %Axial displacement distribution
        sigmazr    %Matrix of sigmzr distribution
        sigmazz    %Matrix of sigmzz distribution
        sigmarr    %Matrix of sigmrr distribution
        sigmaff    %Matrix of sigmff distribution
    end
    methods
        function obj=plt(solver)
            obj.N=solver.N;
            obj.M=solver.M;
            obj.h=solver.z(solver.N);
            obj.R1=solver.r(1);
            obj.R2=solver.r(solver.M);
            obj.v=solver.v;
            obj.w=solver.w;
            obj.sigmazr=solver.sigmazr;
            obj.sigmazz=solver.sigmazz;
            obj.sigmarr=solver.sigmarr;
            obj.sigmaff=solver.sigmaff;
            obj.r1=zeros(obj.N,obj.M);
            for j=1:obj.N
                obj.r1(j,:)=solver.r;
            end
            obj.z1=zeros(obj.N,obj.M);
            for k=1:obj.M
                obj.z1(:,k)=(solver.z)';
            end
        end
        %Graph of mesh
        function mesh_plot(obj)
            hold on;
            for j=1:obj.M
                plot(obj.r1(1,j)*ones(1,obj.N),0:obj.h/(obj.N-1):obj.h,'r');
            end
            for k=1:obj.N
                plot(obj.R1:(obj.R2-obj.R1)/(obj.M-1):obj.R2,obj.z1(k,1)*ones(1,obj.M),'r');
            end
            xlabel("r");
            ylabel("z");
        end
        %Graph of the difference in the grid
        function dmesh_plot(obj)
            hold on;
            subplot(1,2,1);
            plot(1:obj.M-1,obj.r1(1,2:obj.M)-obj.r1(1,1:obj.M-1));
            subplot(1,2,2);
            plot(1:obj.N-1,obj.z1(2:obj.N,1)-obj.z1(1:obj.N-1,1));
        end
        %Graph of radius sestion
        function r_plot(obj,N)
            hold on;
            s=size(N);
            for i=1:s(2)
                plot(obj.r1(N(i),:)+obj.v(N(i),:),obj.z1(N(i),:)+obj.w(N(i),:));
            end
            xlabel("r");
            ylabel("z");
        end
        %Graph of heigh sestion
        function z_plot(obj,M)
            hold on;
            s=size(M);
            for i=1:s(2)
                plot(obj.r1(:,M(i))+obj.v(:,M(i)),obj.z1(:,M(i))+obj.w(:,M(i)));
            end
            xlabel("r");
            ylabel("z");
        end
        %Graph of left and right Heigh sestions
        function z_radius_plot(obj)
            hold on
            mx1=max(obj.z1(:,1)+obj.w(:,1));
            mx2=max(obj.z1(:,obj.M)+obj.w(:,obj.M));
            mn=min(mx1,mx2);
            xx=(0:mn/1000:mn);
            y1=spline(obj.z1(:,1)+obj.w(:,1),obj.r1(:,1)+obj.v(:,1),xx);
            y2=spline(obj.z1(:,obj.M)+obj.w(:,obj.M),obj.r1(:,obj.M)+obj.v(:,obj.M),xx);
            plot(xx,y2-y1);
            xlabel("z");
            ylabel("h");
        end
        %v- and w-displacement Graph of Heigh M
        function z_vw_plot(obj,M)
            s=size(M);
            for i=1:s(2)
                subplot(s(2),2,s(2)*2-1);
                plot(obj.z1(:,M(i)),obj.v(:,M(i))/max(abs(obj.v(:,M(i)))));
                xlabel("z");
                ylabel("v/v0");
                subplot(s(2),2,s(2)*2);
                plot(obj.z1(:,M(i)),obj.w(:,M(i))/max(abs(obj.v(:,M(i)))));
                xlabel("z");
                ylabel("w/v0");
            end
        end
        %v- and w-displacement Graph of Radius N
        function r_vw_plot(obj,N)
            s=size(N);
            for i=1:s(2)
                subplot(s(2),2,s(2)*2-1);
                plot(obj.r1(N(i),:),obj.v(N(i),:)/max(abs(obj.v(N(i),:))));
                xlabel("z");
                ylabel("v/v_0");
                subplot(s(2),2,s(2)*2);
                plot(obj.r1(N(i),:),obj.w(N(i),:)/max(abs(obj.v(N(i),:))));
                xlabel("z");
                ylabel("w/v_0");
            end
        end
        % Gray Heatmap of sigma
        function sigma_gray_plot(obj,N1,t)
            s=size(obj.sigmazr);
            switch t
                case 1
                    st="ZZ";
                    c=obj.sigmazz;
                case 2
                    st="RR";
                    c=obj.sigmarr;
                case 3
                    st="FF";
                    c=obj.sigmaff;
                case 4
                    st="ZR";
                    c=obj.sigmazr;
            end
            cmap1=zeros(N1,3);
            a=10000;
            for i=1:N1
                cmap1(i,1)=log((i+a-1)/a)/log((N1+a-1)/a);
                cmap1(i,2)=log((i+a-1)/a)/log((N1+a-1)/a);
                cmap1(i,3)=log((i+a-1)/a)/log((N1+a-1)/a);
            end
            subplot(1,2,1);
            c1=flip((c(2:s(1)-1,2:s(2)-1)),1);
            c1=(c1>0).*c1;
            h1=heatmap(string(obj.r1(1,2:s(2)-1)),flip(string((obj.z1(2:s(1)-1,1))')),c1);
            title(strcat("sigma",st,">0"));
            xlabel('r');
            ylabel('z');
            colormap(subplot(1,2,1),flip(cmap1));
            a=10000;
            for i=1:N1
                cmap1(i,1)=log((i+a-1)/a)/log((N1+a-1)/a);
                cmap1(i,2)=log((i+a-1)/a)/log((N1+a-1)/a);
                cmap1(i,3)=log((i+a-1)/a)/log((N1+a-1)/a);
            end
            subplot(1,2,2);
            c2=flip((c(2:s(1)-1,2:s(2)-1)),1);
            f=c2<0;
            c2=f.*c2;
            h2=heatmap(string(obj.r1(1,2:s(2)-1)),flip(string((obj.z1(2:s(1)-1,1))')),c2);
            title(strcat("sigma",st,"<0"));
            xlabel('r');
            ylabel('z');
            colormap(subplot(1,2,2),cmap1)
            h1.ColorScaling = 'log';
            h2.ColorScaling = 'log';
        end
        % Turbo Heatmap of sigma
        function sigma_turbo_plot(obj,t)
            s=size(obj.sigmazr);
            switch t
                case 1
                    heatmap(string(obj.r1(1,2:s(2)-1)),flip(string((obj.z1(2:s(1)-1,1))')),flip((obj.sigmazz(2:s(1)-1,2:s(2)-1)),1));
                    title("zz");
                case 2
                    heatmap(string(obj.r1(1,2:s(2)-1)),flip(string((obj.z1(2:s(1)-1,1))')),flip(obj.sigmarr(2:s(1)-1,2:s(2)-1),1));
                    title("rr");
                case 3
                    heatmap(string(obj.r1(1,2:s(2)-1)),flip(string((obj.z1(2:s(1)-1,1))')),flip(obj.sigmaff(2:s(1)-1,2:s(2)-1),1));
                    title("ff");
                case 4
                    heatmap(string(obj.r1(1,2:s(2)-1)),flip(string((obj.z1(2:s(1)-1,1))')),flip(obj.sigmazr(2:s(1)-1,2:s(2)-1),1));
                    title("zr");
            end
            xlabel('r');
            ylabel('z');
            colormap turbo
        end
        %Heatmap of plastic zones
        function plastic(obj,sigmaS)
            s=size(obj.sigmazr);
            c=zeros(obj.N,obj.M);
            mx=0;
            for j=2:obj.N-1
                for k=2:obj.M-1
                    f=1/sqrt(2)*sqrt((obj.sigmarr(j,k)-obj.sigmaff(j,k))^2+(obj.sigmarr(j,k)-obj.sigmazz(j,k))^2+(obj.sigmazz(j,k)-obj.sigmaff(j,k))^2+6*obj.sigmazr(j,k)^2);
                    if f>mx
                        mx=f;
                    end
                    if f>sigmaS
                        c(j,k)=1;
                    else
                        c(j,k)=0;
                    end
                end
            end
            heatmap(string(obj.r1(1,2:s(2)-1)),flip(string((obj.z1(2:s(1)-1,1))')),flip(c(2:s(1)-1,2:s(2)-1),1));
            title(strcat(string(round(mx,-6)/1000000)," MPa"));
            colormap gray
        end
    end
end