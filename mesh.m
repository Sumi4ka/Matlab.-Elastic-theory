classdef mesh
    properties
        r          %Radius distribution
        z          %Height distribution
        N          %Count of nodes in the axial direction
        M          %Count of nodes in the radial direction
        H          %Height
        R1         %Inner radius
        R2         %Outter radius
    end
    methods
        %Constructor
        function obj=mesh(N,M,H,R1,R2)
            obj.r=linspace(R1,R2,M);
            obj.z=linspace(0,H,N);
            obj.N=N;
            obj.M=M;
            obj.H=H;
            obj.R1=R1;
            obj.R2=R2;
        end
        %Right concentrating in radial direction
        function obj=r_right(obj,n1,n2,si)
            dro=(log(si)-log(obj.r(n2)-obj.r(n1)+si))/(n2-n1);
            ro=zeros(1,n2-n1-1);
            for j=1:n2-n1
                ro(j)=log(obj.r(n2)-obj.r(n1)+si)+(j-1)*dro;
            end
            obj.r(n1+1:n2-1)=obj.r(n2)+si-exp(ro(2:n2-n1));
        end
        %Left concentrating in radial direction
        function obj=r_left(obj,n1,n2,si)
            dro=(log(obj.r(n2)-obj.r(n1)+si)-log(si))/(n2-n1);
            ro=zeros(1,n2-n1-1);
            for j=1:n2-n1
                ro(j)=log(si)+(j-1)*dro;
            end
            obj.r(n1+1:n2-1)=obj.r(n1)-si+exp(ro(2:n2-n1));
        end
        %Right concentrating in axial direction
        function obj=z_right(obj,n1,n2,si)
            dro=(log(si)-log(obj.z(n2)-obj.z(n1)+si))/(n2-n1);
            ro=zeros(1,n2-n1-1);
            for j=1:n2-n1
                ro(j)=log(obj.z(n2)-obj.z(n1)+si)+(j-1)*dro;
            end
            obj.z(n1+1:n2-1)=obj.z(n2)+si-exp(ro(2:n2-n1));
        end
        %Left concentrating in axial direction
        function obj=z_left(obj,n1,n2,si)
            dro=(log(obj.z(n2)-obj.z(n1)+si)-log(si))/(n2-n1);
            ro=zeros(1,n2-n1-1);
            for j=1:n2-n1
                ro(j)=log(si)+(j-1)*dro;
            end
            obj.z(n1+1:n2-1)=obj.z(n1)-si+exp(ro(2:n2-n1));
        end
    end
end