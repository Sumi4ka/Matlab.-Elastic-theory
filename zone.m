classdef zone
    properties
        N          %Count of nodes in the axial direction
        M          %Count of nodes in the radial direction
        z_left_v   %Condition for Radial displacement in the left side
        z_left_w   %Condition for Axial displacement in the left side
        z_right_v  %Condition for Radial displacement in the right side
        z_right_w  %Condition for Axial displacement in the right side
        z_bot_v    %Condition for Radial displacement in the bot side
        z_bot_w    %Condition for Axial displacement in the bot side
        z_top_v    %Condition for Radial displacement in the top side
        z_top_w    %Condition for Axial displacement in the top side
    end
    methods
        %Constructor 
        function obj=zone(N,M)
            obj.N=N;
            obj.M=M;
            obj.z_left_v=zeros(1,M);
            obj.z_left_w=zeros(1,M);
            obj.z_right_v=zeros(1,N);
            obj.z_right_w=zeros(1,N);
            obj.z_bot_v=zeros(1,M);
            obj.z_bot_w=zeros(1,M);
            obj.z_top_v=zeros(1,N);
            obj.z_top_w=zeros(1,N);
        end
        % Free surface conditions for:
        % Radial displacement on the left side
        function obj=z_left_v_free(obj,n1,n2)
            obj.z_left_v(n1:n2)=-1;
        end
        % Axial  displacement on the left side
        function obj=z_left_w_free(obj,n1,n2)
            obj.z_left_w(n1:n2)=-1;
        end
        % Radial displacement on the right side
        function obj=z_right_v_free(obj,n1,n2)
            obj.z_right_v(n1:n2)=-1;
        end
        % Axial  displacement on the right side
        function obj=z_right_w_free(obj,n1,n2)
            obj.z_right_w(n1:n2)=-1;
        end
        % Radial displacement on the bot side
        function obj=z_bot_v_free(obj,n1,n2)
            obj.z_bot_v(n1:n2)=-1;
        end
        % Axial  displacement on the bot side
        function obj=z_bot_w_free(obj,n1,n2)
            obj.z_bot_w(n1:n2)=-1;
        end
        % Radial displacement on the top side
        function obj=z_top_v_free(obj,n1,n2)
            obj.z_top_v(n1:n2)=-1;
        end
        % Axial  displacement on the top side
        function obj=z_top_w_free(obj,n1,n2)
            obj.z_top_w(n1:n2)=-1;
        end
        
        % Static press condition for
        % Radial displacement on the left side
        function obj=z_left_v_press(obj,n1,n2,depth)
            obj.z_left_v(n1:n2)=depth;
        end
        % Radial displacement on the right side
        function obj=z_right_v_press(obj,n1,n2,depth)
            obj.z_right_v(n1:n2)=depth;
        end
        % Axial  displacement on the top side
        function obj=z_top_w_press(obj,n1,n2,depth)
            obj.z_top_w(n1:n2)=depth;
        end
        % Radial displacement on the bot side
        function obj=z_bot_v_press(obj,n1,n2,depth)
            obj.z_bot_v(n1:n2)=depth;
        end
    end
end