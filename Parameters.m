function Para=Parameters()
    global Para
    
    added_mass;
    
    %% Initial Speed and position in Earth-fixed frame
    
    Para.ICPos = [0 0 2 0 0 0];
    Para.ICSpeed = [0 0 0 0 0 0] ;
    
    %% General parameters
    Para.rho_water = 1000 ;                     % Masse volumique de l'eau (kg/m^3)
    Para.R = 0.115 ;                             % Sparus Radius (m)
    Para.L = 1.6;  	                            % Sparus length (m)
    Para.m = 52 ; 	                            % Sparus mass (kg)
    Para.mb = 52.1;                           	% Sparus buoyancy mass  (kg) 
    Para.g = 9.81 ;                             % Earth Gravity (m*s^(-2))
    Para.P = Para.m * Para.g;	                % Sparus weight (N)
    Para.B = Para.mb * Para.g;	                % Buoyancy (N)
    
    %% Center of gravity and Buoyancy position in body-fixed frame
    
    Para.xg = 0 ;    %x-positon of Center of gravity
    Para.yg = 0 ;    %y-positon of Center of gravity
    Para.zg = 0 ;    %z-positon of Center of gravity
    
    Para.rg = [Para.xg Para.yg Para.zg]' ;
    
    
    Para.xb = 0      ;    % x-positon of Center of Buoyancy
    Para.yb = 0      ;    % y-positon of Center of Buoyancy
    Para.zb = -0.02  ;    % z-positon of Center of Buoyancy
    
    Para.rb = [Para.xb Para.yb Para.zb]' ;
    
    %% Body positions
    
    % Main Body S0;
    Para.S0.r=[0,0,0]'; % Position (m)
    % First Body S1 (Thruster 1);
    Para.S1.r=[-0.5,-0.159,0]'; % Position (m)
    % Second Body S2 (Thruster 2);
    Para.S2.r=[-0.5,0.159,0]'; % Position (m)
    % Second Body S3 (Antenna);
    Para.S3.r=[0.39,0,-0.24]'; % Position (m)
    % Thruster 3 (middle);
    Para.S4.r=[0,0,0.08]'; % Position (m)
   
    Para.DVL.r = [-0.4145; 0; 0.11];
    
    %% Body Mass matrices
    
    % Global real mass matrix
    Para.Mb = [52, 0, 0, 0, -0.1, 0;
               0, 52, 0, 0.1, 0, -1.3;
               0, 0, 52, 0, 1.3, 0;
               0, 0.1, 0, 0.5, 0, 0;
               -0.1, 0, 1.3, 0, 9.4, 0;
               0, -1.3, 0, 0, 0, 9.5];

    %% Body added Mass matrices
    
    % Main Body S0;
    Para.S0.Ma = Ma_hull;
    
    % First Body S1 (Thruster 1);
    Para.S1.Ma = transformMatrix(-Para.rb, Ma_thruster1_cb);
    
    % Second Body S2 (Thruster 2);
    Para.S2.Ma = transformMatrix(-Para.rb, Ma_thruster2_cb);
    
    % Third Body S3 (Antenna);
    Para.S3.Ma = transformMatrix(-Para.rb, Ma_antenna_cb);
    
    %% Generalized mass matrix
    
    % Para.S0.Mg = Para.S0.Mb + Para.S0.Ma ; 
    % Para.S1.Mg = Para.S1.Mb + Para.S1.Ma ;
    % Para.S2.Mg = Para.S2.Mb + Para.S2.Ma ;
    
    Para.Ma = Ma_total_cg;
    Para.Mg = Para.Mb + Para.Ma;
    
    
    %% Generalized coriolis matrix
    
    % Computed in RovModel.m
    
    %% Friction matrices
    
    % Main Body S0;
    Para.S0.Kq = [2.08 0 0 0 0 0;
                  0 55.2 0 0 0 0;
                  0 0 55.2 0 0 0;
                  0 0 0 0 0 0;
                  0 0 0 0 7.07 0;
                  0 0 0 0 0 7.07];
    
    % First Body S1 (Thruster 1);
    Para.S1.Kq = [2.66 0 0 0 0 0;
                  0 0 0 0 0 0;
                  0 0 3.09 0 0 0;
                  0 0 0 0 0 0;
                  0 0 0 0 0 0;
                  0 0 0 0 0 0];
    
    % Second Body S2 (Thruster 2);
    Para.S2.Kq = [2.66 0 0 0 0 0;
                  0 0 0 0 0 0;
                  0 0 3.09 0 0 0;
                  0 0 0 0 0 0;
                  0 0 0 0 0 0;
                  0 0 0 0 0 0];
    
    % Third Body S3 (Antenna);
    Para.S3.Kq = [4.03 0 0 0 0 0;
                  0 20.96 0 0 0 0;
                  0 0 0 0 0 0;
                  0 0 0 0.01 0 0;
                  0 0 0 0 0 0;
                  0 0 0 0 0 0];
    
    %% Thruster modelling
    
    %Thruster positions in body-fixed frame
    
    Para.d1x = Para.S4.r(1); 
    Para.d1y = Para.S4.r(2);
    Para.d1z = Para.S4.r(3);
    Para.d2x = Para.S2.r(1);
    Para.d2y = Para.S2.r(2);
    Para.d2z = Para.S2.r(3);
    Para.d3x = Para.S1.r(1);
    Para.d3y = Para.S1.r(2);
    Para.d3z = Para.S1.r(3);
    
    
    Para.rt1 = [Para.d1x, Para.d1y, Para.d1z]' ;
    Para.rt2 = [Para.d2x, Para.d2y, Para.d2z]' ;
    Para.rt3 = [Para.d3x, Para.d3y, Para.d3z]' ;
    
    
    Para.rt = [Para.rt1 Para.rt2 Para.rt3] ;
    
    %Thruster gains
    
    % Forward
    Para.kt1f = 55    ;
    Para.kt2f = 71.5    ;
    Para.kt3f = 71.5    ;

    % Rearward
    Para.kt1r = 28.5    ;
    Para.kt2r = 30    ;
    Para.kt3r = 30    ;

    Para.kt1 = 55    ;
    Para.kt2 = 71.5    ;
    Para.kt3 = 71.5    ;

    %%%%%%%%%%%%%%%

    Para.Rkt1 = 28.5;
    Para.Rkt2 = 30;
    Para.Rkt3 = 30;
    
    Para.Fkt1 = 55;
    Para.Fkt2 = 71.5;
    Para.Fkt3 = 71.5;
    
    
    
    %Thruster gain vectors
    
    Para.Kt=[Para.kt1f;Para.kt2f;Para.kt3f];
    
    %%%%%%%%%%%%%%%%%

    Para.RKt=[Para.Rkt1;Para.Rkt2;Para.Rkt3];

    Para.FKt=[Para.Fkt1;Para.Fkt2;Para.Fkt3];

    %Thruster time constants
    
    Para.Tau1 = 0.4 ;
    Para.Tau2 = 0.8 ;
    Para.Tau3 = 0.8 ;
    
    
    %Thruster time constant vectors
    
    Para.Tau = [Para.Tau1;Para.Tau2;Para.Tau3] ;
    
    % Mapping of thruster
    
    % Para.Eb_F = Para.Kt'.*[0,1,1;0,0,0;1,0,0];
    % 
    % Para.Eb_M = [S_(Para.S4.r) * [0;0;1], S_(Para.S2.r) * [1;0;0], S_(Para.S1.r) * [1;0;0]];
    % 
    % Para.Eb = [ Para.Eb_F ; Para.Eb_M ] ;

    Para.Eb_F = [0 1 1; 0 0 0; 1 0 0];

    Para.Eb_M = [0 0 0; 0 0 0; 0 -Para.d2y -Para.d3y];

    Para.Eb = [ Para.Eb_F ; Para.Eb_M ] ;
    
    % Inverse Mapping of thruster
    Para.Ebinv = pinv(Para.Eb) ;

end





 
           

