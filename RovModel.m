function [AccG] = RovModel(Thrust,PosE,VitB)

    global Para
    
    %% Attitudes in earth frame
    % z=PosE(3,1);

    phi     = PosE(4,1)	;
    theta   = PosE(5,1)	;
    
    %% Gravity Force
    
    Fg = 1* [-Para.P * sin(theta) ;
            Para.P * cos(theta)*sin(phi) ;
            Para.P * cos(theta)*cos(phi) ;
            0 ;
            0 ;
            0 ];
        
    % Expressed in b and computed in G
        
    %% Force d'Archimède
    
    Fa_F = [Para.B * sin(theta) ;
            -Para.B * cos(theta)*sin(phi) ;
            -Para.B * cos(theta)*cos(phi) ;
            ];
    %  Expressed in b
    
    
    Fa_M = S_(Para.rb-Para.rg) * Fa_F ; % Computed in G
    
    Fa = [ Fa_F ; Fa_M ] ;
    %  Expressed in b and computed in G
    
    %% Force de Coriolis
    
    u = VitB(1,1)   ;
    v = VitB(2,1)   ;
    w = VitB(3,1)   ;
    V_ = [u;v;w];
    p = VitB(4,1)   ;   %Body fixed velocity roll (rad*s^(-1))
    q = VitB(5,1)   ;   %Body fixed velocity pitch (rad*s^(-1))
    r = VitB(6,1)   ;   %Body fixed velocity yaw (rad*s^(-1))
    W_ = [p;q;r]     ;  %General vector
    
    
    
    % General coriolis matrix :
    C_g = [Para.m * S_(W_), zeros(3);
            zeros(3), -S_(Para.Mb(4:6,4:6)*W_)];

    C_a = [zeros(3), -S_(Para.Ma(1:3,1:3)*V_ + Para.Ma(1:3,4:6)*W_);
           -S_(Para.Ma(1:3,1:3)*V_ + Para.Ma(1:3,4:6)*W_), -S_(Para.Ma(4:6,1:3)*V_ + Para.Ma(4:6,4:6)*W_)];
    
    % C_all = C_a + C_g;
    C_all = C_a;

    %coriolis Force :
    Fc = C_all * VitB;

    function H = Hmatrix(r)
        H = [eye(3), (S_(r))'; zeros(3, 3), eye(3)];
    end
    
    %% Friction forces
    Vit_0=VitB;
    r_b_dvl = Para.DVL.r - Para.rb;
    Vit_0b = Hmatrix(r_b_dvl)*Vit_0;
    % Vit_0b = Hmatrix(Para.S0.r)'*Vit_0;
    Ff_0 =  Para.S0.Kq * abs(Vit_0b).*Vit_0b ;
    Ff_0g = Hmatrix(Para.rb)' * Ff_0; % Friction force on the main hull at CG of the body
    
    Vit_1=VitB;
    r_b_dvl = Para.DVL.r - Para.S1.r;
    Vit_1b = Hmatrix(r_b_dvl)*Vit_1;
    % Vit_1b = Hmatrix(Para.S1.r)'*Vit_1;
    Ff_1 =  Para.S1.Kq * abs(Vit_1b).*Vit_1b ; % Friction force on the buoyancy center of the Thruster 1
    Ff_1g = Hmatrix(Para.S1.r)' * Ff_1; % Friction force on Thruster 1 at CG of the body
    
    Vit_2=VitB;
    r_b_dvl = Para.DVL.r - Para.S2.r;
    Vit_2b = Hmatrix(r_b_dvl)*Vit_2;
    % Vit_2b = Hmatrix(Para.S2.r)'*Vit_2;
    Ff_2 =  Para.S2.Kq * abs(Vit_2b).*Vit_2b ; % Friction force on the buoyancy center of the Thruster 2
    Ff_2g = Hmatrix(Para.S2.r)' * Ff_2; % Friction force on Thruster 2 at CG of the body
    
    Vit_3=VitB;
    r_b_dvl = Para.DVL.r - Para.S3.r;
    Vit_3b = Hmatrix(r_b_dvl)*Vit_3;
    % Vit_3b = Hmatrix(Para.S3.r)'*Vit_3;
    Ff_3 =  Para.S3.Kq * abs(Vit_3b).*Vit_3b ; % Friction force on the buoyancy center of the antenna
    Ff_3g = Hmatrix(Para.S3.r)' * Ff_3; % Friction force on antenna at CG of the body

    Ff = Ff_0g + Ff_1g + Ff_2g + Ff_3g;
    
    %% Propulsions Forces
    Fpt = Thrust'.*Para.Eb;
    Fp = sum([Hmatrix(Para.S4.r)' * Fpt(:,1), Hmatrix(Para.S2.r)' * Fpt(:,2), Hmatrix(Para.S1.r)' * Fpt(:,3)], 2);
    % Fp = Para.Eb * Thrust ;

    % force for applying linear acceleration and testing the impact of the different coefficients in the global mass matrix
    Fext = [0;
            0;
            0;
            0;
            0;
            0];

    %% Accelearion computation :
    AccG = Para.Mg\ (-Ff + Fa + Fg + Fp - Fc) ; % Mg\ = Mg^-1 computed at the gravity center of the Sparus
    % AccG = Para.Mg\ (Fext);
end
