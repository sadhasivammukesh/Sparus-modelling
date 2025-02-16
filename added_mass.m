function Ma = calculate_hull_added_mass(L0, L1, L2, L3, L, R, rho)
    syms x;

    x0 = -L0;
    x1 = x0 + L1;
    x2 = x1 + L3;
    x3 = L - L0;

    % Integration intervals
    interval1 = [x0, x1];
    interval2 = [x1, x2];
    interval3 = [x2, x3];

    q11 = int(1 - (((x - x1)^2) / L1^2), interval1);
    q12 = int(x / x, interval2);
    q13 = int(1 - (((x - x2)^2) / L2^2), interval3);

    m22 = rho * pi * R^2 * (q11 + q12 + q13);
    m33 = m22;
    m44 = 0;

    q21 = int(x*(1 - (((x - x1)^2) / L1^2)), interval1);
    q22 = int(x, interval2);
    q23 = int(x*(1 - (((x - x2)^2) / L2^2)), interval3);

    m26 = rho * pi * R^2 * (q21 + q22 + q23);
    m35 = -m26;
    m53 = m35;
    m62 = m26;

    q31 = int((x^2)*(1 - (((x - x1)^2) / L1^2)), interval1);
    q32 = int(x^2, interval2);
    q33 = int((x^2)*(1 - (((x - x2)^2) / L2^2)), interval3);

    m55 = rho * pi * R^2 * (q31 + q32 + q33);
    m66 = m55;

    % Added mass matrix
    Ma = [
        0, 0, 0, 0, 0, 0;
        0, double(m22), 0, 0, 0, double(m26);
        0, 0, double(m33), 0, double(m35), 0;
        0, 0, 0, double(m44), 0, 0;
        0, 0, double(m53), 0, double(m55), 0;
        0, double(m62), 0, 0, 0, double(m66);
        ];
end

function Ma = calculate_thruster_added_mass(L, R, rho)
    syms x;

    % Integration intervals
    interval = [-L/2, L/2];

    m22 = 0; % projected area on the main body
    m33 = rho * pi * R^2 * int(x / x, interval);
    m44 = 0;

    m26 = rho * pi * R^2 * int(x, interval);
    m35 = -m26;
    m53 = m35;
    m62 = m26;

    m55 = rho * pi * R^2 * int(x^2, interval);
    m66 = m55;

    % Added mass matrix
    Ma = [
        0, 0, 0, 0, 0, 0;
        0, double(m22), 0, 0, 0, double(m26);
        0, 0, double(m33), 0, double(m35), 0;
        0, 0, 0, double(m44), 0, 0;
        0, 0, double(m53), 0, double(m55), 0;
        0, double(m62), 0, 0, 0, double(m66);
        ];
end

function Ma = calculate_antenna_added_mass(W, H, L, rho)
    syms x;

    CA1 = 1.7; % a/b = 0.37
    CA2 = 1.36; % a/b = 2.7
    b1 = 0.15; % a/b = 2.7

    % Integration intervals
    interval = [-H/2, H/2];

    m11 = int((x/ x) * rho * CA1 * (pi * (W/2)^2), interval);
    m22 = int((x/ x) * rho * CA2 * (pi * (L/2)^2), interval);
    m33 = 0;
    m44 = int((x^2) * rho * CA2 * (pi * (L/2)^2), interval);
    m55 = int((x^2) * rho * CA1 * (pi * (W/2)^2), interval);
    m66 = int((x/ x) * b1 * rho * pi * (L/2)^4, interval);

    m26 = int(x * rho * CA2 * (pi * (L/2)^2), interval);
    m35 = 0;
    m53 = m35;
    m62 = m26;

    % Added mass matrix
    Ma = [
        double(m11), 0, 0, 0, 0, 0;
        0, double(m22), 0, 0, 0, double(m26);
        0, 0, double(m33), 0, double(m35), 0;
        0, 0, 0, double(m44), 0, 0;
        0, 0, double(m53), 0, double(m55), 0;
        0, double(m62), 0, 0, 0, double(m66);
        ];
end

Ma_hull = calculate_hull_added_mass(0.7465, 0.26, 0.0974, 1.2426, 1.6, 0.115, 1000);
Ma_thruster = calculate_thruster_added_mass(0.2365, 0.0435, 1000);
Ma_antenna = calculate_antenna_added_mass(0.0244, 0.254, 0.066, 1000);

r_hull = [0; 0; -0.02];
r_thuster1 = [-0.5; -0.1585; 0.02];
r_thuster2 = [-0.5; 0.1585; 0.02];
r_antenna = [-0.389; 0; -0.222];

Ma_hull_cb = transformMatrix(r_hull, Ma_hull);
Ma_thruster1_cb = transformMatrix(r_thuster1, Ma_thruster);
Ma_thruster2_cb = transformMatrix(r_thuster2, Ma_thruster);
Ma_antenna_cb = transformMatrix(r_antenna, Ma_antenna);

Ma_total_cb = Ma_hull_cb + Ma_thruster1_cb + Ma_thruster2_cb + Ma_antenna_cb;

Ma_total_cg = transformMatrix([0;0;0.02], Ma_total_cb);