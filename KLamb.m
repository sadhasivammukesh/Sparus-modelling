function M = KLamb(a, b, rho)
    M = zeros(6,6);

    e = sqrt(1 - (b^2)/(a^2));

    a0 = 2*((1-e^2)/(e^3))*((1/2)*log((1+e)/(1-e)) - e);
    b0 = 1/(e^2) - ((1-e^2)/(2*e^3))*log((1+e)/(1-e));

    k1 = a0/(2-a0);
    k2 = b0/(2-b0);
    k4 = 0;
    k5 = ((e^4)*(b0-a0))/((2-e^2)*((2*e^2)-(2-e^2)*(b0-a0)));

    mdf = (4/3)*rho*pi*a*b^2;

    Idf = [
        (mdf/5)*(2*b^2), 0, 0;
        0, (mdf/5)*(a^2 + b^2), 0;
        0, 0, (mdf/5)*(a^2 + b^2)
    ];

    M(1,1) = k1*mdf;
    M(2,2) = k2*mdf;
    M(3,3) = k2*mdf;
    M(4,4) = k4*Idf(1,1);
    M(5,5) = k5*Idf(2,2);
    M(6,6) = k5*Idf(3,3);
end