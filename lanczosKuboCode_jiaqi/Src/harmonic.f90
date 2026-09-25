subroutine harmonicApprox(dx,dy,A,B,C, Hjj)
    use constants

    double precision, intent(in) :: dx, dy, A, B, C
    double precision, intent(out) :: Hjj
    double precision :: x1, y1, x2, y2, G1, m1, m2, a1, a2, a3, a4, alpha, beta, gammma
    double precision :: delta, Ax, Ay, acc, f1, D, c0, c1, phi
    
    if (abs(B-C) < 0.0000001) then
        print*, "we have a singularity"
        !D = (A-B)/(10^(-16))
    else
        D = (A-B)/(B-C)
    end if

    x1 = 0.0_dp
    y1 = 1.0_dp/sqrt(3.0_dp)
    x2 = 0.0_dp
    y2 = 0.0_dp
    G1 = 4.0_dp*pi/(sqrt(3.0_dp))
    m1 = cos(sqrt(3.0_dp)*G1*x1/2.0_dp)
    m2 = cos(sqrt(3.0_dp)*G1*x2/2.0_dp)
    a1 = cos(G1*y1) - cos(G1*y2)
    a3 = 2.0_dp*cos(G1*y1/2.0_dp) * m1 - 2.0_dp*cos(G1*y2/2.0_dp)* m2
    alpha = a1 + a3
    a2 = sin(G1*y1) - sin(G1*y2)
    a4 = 2.0_dp*sin(G1*y1/2.0_dp) * m1 - 2.0_dp*sin(G1*y2/2.0_dp)* m2
    beta = a2 - a4
    x1 = 0.0_dp
    y1 = 0.0_dp
    x2 = 0.0_dp
    y2 = 2.0_dp/sqrt(3.0_dp)
    G1 = 4.0_dp*pi/(sqrt(3.0_dp))
    m1 = cos(sqrt(3.0_dp)*G1*x1/2.0_dp)
    m2 = cos(sqrt(3.0_dp)*G1*x2/2.0_dp)
    a1 = cos(G1*y1) - cos(G1*y2)
    a3 = 2.0_dp*cos(G1*y1/2.0_dp) * m1 - 2.0_dp*cos(G1*y2/2.0_dp)* m2
    gammma = a1 + a3
    a2 = sin(G1*y1) - sin(G1*y2)
    a4 = 2.0_dp*sin(G1*y1/2.0_dp) * m1 - 2.0_dp*sin(G1*y2/2.0_dp)* m2
    delta = a2 - a4
    phi = atan((1.0_dp/(delta/beta*D-1.0_dp)*((delta*alpha-beta*gammma)/(beta*delta)))-(gammma)/(delta))
    if (abs(B-C) < 0.0000001_dp) then
        !c1 = (10^(-16))/(2.0_dp*(gammma*cos(phi)+delta*sin(phi)))
        print*, "we have a singularity2"
    else
        c1= (B-C)/(2.0*(gammma*cos(phi)+delta*sin(phi)))
    end if
    Ax = 0.0_dp
    Ay = 1_dp/sqrt(3.0_dp)
    G1 = 4.0_dp*pi/(sqrt(3.0_dp))
    c0 = A - 2.0_dp * c1 * cos(phi - G1 * Ay) - 4.0_dp * c1 * cos(G1 * Ay / 2.0_dp + phi) * cos(sqrt(3.0) * G1 * Ax / 2.0)
    acc = 1.0_dp
    G1 = 4.0_dp*pi/(sqrt(3.0_dp)*acc)
    f1 = 2.0_dp*c1*cos(phi-G1*dy) + 4.0_dp*c1*cos(G1*dy/2.0_dp + phi)*cos(sqrt(3.0_dp)*G1*dx/2.0_dp)
    Hjj = c0 + f1

    return

end subroutine
