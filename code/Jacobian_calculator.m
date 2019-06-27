% Function that computes the Jacobian matrix 

function J = Jacobian_calculator(N, P, F0, mu, a, alpha, A)

J = zeros(2*N+1,2*N+1) ;
    
    for i = 1 : N
        for j = 1 : N
            if j == i % diagonal and (2N+1) blocks of Jacobian
                for k = 1 : N
                    if k ~= j
                        
                        r = sqrt ( (P(2*k-1)-P(2*i-1))^2 + (P(2*k)-P(2*i))^2 ) ;
                        psi = P(2*N+1) ;
                        theta = atan2( P(2*k)-P(2*i) , P(2*k-1)-P(2*i-1) ) ;
                        phi = psi - theta ;
                        del_x = P(2*k-1) - P(2*i-1) ;
                        del_y = P(2*k) - P(2*i) ;
                        
                        J(2*i-1,2*i-1) = J(2*i-1,2*i-1) + ( F0/2/mu/r^7 ) * ( -2*a^alpha*A * (alpha*del_x^2-del_y^2) * r^(4-alpha) + ...
                            a^4 * (-2*del_x-del_y)*(-2*del_x+del_y) + a^4 * ( cos(2*theta) * ((12*del_x^2-7*del_y^2)*cos(2*psi)+16*del_x*del_y*sin(2*psi)) + ...
                            ( -16*del_x*del_y*cos(2*psi) + (12*del_x^2-7*del_y^2)*sin(2*psi) ) * sin(2*theta) ) ) ;
                        
                        J(2*i-1,2*i) = J(2*i-1,2*i) + ( F0/2/mu/r^7 ) * ( del_x*del_y * ( 5*a^4-2*a^alpha*A*(1+alpha)*r^(4-alpha) ) + ...
                            a^4 * ( 19*del_x*del_y*cos(2*phi) - 8*(-del_x-del_y)*(-del_x+del_y)*sin(2*phi) ) ) ;
                        
                        J(2*i,2*i) = J(2*i,2*i) + ( F0/mu/r^7 ) * ( a^alpha*A * (del_x^2-alpha*del_y^2) * r^(4-alpha) - ...
                            a^4/2 * (-del_x-2*del_y)*(-del_x+2*del_y) + a^4/2 * ( cos(2*theta) * ((12*del_y^2-7*del_x^2)*cos(2*psi)-16*del_x*del_y*sin(2*psi)) + ...
                            ( 16*del_x*del_y*cos(2*psi) + (12*del_y^2-7*del_x^2)*sin(2*psi) ) * sin(2*theta) ) ) ;
                        
                        J(2*i,2*i-1) = J(2*i,2*i-1) + ( F0/2/mu/r^7 ) * ( del_x*del_y * ( 5*a^4-2*a^alpha*A*(1+alpha)*r^(4-alpha) ) + ...
                            a^4 * ( 19*del_x*del_y*cos(2*phi) - 8*(-del_x-del_y)*(-del_x+del_y)*sin(2*phi) ) ) ;
                        
                        J(2*i-1,2*N+1) = J(2*i-1,2*N+1) + ( F0*a^4/mu/r^5 ) * ( 2*del_y*cos(2*phi) - 3*del_x*sin(2*phi) ) ;
                        
                        J(2*i,2*N+1) = J(2*i,2*N+1) + ( F0*a^4/mu/r^5 ) * ( -2*del_x*cos(2*phi) - 3*del_y*sin(2*phi) ) ;
                        
                    end
                end
            else % off-diagonal blocks of Jacobian
                
                r = sqrt ( (P(2*j-1)-P(2*i-1))^2 + (P(2*j)-P(2*i))^2 ) ;
                psi = P(2*N+1) ;
                theta = atan2( P(2*j)-P(2*i) , P(2*j-1)-P(2*i-1) ) ;
                phi = psi - theta ;
                del_x = P(2*j-1) - P(2*i-1) ;
                del_y = P(2*j) - P(2*i) ;
                
                J(2*i-1,2*j-1) = J(2*i-1,2*j-1) + ( F0/mu/r^7 ) * ( a^alpha*A * (alpha*del_x^2-del_y^2) * r^(4-alpha) + ...
                    a^4 * (-2*del_x-del_y)*(-2*del_x+del_y) + a^4 * ( 3*(-4*del_x^2+del_y^2) * (cos(phi))^2 - 2*del_y * ...
                    ( -del_y*cos(2*phi) + 4*del_x*sin(2*phi) ) ) ) ;
                
                J(2*i-1,2*j) = J(2*i-1,2*j) + ( F0/mu/r^7 ) * ( del_x*del_y * ( 5*a^4+a^alpha*A*(1+alpha)*r^(4-alpha) ) + ...
                    a^4 * ( -del_x*del_y/2 * (15+19*cos(2*phi)) + 4*(-del_x-del_y)*(-del_x+del_y)*sin(2*phi) ) ) ;
                
                J(2*i,2*j-1) = J(2*i,2*j-1) + ( F0/mu/r^7 ) * ( del_x*del_y * ( 5*a^4+a^alpha*A*(1+alpha)*r^(4-alpha) ) + ...
                    a^4 * ( -del_x*del_y/2 * (15+19*cos(2*phi)) + 4*(-del_x-del_y)*(-del_x+del_y)*sin(2*phi) ) ) ;
                
                J(2*i,2*j) = J(2*i,2*j) + ( F0/2/mu/r^7 ) * ( -2*a^alpha*A * (del_x^2-alpha*del_y^2) * r^(4-alpha) + ...
                    a^4 * (-del_x-2*del_y)*(-del_x+2*del_y) + a^4 * ( cos(2*theta) * ( (7*del_x^2-12*del_y^2)*cos(2*psi)+16*del_x*del_y*sin(2*psi) ) + ...
                    ( -16*del_x*del_y*cos(2*psi) + (-12*del_y^2+7*del_x^2)*sin(2*psi) ) * sin(2*theta) ) ) ;
                
            end
        end
    end
