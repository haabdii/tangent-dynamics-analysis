% Function that computes the force on each particle 
function P_dot = velocity_calculator (N, w, P, C1, C2, alpha)
  P_dot = zeros(2*N+1,1) ;
    
    P_dot(2*N+1,1) = w ;
    
    for i = 1 : N
        for j = 1: N
            if j ~= i
                r = sqrt ( (P(2*j-1)-P(2*i-1))^2 + (P(2*j)-P(2*i))^2 ) ;
                phi = P(2*N+1) - atan2( P(2*j)-P(2*i) , P(2*j-1)-P(2*i-1) ) ;
                del_x = P(2*j-1) - P(2*i-1) ;
                del_y = P(2*j) - P(2*i) ;
                F_r = 3 * ( cos(phi) )^2 -1 ;
                F_t = sin(2*phi) ;
                P_dot(2*i-1) = P_dot(2*i-1) + (C1 / r^5) * ( del_x * F_r + del_y * F_t ) - C2 * del_x / r^(alpha+1) ;
                P_dot(2*i) = P_dot(2*i) + (C1 / r^5) * ( del_y * F_r - del_x * F_t ) - C2 * del_y / r^(alpha+1) ;
            end
        end
    end
end
