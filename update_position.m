% Function that updates the position vector 
function P = update_position (N, TimeStep, P, P_dot)
    for i = 1 : N
        P(2*i-1) = P(2*i-1) + P_dot(2*i-1) * TimeStep ;
        P(2*i) = P(2*i) + P_dot(2*i) * TimeStep ;
    end
    P(2*N+1) = P(2*N+1) + P_dot(2*N+1) * TimeStep ;
end