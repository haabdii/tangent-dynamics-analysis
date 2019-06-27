% Function that updates separation vector
function E = separation_update(N, E, E_dot, TimeStep)
    for i = 1 : N
        E(2*i-1) = E(2*i-1) + E_dot(2*i-1) * TimeStep ;
        E(2*i) = E(2*i) + E_dot(2*i) * TimeStep ;
    end
    E(2*N+1) = E(2*N+1) + E_dot(2*N+1) * TimeStep ;
end
