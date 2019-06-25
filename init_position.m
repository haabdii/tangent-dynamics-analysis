% Function that initializes the position of the particles 

function P = init_position (N, a, IC)
P = zeros(2*N+1,1) ; % including the x and y components 
% of the particle positions and the field angle

for i = 1 : N
    P(2*i-1,1) = (i-1)*a ;
    P(2*i) = 0 ;
end

P(2*N+1,1) = IC ;

end