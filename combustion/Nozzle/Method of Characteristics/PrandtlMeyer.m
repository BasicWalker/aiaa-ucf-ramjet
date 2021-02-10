function nu = PrandtlMeyer (g, M)
    nu = ( sqrt((g+1)/(g-1)) ) * atand( sqrt( ((g-1)/(g+1))* (M^2-1) ) ) - ...
            atand( sqrt( M^2 - 1) ); 
end

    