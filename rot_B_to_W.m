%% Rotate from body frame to world frame (FRD to NED matrix from RPY in radians)
function R_B_to_W = rot_B_to_W(roll, pitch, yaw)
    sphi = sin(roll);
    cphi = cos(roll);
    stheta = sin(pitch);
    ctheta = cos(pitch);
    spsi = sin(yaw);
    cpsi = cos(yaw);

    R_B_to_W = zeros(3,3);
    R_B_to_W(1, 1) = ctheta * cpsi;
    R_B_to_W(1, 2) = sphi * stheta * cpsi - cphi * spsi;
    R_B_to_W(1, 3) = cphi * stheta * cpsi + sphi * spsi;
    
    R_B_to_W(2, 1) = ctheta * spsi;
    R_B_to_W(2, 2) = sphi * stheta * spsi + cphi * cpsi;
    R_B_to_W(2, 3) = cphi * stheta * spsi - sphi * cpsi;
    
    R_B_to_W(3, 1) = -stheta;
    R_B_to_W(3, 2) = sphi * ctheta;
    R_B_to_W(3, 3) = cphi * ctheta;
end
