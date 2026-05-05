function [r_dN, v_dN] = hill2rv(r_cN, v_cN, rho_H, rhoP_H)
    [ON,~,omega_HN_H] = N2H(r_cN, v_cN);
    r_dN = r_cN + ON' * rho_H;
    v_dN = v_cN + ON' * (rhoP_H + cross(omega_HN_H,rho_H));
    
    function [ON,h,omega_HN_H] = N2H(r_cN, v_cN)
        h_vec = cross(r_cN, v_cN);
        h = norm(h_vec);
    
        o_cap_r = r_cN/norm(r_cN);
        o_cap_h = h_vec/h;
        o_cap_theta = cross(o_cap_h, o_cap_r);
    
        ON = [o_cap_r';o_cap_theta';o_cap_h'];
        fDot = h/(norm(r_cN)^2);
        omega_HN_H = [0,0,fDot]';
    end
end