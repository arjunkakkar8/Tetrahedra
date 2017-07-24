% Determinant that must vanish for dihedral angles
% We use x = [theta12, theta13, theta14, theta23, theta24, theta34]
mat = @(x) [
    -1, cos(x(6)), cos(x(5)), cos(x(4));
    cos(x(6)), -1, cos(x(3)), cos(x(2));
    cos(x(5)), cos(x(3)), -1, cos(x(1));
    cos(x(4)), cos(x(2)), cos(x(1)), -1];

matdet = @(x) det(mat(x));

% Three dihedral angles 2pi/a, 2pi/b, 2pi/c must sum to more than pi
violate3 = @(a, b, c) 2*(a*b + b*c + c*a) <= a*b*c;

% Four dihedral angles 2pi/a, 2pi/b, 2pi/c, 2pi/d must sum to less than 2pi
violate4 = @(a, b, c, d) a*b*c + b*c*d + c*d*a + d*a*b >= a*b*c*d;

result = [];

% Search with theta_ij = 2pi/n_ij
% where 3 <= n_ij <= 41

% We enforce the following symmetries:
% - n12 is the least among all n_ij
% - n13 is the least among n13, n14, n23, n24
% - if n12 = n13, then n24 <= n34
% - if n13 = n24, then n14 <= n23
% - if n13 = n14, then n23 <= n24
for n12 = 3:41
    for n13 = n12:41
        for n14 = n13:41
            if violate3(n12, n13, n14)
                continue 
            end
            
            for n23 = n13:41
                for n24 = n13:41
                    if violate3(n12, n23, n24)
                        continue
                    end
                    
                    for n34 = n12:41
                        if n12 == n13 && n24 > n34
                            continue
                        end
                        if n13 == n24 && n14 > n23
                            continue
                        end
                        if n13 == n14 && n23 > n24
                            continue
                        end
                        
                        if violate3(n13, n23, n34) || violate3(n14, n24, n34)
                            continue
                        end
                        
                        if violate4(n12, n34, n13, n24) || ...
                           violate4(n12, n34, n14, n23) || ...
                           violate4(n13, n24, n14, n23)
                            continue
                        end
                        
                        l = [n12, n13, n14, n23, n24, n34];
                        
                        if abs(matdet(2*pi ./ l)) < 1e-9
                            result = [result; l];
                        end
                    end
                end
            end
        end
    end
end

result

areas = [];

% Compute face areas from dihedral angles
for i=1:length(result)
    [V D] = eig(mat(2*pi./result(i,:)));
    areas = [areas; vpa(V(:, 4)')];
end

areas