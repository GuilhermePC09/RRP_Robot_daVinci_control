function J = GeoJac( Robot )
%GEOJAC Computes the Geometric Jacobian

num_joints = size(Robot,1);      %get number of robot joints

T = sym(eye(4));
T_matrices = cell(1, num_joints);

for i = 1:num_joints
    T_i = DHTransf(Robot(i,:));
    T = T * T_i;
    T_matrices{i} = T;
end

% Compute Jacobian
J = sym(zeros(6, num_joints));
for i = 1:num_joints
    if i == 1
        T_prev = sym(eye(4));
    else
        T_prev = T_matrices{i-1};
    end

    z = T_prev(1:3, 3);
    o_n = T_matrices{end}(1:3, 4);
    o_i = T_prev(1:3, 4);
    
    p = Robot(i,:);

    if ~isempty(symvar(p(2))) % rotational joint 
        Jv = cross(z, (o_n - o_i));
        Jw = z;
    else 
        Jv = z;
        Jw = sym([0; 0; 0]);
    end

    J(:, i) = [Jv; Jw];

end


end

