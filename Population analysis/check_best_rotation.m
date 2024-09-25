function [out_dim_deg, rotated_ctx1_to_ctx2] = check_best_rotation(Ctx1, Ctx2)

% calculate cross and dot products
% find optimal rotation between contexts

get_angle = @(x,y) rad2deg((acos(dot(x,y) / (norm(x)*norm(y)))));

out_dim_deg =[];

% rotate Ctx1 to fit Ctx2
[ret_R, ret_t] = rigid_transform_3D(Ctx1', Ctx2');

% Compare the recovered rotation and t with the original
rotated_ctx1_to_ctx2 = ((ret_R*Ctx1') + repmat(ret_t, 1, size(Ctx2,1)))';

% figure; 
% hold on
% plot3(rotated_ctx1_to_ctx2(:,1),rotated_ctx1_to_ctx2(:,2), rotated_ctx1_to_ctx2(:,3), 'LineWidth', 3);
% plot3(Ctx2(:,1),Ctx2(:,2), Ctx2(:,3), 'LineWidth', 3);


% find out how many degrees of rotation were necessary
for i = 1:size(Ctx2,1)
    
    dim_deg(i) = get_angle(Ctx1(i,:),rotated_ctx1_to_ctx2(i,:));
    
end % of looping 

out_dim_deg = mean(dim_deg(1:2));

end % of function 