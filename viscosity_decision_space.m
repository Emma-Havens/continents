function V_space = viscosity_decision_space(H_min, H_max)

res = 200;
% res_step = (H_max - H_min) / res;
% H_step = H_min;
temp_space = linspace(0, 1, res); % temperature
H_space = linspace(H_min, H_max, res);
fug_space = linspace(1, 0, res);

V_space = log10(reshape(cell2mat(arrayfun(@v_row, temp_space, 'UniformOutput', false)), res, res));

%disp(V_space);

figure;
image(temp_space, fug_space, V_space, 'CDataMapping','scaled');
c = colorbar;
c.Label.String = 'Log(Viscosity)';
clim([-2 2.5]);
xlabel("Nondimensionalized temp");
ylabel("Nondimensionalized water fugacity");

figure; hold on;
colormap jet;
cmap = colormap;
cmap_res = floor(size(cmap,1) / res);
xlabel("Nondimensionalized temp");
xlim([0 1]);
ylabel("Log(Nondimensionalized Viscosity)");
ylim([-2 3]);
for i = 1:res
    plot(temp_space, V_space(i, :), 'Color', cmap(i*cmap_res, :));
end

% for i = 1:res
% 
%     H = ones(1,res)*H_step;
%     V_space(i, :) = arrayfun(@thermal_profile_perturb_activ_fug, temp_space, H, zeros(1, res));
%     y2 = arrayfun(@thermal_profile_perturb_activ_fug, temp_space, H2, ones(1, res));
% 
%     H_step = H_step + res_step;
% end
% 
% figure; hold on;
% plot(temp_space, log10(y1), "red");
% plot(temp_space, log10(y2), "blue");
% xlabel("Nondimensionalized temp");
% ylabel("Log(MU_{temp perturb})");
% legend('dry', 'wet');

function row = v_row(temp)

row = arrayfun(@thermal_profile_perturb_activ_fug, temp*ones(1, res), H_space, fug_space);

end

end