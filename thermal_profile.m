function x = thermal_profile(H1_val, H2_val)

res = 2000;
x = linspace(0, 1, res);
H1 = ones(1,res)*H1_val;
H2 = ones(1,res)*H2_val;
y1 = arrayfun(@thermal_profile_perturb_activ_fug, x, H1, zeros(1, res));
y2 = arrayfun(@thermal_profile_perturb_activ_fug, x, H2, ones(1, res));

figure; hold on;
plot(x, log10(y1), "red");
plot(x, log10(y2), "blue");
xlabel("Nondimensionalized temp");
ylabel("Log(MU_{temp perturb})");
legend('dry', 'wet');

end