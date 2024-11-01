function [species_structs, gupta_structs, gupta_mu_structs, gupta_kappa_structs] = transport()
N_blottner  = struct("A", 0.0115572, "B", 0.6031679,  "C",-12.4327495);
O_blottner  = struct("A", 0.0203144, "B", 0.4294404,  "C", -11.6021403);
NO_blottner = struct("A", 0.0436378, "B", -0.0335511, "C", -9.5767430);
N2_blottner = struct("A", 0.0268142, "B", 0.3177838,  "C", -11.3155513);
O2_blottner = struct("A", 0.0449290, "B", -0.0826158, "C", -9.2019475);

guptaCoeffs = [
    0.0112 1.6182 -11.3091
    0.0465 0.9271 -8.1137
    0.0410 1.0023 -8.3597
    0.0195 1.4880 -10.3654
    0.0179 1.4848 -10.2810
    0.0033 1.5572 -11.1616
    0.0140 1.5824 -10.8819
    0.0226 1.3700 -9.6631
    -0.0048 1.9195 -11.9261
    0.0034 1.5572 -11.1729
    0.0291 1.2676 -9.6878
    0.0438 0.9647 -8.2380
    0.0185 1.4882 -10.3301
    0.0179 1.4848 -10.3155
    0.0364 1.1176 -8.9695
];

N2_N2 = struct("A", 0.0, "B", guptaCoeffs(1,1), "C", guptaCoeffs(1,2),  "D",guptaCoeffs(1,3));
O2_N2 = struct("A", 0.0, "B", guptaCoeffs(2,1), "C", guptaCoeffs(2,2),  "D",guptaCoeffs(2,3));
O2_O2 = struct("A", 0.0, "B", guptaCoeffs(3,1), "C", guptaCoeffs(3,2),  "D",guptaCoeffs(3,3));
N_N2  = struct("A", 0.0, "B", guptaCoeffs(4,1), "C", guptaCoeffs(4,2),  "D",guptaCoeffs(4,3));
N_O2  = struct("A", 0.0, "B", guptaCoeffs(5,1), "C", guptaCoeffs(5,2),  "D",guptaCoeffs(5,3));
N_N   = struct("A", 0.0, "B", guptaCoeffs(6,1), "C", guptaCoeffs(6,2),  "D",guptaCoeffs(6,3));
O_N2  = struct("A", 0.0, "B", guptaCoeffs(7,1), "C", guptaCoeffs(7,2),  "D",guptaCoeffs(7,3));
O_O2  = struct("A", 0.0, "B", guptaCoeffs(8,1), "C", guptaCoeffs(8,2),  "D",guptaCoeffs(8,3));
O_N   = struct("A", 0.0, "B", guptaCoeffs(9,1), "C", guptaCoeffs(9,2),  "D",guptaCoeffs(9,3));
O_O   = struct("A", 0.0, "B", guptaCoeffs(10,1),"C",  guptaCoeffs(10,2),"D",  guptaCoeffs(10,3));
NO_N2 = struct("A", 0.0, "B", guptaCoeffs(11,1), "C", guptaCoeffs(11,2), "D", guptaCoeffs(11,3));
NO_O2 = struct("A", 0.0, "B", guptaCoeffs(12,1), "C", guptaCoeffs(12,2), "D", guptaCoeffs(12,3));
NO_N  = struct("A", 0.0, "B", guptaCoeffs(13,1), "C", guptaCoeffs(13,2), "D", guptaCoeffs(13,3));
NO_O  = struct("A", 0.0, "B", guptaCoeffs(14,1), "C", guptaCoeffs(14,2), "D", guptaCoeffs(14,3));
NO_NO = struct("A", 0.0, "B", guptaCoeffs(15,1), "C", guptaCoeffs(15,2), "D",guptaCoeffs(15,3));


% For Blottner viscosity
species_structs = {N_blottner, O_blottner, NO_blottner, N2_blottner, O2_blottner};
% For Gupta Dij fits
gupta_structs = {
    N_N  O_N  NO_N  N_N2  N_O2
    O_N  O_O  NO_O  O_N2  O_O2
    NO_N NO_O NO_NO NO_N2 NO_O2
    N_N2 O_N2 NO_N2 N2_N2 O2_N2
    N_O2 O_O2 NO_O2 O2_N2 O2_O2
};

% data_mu= CSV.File("/Users/rloekvh/Downloads/tabula-gupta-mu.csv", header=true) |> Tables.matrix
% data_mu = data_mu[[3,4,5,1,2],2:end]
% data_kappa= CSV.File("/Users/rloekvh/Downloads/tabula-gupta-kappa.csv", header=true) |> Tables.matrix
% data_kappa = data_kappa[[3,4,5,1,2],2:end]
data_mu = [ 
0.012    0.593   -12.3805
0.0205   0.4257  -11.5803
0.0452  -0.0609   -9.4596
0.0203   0.4329  -11.8153
0.0484  -0.1455   -8.9231];

data_kappa = [
    0.0       0.0       0.01619     0.55022  -12.9219
    0.0       0.0       0.0331      0.22834  -11.5812
    0.02792  -0.87133  10.1797    -52.0347    88.6706
    0.03607  -1.07503  11.9503    -57.9006    93.2178
    0.07987  -2.58428  31.2596   -166.763    321.698
];

N_mu_struct  = struct("A", data_mu(1, 1), "B", data_mu(1, 2), "C", data_mu(1, 3));
O_mu_struct  = struct("A", data_mu(2, 1), "B", data_mu(2, 2), "C", data_mu(2, 3));
NO_mu_struct = struct("A", data_mu(3, 1), "B", data_mu(3, 2), "C", data_mu(3, 3));
N2_mu_struct = struct("A", data_mu(4, 1), "B", data_mu(4, 2), "C", data_mu(4, 3));
O2_mu_struct = struct("A", data_mu(5, 1), "B", data_mu(5, 2), "C", data_mu(5, 3));
gupta_mu_structs = {N_mu_struct, O_mu_struct, NO_mu_struct, N2_mu_struct, O2_mu_struct};

N_kappa_struct  = struct("A", data_kappa(1, 1), "B", data_kappa(1, 2), "C", data_kappa(1, 3), "D", data_kappa(1, 4), "E", data_kappa(1, 5));
O_kappa_struct  = struct("A", data_kappa(2, 1), "B", data_kappa(2, 2), "C", data_kappa(2, 3), "D", data_kappa(2, 4), "E", data_kappa(2, 5));
NO_kappa_struct = struct("A", data_kappa(3, 1), "B", data_kappa(3, 2), "C", data_kappa(3, 3), "D", data_kappa(3, 4), "E", data_kappa(3, 5));
N2_kappa_struct = struct("A", data_kappa(4, 1), "B", data_kappa(4, 2), "C", data_kappa(4, 3), "D", data_kappa(4, 4), "E", data_kappa(4, 5));
O2_kappa_struct = struct("A", data_kappa(5, 1), "B", data_kappa(5, 2), "C", data_kappa(5, 3), "D", data_kappa(5, 4), "E", data_kappa(5, 5));
gupta_kappa_structs = {N_kappa_struct, O_kappa_struct, NO_kappa_struct, N2_kappa_struct, O2_kappa_struct};

end
% 5×4 Matrix{Any}:
%  "N2"  0.0203   0.4329  -11.8153
%  "O2"  0.0484  -0.1455   -8.9231
%  "N"   0.012    0.593   -12.3805
%  "O"   0.0205   0.4257  -11.5803
%  "NO"  0.0452  -0.0609   -9.4596



% % WhAT IS ETA_i in mutationpp
% Ptest = mix.standardStateP()
% % Ptest = 1000.0
% Trange = 300.0:50:12000.0
% mu_plt = zeros(length(Trange), 5)
% lam_plt = zeros(length(Trange), 5)
% Dij_plt = zeros(length(Trange), 5)
% Dij_plt2 = zeros(length(Trange), 5)
% Dij_mpp = zeros(length(Trange), 5)
% mu_tot_mpp = zeros(length(Trange))
% lambda_tot_mpp = zeros(length(Trange))
% mu_tot_jl = zeros(length(Trange))
% lambda_tot_jl = zeros(length(Trange))
% Mw = mix.speciesMw()
% for i = 1:length(Trange)
%     % mix.equilibrate(Trange[i], Ptest)
%     mix.setState(rho_inf * Y_vec, Trange[i], 1)
%     Ptest = mix.P()
%     mu_curr = speciesViscosities(Trange[i], gupta_mu_structs)
%     % lambda_curr = mu_curr * 3/2 .* getCpsMass(Trange[i], mix.speciesMw())
%     % lambda_curr = 3.75 * mu_curr * RU ./ Mw
%     lambda_curr = speciesViscosities(Trange[i], gupta_kappa_structs)
%     phi_curr = euckenPhi(mu_curr, Mw, mix.X())
%     mu_plt[i,:] = mu_curr
%     lam_plt[i,:] = lambda_curr
%     Dij_plt[i,:] = averageDiffusionCoeffs(mix.T(), mix.X(), mix.Y(), mix.speciesMw(), mix.P())
%     Dij_mpp[i,:] = mix.averageDiffusionCoeffs()
%     Dij_plt2[i,:] = averageDiffusionCoeffs(mix.T(), mix.X(), mix.X(), mix.speciesMw(), mix.P())

%     mu_tot_mpp[i] = mix.viscosity()
%     mu_tot_jl[i] = viscosityWilke(mu_curr, mix.X(), phi_curr)
%     lambda_tot_mpp[i] = mix.frozenThermalConductivity()
%     lambda_tot_jl[i] = viscosityWilke(lambda_curr, mix.X(), phi_curr)
% end
% % plot(Trange, Dij_plt)
% iplt = 5
% % plot(Trange, Dij_plt[:,iplt], xlims=(0.0, 5000.0), ylims=(0, 0.004))
% % % plot!(Trange, Dij_plt2[:,iplt], xlims=(0.0, 5000.0), ylims=(0, 0.004))
% % plot!(Trange, Dij_mpp[:,iplt], xlims=(0.0, 5000.0), ylims=(0, 0.004))
% plot(Trange, Dij_plt[:,iplt])
% % plot!(Trange, Dij_plt2[:,iplt], xlims=(0.0, 5000.0), ylims=(0, 0.004))
% plot!(Trange, Dij_mpp[:,iplt])

% plot(Trange, Dij_mpp, xlims=(0.0, 5000.0), ylims=(0, 0.004))
% plot(Trange, mu_plt, xlims=(300, 5050.0), ylims=(0.0, 0.0002))
% plot(Trange, lam_plt, xlims=(300, 5050.0), ylims=(0.0, 0.5))

% plot(Trange, mu_tot_mpp)
% plot!(Trange, mu_tot_jl)

% plot(Trange, lambda_tot_mpp)
% plot!(Trange, lambda_tot_jl)

% plot(Trange, abs.((mu_tot_jl - mu_tot_mpp))./mu_tot_mpp)
% plot!(Trange, abs.((lambda_tot_jl - lambda_tot_mpp))./lambda_tot_mpp)

% #Mu seems alright but lambda seems a bit off....also i feel like it should be UNDERShooting lambda but it's overshooting
% #mu: <4%
% #kappa: up to around 8% toward the higher ranges of temperature
    % but the weird thing about kappa is that i feel like it should not be overshooting
    % well I figured out a certain formluation that is more consistent with the code AND the Gupta paper. 
    % The errors are worse but it's undershooting by a lot...i wonder if the internal cp thing can make up the difference. Let's comment it out! 