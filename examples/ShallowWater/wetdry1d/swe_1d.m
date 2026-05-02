function swe_1d
% 1D shallow-water on a sloped beach.
%   Left  boundary : reflecting wall (drying intertidal).
%   Right boundary : prescribed sinusoidal free-surface (tidal forcing).
%   Scheme         : finite volume, conservative Rusanov flux with
%                    Audusse (2004) hydrostatic reconstruction,
%                    SSP-RK2 (Heun) in time.
%   Lake-at-rest and wet/dry over bed steps are handled via interface
%   reconstruction; no separate bed source term is needed.

    g = 9.81;

    % -------- domain & bathymetry --------
    L  = 10000;
    N  = 400;
    dx = L/N;
    x  = (0.5*dx : dx : L-0.5*dx).';

    % bed elevation zb(x): +3 m at x=0 (above MSL, dries) -> -7 d m at x=L (deep)
    z_beach = 3.0;  z_deep = -7.0;
    zb = z_beach + (z_deep - z_beach) * x / L;

    % -------- forcing --------
    A_tide = 0.4;                     % m  (small forcing)
    T_tide = 1800;                    % s
    omega  = 2*pi / T_tide;
    eta_bc = @(t) A_tide * sin(omega*t);

    % -------- numerics --------
    H_tol = 1e-4;                     % wet/dry threshold
    CFL   = 0.4;
    t_end = 8 * T_tide;               % longer run

    % -------- internal weir as a rectangular berm --------
    has_weir       = true;           % <-- flip to true for the weir case
    x_weir         = 4000;            % m — berm center
    z_c            = 0.2;             % crest elevation [m]
    weir_halfwidth = 75;              % m — full width 150 m
    if has_weir
        in_weir     = abs(x - x_weir) <= weir_halfwidth;
        zb(in_weir) = z_c;            % flat top with vertical sides
    end
    iw_face  = round(x_weir/dx) + 1;  % center face index — proxy for crest flux
    x_weir   = (iw_face-1)*dx;        % snap label to face

    % case tag drives figure and movie filenames (swe_noweir_* vs swe_weir_*)
    if has_weir, case_tag = 'weir'; else, case_tag = 'noweir'; end

    % -------- live plot & movie controls --------
    live_plot   = true;
    plot_stride = 20;                 % update every N steps
    make_movie    = false;
    movie_profile = 'Motion JPEG AVI';                 % robust; 'MPEG-4' for smaller files
    movie_file    = sprintf('swe_%s.avi', case_tag);   % per-case filename
    movie_fps     = 30;
    if has_weir
        zoom_xlim = [2.5 4.8];        % km — covers shoreline + weir
    else
        zoom_xlim = [2.0 4.0];        % km — centered on shoreline motion
    end
    zoom_ylim   = [-1.5 1.5];         % m

    % -------- initial state: still water at MSL (eta = 0) --------
    H  = max(0 - zb, H_tol);
    Hu = zeros(N,1);

    % -------- diagnostics storage (grow, trim at end) --------
    cap = 200000;
    times       = zeros(cap,1);
    mass_hist   = zeros(cap,1);
    mom_hist    = zeros(cap,1);
    max_eta_h   = zeros(cap,1);
    min_eta_h   = zeros(cap,1);
    shore_hist  = zeros(cap,1);
    flux_right  = zeros(cap,1);
    flux_weir   = zeros(cap,1);
    dt_hist     = zeros(cap,1);
    cfl_hist    = zeros(cap,1);

    figure(1); clf; plot(x, H); 
    surf0 = H + zb; %max(H + zb, zb);
    figure(2); clf; plot(x, surf0); 

    % -------- live figure (full view + shoreline zoom) --------
    if live_plot
        hfig = figure('Name','SWE live','Position',[80 80 1100 700], ...
                      'Resize','off','Color','w');
        set(hfig,'PaperPositionMode','auto');
        Xw    = [x; flipud(x)] / 1000;
        surf0 = max(H + zb, zb);
        Yw    = [zb; flipud(surf0)];

        % --- full domain ---
        ax1 = subplot(2,1,1,'Parent',hfig); hold(ax1,'on');
        fill(ax1, [x; flipud(x)]/1000, [zb; -10*ones(N,1)], ...
             [0.55 0.27 0.07], 'EdgeColor','none','FaceAlpha',0.6);
        plot(ax1, x/1000, zb, 'k-','LineWidth',1);
        yline(ax1, 0,'--','Color',[0.5 0.5 0.5]);
        hwat1 = fill(ax1, Xw, Yw, [0.27 0.51 0.71], ...
                     'EdgeColor',[0.1 0.25 0.45],'FaceAlpha',0.75,'LineWidth',0.6);
        hbc   = plot(ax1, L/1000, eta_bc(0),'ro','MarkerFaceColor','r');
        % rectangle showing the zoom window on the full view
        rectangle(ax1,'Position',[zoom_xlim(1) zoom_ylim(1) ...
                    diff(zoom_xlim) diff(zoom_ylim)], ...
                    'EdgeColor',[0.8 0.1 0.1],'LineStyle','--','LineWidth',1);
        if has_weir
            plot(ax1,[x_weir-weir_halfwidth x_weir+weir_halfwidth]/1000, ...
                 [z_c z_c],':','Color',[0.2 0.2 0.2],'LineWidth',1);
        end
        xlim(ax1,[0 L/1000]); ylim(ax1,[-9 3]);
        xlabel(ax1,'x [km]'); ylabel(ax1,'elevation [m]');
        htitle = title(ax1, sprintf('t = 0.000 T   (A = %.2f m, T = %.0f s)', A_tide, T_tide));
        grid(ax1,'on');

        % --- zoom on moving shoreline ---
        ax2 = subplot(2,1,2,'Parent',hfig); hold(ax2,'on');
        fill(ax2, [x; flipud(x)]/1000, [zb; -10*ones(N,1)], ...
             [0.55 0.27 0.07], 'EdgeColor','none','FaceAlpha',0.6);
        plot(ax2, x/1000, zb, 'k-','LineWidth',1.2);
        yline(ax2, 0,'--','Color',[0.5 0.5 0.5]);
        hwat2    = fill(ax2, Xw, Yw, [0.27 0.51 0.71], ...
                        'EdgeColor',[0.1 0.25 0.45],'FaceAlpha',0.75,'LineWidth',0.8);
        hshore2  = plot(ax2, NaN, NaN,'o','Color',[0.85 0.1 0.1], ...
                        'MarkerFaceColor',[0.85 0.1 0.1],'MarkerSize',6);
        if has_weir
            yline(ax2, z_c,':','Color',[0.2 0.2 0.2],'LineWidth',0.8);
        end
        xlim(ax2, zoom_xlim); ylim(ax2, zoom_ylim);
        xlabel(ax2,'x [km]'); ylabel(ax2,'elevation [m]');
        title(ax2,sprintf('shoreline + weir (z_c = %.2f m)', z_c)); grid(ax2,'on');
        drawnow;

        % --- movie writer ---
        if make_movie
            vw = VideoWriter(fullfile(pwd,movie_file), movie_profile);
            vw.FrameRate = movie_fps;
            if isprop(vw,'Quality'), vw.Quality = 90; end
            open(vw);
        end
    end

    pause

    % -------- time loop --------
    t = 0;  k = 0;
    while t < t_end
        Hs   = max(H, H_tol);
        u    = Hu ./ Hs;
        c    = sqrt(g*Hs);
        vmax = max(abs(u) + c);
        dt   = CFL * dx / vmax;
        if t + dt > t_end, dt = t_end - t; end

        % if mod(k, 30) == 0
        %   figure(1); clf; plot(x, u); pause(0.01);
        % end
       
        % SSP-RK2 (Heun):  U_new = U + 0.5*dt*(L(U,t) + L(U+dt*L(U,t), t+dt))
        [dH1, dHu1, FH1] = rhs_step(H, Hu, t);
        H_s  = H  + dt*dH1;
        Hu_s = Hu + dt*dHu1;
        dry_s          = H_s < H_tol;
        H_s(dry_s)     = H_tol;
        Hu_s(dry_s)    = 0;

        [dH2, dHu2, FH2] = rhs_step(H_s, Hu_s, t + dt);
        H_new  = H  + 0.5*dt*(dH1  + dH2);
        Hu_new = Hu + 0.5*dt*(dHu1 + dHu2);

        dry          = H_new < H_tol;
        H_new(dry)   = H_tol;
        Hu_new(dry)  = 0;

        % step-averaged mass flux for diagnostics
        F_H = 0.5*(FH1 + FH2);

        H = H_new;  Hu = Hu_new;
        t = t + dt;
        k = k + 1;

        % record
        times(k)      = t;
        mass_hist(k)  = sum(H)*dx;
        mom_hist(k)   = sum(Hu)*dx;
        eta_now       = H + zb;
        wet           = H > 2*H_tol;
        if any(wet)
            max_eta_h(k) = max(eta_now(wet));
            min_eta_h(k) = min(eta_now(wet));
            shore_hist(k) = x(find(wet,1,'first'));
        else
            max_eta_h(k) = NaN; min_eta_h(k) = NaN; shore_hist(k) = NaN;
        end
        flux_right(k) = F_H(N+1);
        if has_weir, flux_weir(k) = F_H(iw_face); end
        dt_hist(k)    = dt;
        cfl_hist(k)   = vmax*dt/dx;

        if live_plot && mod(k, plot_stride) == 0
            surf_now = H + zb;
            surf_now(H <= 2*H_tol) = zb(H <= 2*H_tol);
            Ynow = [zb; flipud(surf_now)];
            set(hwat1, 'YData', Ynow);
            set(hwat2, 'YData', Ynow);
            set(hbc,   'YData', eta_bc(t));
            if ~isnan(shore_hist(k))
                xs = shore_hist(k)/1000;
                ys = interp1(x/1000, zb, xs, 'linear');
                set(hshore2,'XData',xs,'YData',ys);
            end
            set(htitle,'String', sprintf('t = %.3f T   (step %d,  CFL_{max} = %.2f)', ...
                t/T_tide, k, cfl_hist(k)));
            if make_movie
                drawnow;                                   % full flush before capture
                img = print(hfig,'-RGBImage','-r100');     % deterministic size
                % force even height & width (H.264 requirement; harmless for MJPEG)
                if mod(size(img,1),2), img = img(1:end-1,:,:); end
                if mod(size(img,2),2), img = img(:,1:end-1,:); end
                if ~exist('frame_size','var')
                    frame_size = [size(img,1) size(img,2)];
                elseif any([size(img,1) size(img,2)] ~= frame_size)
                    img = imresize(img, frame_size);       % guard against HiDPI drift
                end
                writeVideo(vw, img);
            else
                drawnow limitrate;
            end
        end
    end
    if live_plot
        drawnow;
        if make_movie, close(vw); end
    end

    % trim
    times=times(1:k); mass_hist=mass_hist(1:k); mom_hist=mom_hist(1:k);
    max_eta_h=max_eta_h(1:k); min_eta_h=min_eta_h(1:k);
    shore_hist=shore_hist(1:k); flux_right=flux_right(1:k);
    flux_weir=flux_weir(1:k);
    dt_hist=dt_hist(1:k); cfl_hist=cfl_hist(1:k);

    cum_flux_right = cumsum(flux_right .* dt_hist);

    if live_plot, saveas(hfig, fullfile(pwd, sprintf('swe_%s_final.png', case_tag))); end

    % -------- diagnostics figure --------
    figure('Name','SWE diagnostics','Position',[120 120 1000 650]);

    subplot(2,2,1);
    M0 = mass_hist(1);
    plot(times/T_tide, mass_hist - M0); hold on;
    plot(times/T_tide, -cum_flux_right, '--');
    legend('M(t)-M(0)','-\int F_{right} dt','Location','best');
    xlabel('t / T'); ylabel('mass [m^2]'); title('Mass balance'); grid on;

    subplot(2,2,2);
    plot(times/T_tide, (mass_hist - M0) + cum_flux_right);
    xlabel('t / T'); ylabel('residual [m^2]');
    title('Mass residual  (>0: H_{tol} clipping added mass)'); grid on;

    subplot(2,2,3);
    plot(times/T_tide, max_eta_h); hold on;
    plot(times/T_tide, min_eta_h);
    plot(times/T_tide, A_tide*sin(omega*times),'k--','LineWidth',0.6);
    legend('max \eta','min \eta','forcing','Location','best');
    xlabel('t / T'); ylabel('\eta [m]'); title('Free-surface extrema'); grid on;

    subplot(2,2,4);
    if has_weir
        yyaxis left
        plot(times/T_tide, flux_weir);
        ylabel('q_{weir} [m^2/s]');
        yyaxis right
        plot(times/T_tide, shore_hist/1000);
        ylabel('x_{shore} [km]');
        xlabel('t / T');
        title('Weir overtopping flux (L) & shoreline (R)'); grid on;
    else
        plot(times/T_tide, shore_hist/1000);
        xlabel('t / T'); ylabel('x_{shore} [km]');
        title('Leftmost wet cell (shoreline excursion)'); grid on;
    end

    saveas(gcf, fullfile(pwd, sprintf('swe_%s_diag.png', case_tag)));

    % -------- summary --------
    fprintf('steps            : %d\n', k);
    fprintf('final time       : %.1f s  (%.3f T)\n', times(end), times(end)/T_tide);
    fprintf('mean dt          : %.3f s\n', mean(dt_hist));
    fprintf('max CFL observed : %.3f\n',  max(cfl_hist));
    fprintf('initial mass     : %.4f m^2\n', M0);
    fprintf('final mass       : %.4f m^2\n', mass_hist(end));
    fprintf('Dmass            : %+.4e m^2\n', mass_hist(end)-M0);
    fprintf('int F_right dt   : %+.4e m^2\n', cum_flux_right(end));
    fprintf('mass residual    : %+.4e m^2  (clipping artefact)\n', ...
            (mass_hist(end)-M0) + cum_flux_right(end));

    % -------- nested RHS: spatial discretization (Audusse + Rusanov_cons) ------
    function [dH_r, dHu_r, FH_r] = rhs_step(H_in, Hu_in, t_bc)
        eta_r    = H_in + zb;
        zst_r    = max(zb(1:N-1), zb(2:N));
        HL_r     = max(0, eta_r(1:N-1) - zst_r);
        HR_r     = max(0, eta_r(2:N)   - zst_r);
        uL_r     = Hu_in(1:N-1) ./ max(H_in(1:N-1), H_tol);
        uR_r     = Hu_in(2:N)   ./ max(H_in(2:N),   H_tol);
        HuL_r    = HL_r .* uL_r;
        HuR_r    = HR_r .* uR_r;

        [FHi_r, FHui_r] = rusanov_cons(HL_r, HuL_r, HR_r, HuR_r, g, H_tol);
        FHuL_i_r = FHui_r + 0.5*g*(H_in(1:N-1).^2 - HL_r.^2);
        FHuR_i_r = FHui_r + 0.5*g*(H_in(2:N).^2   - HR_r.^2);

        FH_r   = zeros(N+1,1);
        FHuL_r = zeros(N+1,1);
        FHuR_r = zeros(N+1,1);
        FH_r(2:N)   = FHi_r;
        FHuL_r(2:N) = FHuL_i_r;
        FHuR_r(2:N) = FHuR_i_r;

        % left wall (mirror; trivial Audusse — same zb)
        [Fhw_r, Fhuw_r] = rusanov_cons(H_in(1), -Hu_in(1), H_in(1), Hu_in(1), g, H_tol);
        FH_r(1)   = Fhw_r;
        FHuR_r(1) = Fhuw_r;

        % right open BC
        eta_g_r = eta_bc(t_bc);
        H_g_r   = max(eta_g_r - zb(end), H_tol);
        u_int_r = Hu_in(end) / max(H_in(end), H_tol);
        Hu_g_r  = H_g_r * u_int_r;
        [Fhr_r, Fhur_r] = rusanov_cons(H_in(end), Hu_in(end), H_g_r, Hu_g_r, g, H_tol);
        FH_r(N+1)   = Fhr_r;
        FHuL_r(N+1) = Fhur_r;

        dH_r  = -(FH_r(2:N+1)   - FH_r(1:N))   / dx;
        dHu_r = -(FHuL_r(2:N+1) - FHuR_r(1:N)) / dx;
    end
end

% ------------------------------------------------------------------------
function [F_H, F_Hu] = rusanov(HL, HuL, HR, HuR, g, H_tol)
    HLs = max(HL, H_tol);  HRs = max(HR, H_tol);
    uL  = HuL ./ HLs;      uR  = HuR ./ HRs;
    cL  = sqrt(g*HLs);     cR  = sqrt(g*HRs);
    a   = max(abs(uL)+cL, abs(uR)+cR);
    F_H  = 0.5*(HuL       + HuR)       - 0.5*a.*(HR  - HL);
    F_Hu = 0.5*(HuL.*uL   + HuR.*uR)   - 0.5*a.*(HuR - HuL);
end

function [F_H, F_Hu] = rusanov_cons(HL, HuL, HR, HuR, g, H_tol)
    % conservative-form Rusanov: F_Hu = Hu u + 0.5 g H^2 (pressure in the flux)
    HLs  = max(HL, H_tol);  HRs = max(HR, H_tol);
    uL   = HuL ./ HLs;      uR  = HuR ./ HRs;
    cL   = sqrt(g*HLs);     cR  = sqrt(g*HRs);
    a    = max(abs(uL)+cL, abs(uR)+cR);
    fHuL = HuL.*uL + 0.5*g*HL.^2;
    fHuR = HuR.*uR + 0.5*g*HR.^2;
    F_H  = 0.5*(HuL  + HuR)  - 0.5*a.*(HR  - HL);
    F_Hu = 0.5*(fHuL + fHuR) - 0.5*a.*(HuR - HuL);
end
