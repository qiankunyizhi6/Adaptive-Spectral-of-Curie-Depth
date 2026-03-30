import optuna
import numpy as np
import os
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt


def load_measured_heat_flow(measured_file):
    """Load the measured heat flow data"""
    data = np.loadtxt(measured_file)
    lons = data[:, 0]
    lats = data[:, 1]
    heat_flows = data[:, 2]  # The third column represents the heat flow value. (mW/m²)
    return lons, lats, heat_flows


def load_block_polygons(block_dir):
    block_polygons = {}
    block_files = [f for f in os.listdir(block_dir) if f.startswith('block_') and f.endswith('_fixed.txt')]
    for block_file in block_files:
        block_name = block_file.replace('_fixed.txt', '')
        file_path = os.path.join(block_dir, block_file)
        try:
            with open(file_path, 'r') as f:
                content = f.read().strip()
                coords = list(map(float, content.split()))
            polygon = []
            for i in range(0, len(coords), 2):
                if i + 1 < len(coords):
                    polygon.append((coords[i], coords[i + 1]))
            block_polygons[block_name] = polygon
        except Exception as e:
            print(f"Error occurred while reading {block_file}: {e}")
    return block_polygons


def point_in_polygon(point, polygon):
    x, y = point
    n = len(polygon)
    inside = False
    p1x, p1y = polygon[0]
    for i in range(1, n + 1):
        p2x, p2y = polygon[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside


def get_points_in_block(block_polygon, measured_lons, measured_lats, measured_heat_flows):
    block_lons, block_lats, block_heat_flows = [], [], []

    for lon, lat, hf in zip(measured_lons, measured_lats, measured_heat_flows):
        if point_in_polygon((lon, lat), block_polygon):
            block_lons.append(lon)
            block_lats.append(lat)
            block_heat_flows.append(hf)
    return np.array(block_lons), np.array(block_lats), np.array(block_heat_flows)


def find_block_for_point(point, block_polygons):
    for block_name, polygon in block_polygons.items():
        if point_in_polygon(point, polygon):
            return block_name
    return None


def detrend_2d(Z):
    N, M = Z.shape
    X, Y = np.meshgrid(np.arange(1, M + 1), np.arange(1, N + 1))
    Xcol, Ycol, Zcol = X.ravel(), Y.ravel(), Z.ravel()
    A = np.column_stack((Xcol, Ycol, np.ones_like(Xcol)))
    Coefficients, *_ = np.linalg.lstsq(A, Zcol, rcond=None)
    XCoeff, YCoeff, CCoeff = Coefficients
    Z_p = XCoeff * X + YCoeff * Y + CCoeff
    return Z - Z_p


def raas(anom, dl=5.5):
    ny, nx = anom.shape
    N = max(nx, ny)
    if N % 2 != 0:
        N += 1
    TFA = np.fft.fftshift(np.fft.fft2(anom, s=(N, N)))
    EA = 2 * np.abs(TFA) / (ny * nx)
    n = np.arange(-N // 2, N // 2)
    NX, NY = np.meshgrid(n, n)
    rho = np.hypot(NX, NY)
    rho = np.round(rho).astype(int)
    num_bins = N // 2 + 1
    Pf, CF, nm = np.zeros(num_bins), np.zeros(num_bins), np.zeros(num_bins)
    eps = 1e-10
    for r in range(num_bins):
        idx = np.where(rho == r)
        count = idx[0].size
        nm[r] = count
        if count > 0:
            log_EA = np.log(np.where(EA[idx] == 0, eps, EA[idx]))
            Pf[r] = np.nanmean(log_EA)
            CF[r] = 1.96 * np.nanstd(log_EA) / np.sqrt(count)
        else:
            Pf[r] = np.nan
            CF[r] = np.nan
    kn = 1 / (2 * dl)
    k0 = 1 / (N * dl)
    w = np.arange(0, kn + k0 / 2, k0)
    cv2 = np.column_stack((w[1:], Pf[1:], CF[1:], nm[1:]))
    cv3 = cv2[cv2[:, 2] <= 1]
    return cv3


def fractal_correction(x, b1):
    w = x[:, 0]
    ps = (x[:, 1] - np.log(2 * np.pi * w)) + np.log((2 * np.pi * w) ** (0.5 * b1))
    pw = x[:, 1] + np.log((2 * np.pi * w) ** (0.5 * b1))
    return w, pw, ps, x[:, 2]


def compute_curie_depth_and_heat_flow_for_measured_points(b2, window_size, data_file, myrange1, myrange2, measured_lons,
                                                          measured_lats):
    """Calculate the depth of the cavity and the corresponding heat flow value"""
    if len(measured_lons) == 0:
        return np.array([]), np.array([])

    # Heat flow calculation parameters
    Ts = 15  # Surface temperature(°C)
    Tc = 580  # Curie temperature(°C)
    k = 2.62  # Thermal conductivity(W/m·K)
    H0 = 0.00000285  # Heat output rate
    hr = 13800  # Heat output rate attenuation depth(m)

    mag_data = np.loadtxt(data_file)
    lon_grid = np.unique(mag_data[:, 0])
    lat_grid = np.unique(mag_data[:, 1])
    values = mag_data[:, 2]
    value_grid = np.full((len(lat_grid), len(lon_grid)), np.nan)

    for i in range(len(values)):
        lon_idx = np.argmin(np.abs(lon_grid - mag_data[i, 0]))
        lat_idx = np.argmin(np.abs(lat_grid - mag_data[i, 1]))
        value_grid[lat_idx, lon_idx] = values[i]

    pred_zbs = np.full(len(measured_lons), np.nan)
    pred_heat_flows = np.full(len(measured_lons), np.nan)

    for ipt, (lon0, lat0) in enumerate(zip(measured_lons, measured_lats)):
        lon_min = lon0 - window_size / 2
        lon_max = lon0 + window_size / 2
        lat_min = lat0 - window_size / 2
        lat_max = lat0 + window_size / 2

        lon_mask = (lon_grid >= lon_min) & (lon_grid <= lon_max)
        lat_mask = (lat_grid >= lat_min) & (lat_grid <= lat_max)

        if lon_mask.sum() < 2 or lat_mask.sum() < 2:
            continue

        lon_idx1, lon_idx2 = np.where(lon_mask)[0][[0, -1]]
        lat_idx1, lat_idx2 = np.where(lat_mask)[0][[0, -1]]
        window_data = value_grid[lat_idx1:lat_idx2 + 1, lon_idx1:lon_idx2 + 1]

        if window_data.size < 4 or np.isnan(window_data).all():
            continue

        try:
            dt = detrend_2d(window_data)
        except Exception:
            continue

        a2 = raas(dt, 5.5)
        if a2.size == 0:
            continue

        w2, pw2, ps2, ci1 = fractal_correction(a2, b2)
        low1, high1 = myrange1
        low2, high2 = myrange2
        if (high1 >= len(w2) or high2 >= len(w2) or
            low1 >= len(w2) or low2 >= len(w2) or
            low1 > high1 or low2 > high2):
            continue

        xs2 = w2[low1:high1 + 1]
        ys2 = pw2[low1:high1 + 1]
        xs = w2[low2:high2 + 1]
        ys = ps2[low2:high2 + 1]

        if len(xs2) < 2 or len(xs) < 2:
            continue

        try:
            p2 = np.polyfit(xs2, ys2, 1)  # Zt
            pl = np.polyfit(xs, ys, 1)  # Z0
        except Exception:
            continue

        Zt = -p2[0] / (2 * np.pi) 
        Z0 = -pl[0] / (2 * np.pi)  
        Zb = 2 * Z0 - Zt   

        Zb_m = 1000 * Zb
        Z0_m = 1000 * Z0

        term1 = k * (Tc - Ts) / (Zb_m - Z0_m)
        term2 = hr ** 2 * H0 * (np.exp(-Zb_m / hr) / Zb_m)
        term3 = hr * H0 * np.exp(-Z0_m / hr)
        qs = term1 + term2 + term3
        qs_mW = qs * 1000  # 转换为 mW/m²

        if 10 <= Zb <= 60 and 0 < qs_mW < 500:  
            pred_zbs[ipt] = Zb
            pred_heat_flows[ipt] = qs_mW

    return pred_zbs, pred_heat_flows


def objective(trial, b2, window_size, data_file, measured_lons, measured_lats, measured_heat_flows):
    myrange1_start = trial.suggest_int("myrange1_start", 0, 20)
    myrange1_end = trial.suggest_int("myrange1_end", myrange1_start + 1, 25)
    myrange2_start = trial.suggest_int("myrange2_start", 0, 20)
    myrange2_end = trial.suggest_int("myrange2_end", myrange2_start + 1, 25)
    myrange1 = [myrange1_start, myrange1_end]
    myrange2 = [myrange2_start, myrange2_end]

    try:
        pred_zbs, pred_heat_flows = compute_curie_depth_and_heat_flow_for_measured_points(
            b2, window_size, data_file, myrange1, myrange2, measured_lons, measured_lats
        )
        valid_mask = ~np.isnan(pred_heat_flows)

        if valid_mask.sum() == 0:
            raise optuna.TrialPruned("No valid points")

        valid_pred_heat_flows = pred_heat_flows[valid_mask]
        valid_measured_heat_flows = measured_heat_flows[valid_mask]

        if np.any(valid_pred_heat_flows < 0) or np.any(valid_pred_heat_flows > 500):
            raise optuna.TrialPruned("Abnormal heat flow values")

        rmse = mean_squared_error(valid_measured_heat_flows, valid_pred_heat_flows) ** 0.5

        return rmse

    except Exception as e:
        raise optuna.TrialPruned(f"Trial failed: {str(e)}")


def optimize_for_b2_ws(b2, window_size, data_file, measured_lons, measured_lats, measured_heat_flows, n_trials):

    if len(measured_lons) == 0:
        print("⚠️ No actual measurement points found. Skipping optimization.")
        return None, 1000.0
    successful_trials = 0
    failed_trials = 0

    def trial_callback(study, trial):
        nonlocal successful_trials, failed_trials
        if trial.state == optuna.trial.TrialState.COMPLETE:
            successful_trials += 1
        elif trial.state == optuna.trial.TrialState.PRUNED:
            failed_trials += 1

        if (successful_trials + failed_trials) % 10 == 0:
            print(f"    Progress progress: {successful_trials + failed_trials}/{n_trials} experiment "
                  f"(Success: {successful_trials}, Failure: {failed_trials})")

    study = optuna.create_study(
        direction="minimize",
        pruner=optuna.pruners.MedianPruner(n_startup_trials=5, n_warmup_steps=0, interval_steps=1)
    )

    try:
        study.optimize(
            lambda trial: objective(
                trial, b2, window_size, data_file, measured_lons, measured_lats, measured_heat_flows
            ),
            n_trials=n_trials,
            callbacks=[trial_callback]
        )
    except Exception as e:
        print(f"   An abnormality occurred during the optimization process.: {e}")
        return None, 1000.0

    if len(study.trials) == 0:
        print(f"    ⚠️ All the experiments failed and no effective parameters could be found.")
        return None, 1000.0

    best_params = study.best_params
    best_loss = study.best_value

    print(f"    ✅  {successful_trials} , {failed_trials} ")
    print(f"    Best RMSE: {best_loss:.4f} mW/m²")

    return best_params, best_loss


def compute_curie_depth_for_full_map(b2, window_size, data_file, block_params, block_polygons, output_dir,
                                     output_resolution=0.1):
    mag_data = np.loadtxt(data_file)

    lon_grid = np.unique(mag_data[:, 0])
    lat_grid = np.unique(mag_data[:, 1])
    values = mag_data[:, 2]

    lon_min, lon_max = 101, 117
    lat_min, lat_max = 32, 42
    lons_full = np.arange(lon_min, lon_max + output_resolution, output_resolution)
    lats_full = np.arange(lat_min, lat_max + output_resolution, output_resolution)

    print(f"Longitude {lons_full[0]:.1f}-{lons_full[-1]:.1f}, Latitude {lats_full[0]:.1f}-{lats_full[-1]:.1f}")
    print(f"Grid point count: {len(lons_full)} x {len(lats_full)} = {len(lons_full) * len(lats_full)}")

    lon_mesh, lat_mesh = np.meshgrid(lons_full, lats_full)
    zb_full = np.full(lon_mesh.shape, np.nan)

    value_grid = np.full((len(lat_grid), len(lon_grid)), np.nan)
    valid_data_points = 0

    for i in range(len(values)):
        lon_idx = np.argmin(np.abs(lon_grid - mag_data[i, 0]))
        lat_idx = np.argmin(np.abs(lat_grid - mag_data[i, 1]))
        if not np.isnan(values[i]):
            value_grid[lat_idx, lon_idx] = values[i]
            valid_data_points += 1

    print(f"Valid magnetic anomaly data points: {valid_data_points}/{len(values)}")
    print(f"Magnetic anomaly grid shape: {value_grid.shape}")

    block_coverage = {name: 0 for name in block_polygons.keys()}
    no_block_points = []
    valid_points = 0
    failed_points = 0
    failed_reasons = {
        'no_block': 0,
        'no_block_params': 0,
        'insufficient_window': 0,
        'no_data': 0,
        'detrend_failed': 0,
        'spectrum_failed': 0,
        'index_out_of_range': 0,
        'insufficient_freq_bands': 0,
        'fit_failed': 0,
        'constraint_failed': 0
    }

    total_points = len(lats_full) * len(lons_full)
    processed_points = 0

    for i in range(len(lats_full)):
        for j in range(len(lons_full)):
            processed_points += 1
            if processed_points % 1000 == 0:
                print(f"progress: {processed_points}/{total_points} ({processed_points / total_points * 100:.1f}%)")

            lon0, lat0 = lons_full[j], lats_full[i]
            block_name = find_block_for_point((lon0, lat0), block_polygons)

            if block_name is None:
                failed_points += 1
                failed_reasons['no_block'] += 1
                no_block_points.append((lon0, lat0))
                continue

            block_coverage[block_name] += 1

            if block_name not in block_params:
                failed_points += 1
                failed_reasons['no_block_params'] += 1
                continue

            myrange1 = block_params[block_name]['myrange1']
            myrange2 = block_params[block_name]['myrange2']

            lon_min_win = lon0 - window_size / 2
            lon_max_win = lon0 + window_size / 2
            lat_min_win = lat0 - window_size / 2
            lat_max_win = lat0 + window_size / 2

            lon_mask = (lon_grid >= lon_min_win) & (lon_grid <= lon_max_win)
            lat_mask = (lat_grid >= lat_min_win) & (lat_grid <= lat_max_win)

            if lon_mask.sum() < 2 or lat_mask.sum() < 2:
                failed_points += 1
                failed_reasons['insufficient_window'] += 1
                continue

            lon_idx1, lon_idx2 = np.where(lon_mask)[0][[0, -1]]
            lat_idx1, lat_idx2 = np.where(lat_mask)[0][[0, -1]]
            window_data = value_grid[lat_idx1:lat_idx2 + 1, lon_idx1:lon_idx2 + 1]

            if window_data.size < 4 or np.isnan(window_data).all():
                failed_points += 1
                failed_reasons['no_data'] += 1
                continue

            try:
                dt = detrend_2d(window_data)
            except Exception:
                failed_points += 1
                failed_reasons['detrend_failed'] += 1
                continue

            a2 = raas(dt, 5.5)
            if a2.size == 0:
                failed_points += 1
                failed_reasons['spectrum_failed'] += 1
                continue

            w2, pw2, ps2, ci1 = fractal_correction(a2, b2)
            low1, high1 = myrange1
            low2, high2 = myrange2

            if (high1 >= len(w2) or high2 >= len(w2) or
                    low1 >= len(w2) or low2 >= len(w2) or
                    low1 > high1 or low2 > high2):
                failed_points += 1
                failed_reasons['index_out_of_range'] += 1
                continue

            xs2 = w2[low1:high1 + 1]
            ys2 = pw2[low1:high1 + 1]
            xs = w2[low2:high2 + 1]
            ys = ps2[low2:high2 + 1]

            if len(xs2) < 2 or len(xs) < 2:
                failed_points += 1
                failed_reasons['insufficient_freq_bands'] += 1
                continue

            try:
                p2 = np.polyfit(xs2, ys2, 1)
                pl = np.polyfit(xs, ys, 1)
            except Exception:
                failed_points += 1
                failed_reasons['fit_failed'] += 1
                continue

            Zt = -p2[0] / (2 * np.pi)
            Z0 = -pl[0] / (2 * np.pi)
            Zb = 2 * Z0 - Zt

            if 10 <= Zb <= 60:  
                zb_full[i, j] = Zb
                valid_points += 1
            else:
                failed_points += 1
                failed_reasons['constraint_failed'] += 1


    if no_block_points:
        no_block_file = os.path.join(output_dir, "no_block_coverage_points.txt")
        with open(no_block_file, "w") as f:
            f.write("Lon Lat\n")
            for lon, lat in no_block_points:
                f.write(f"{lon:.3f} {lat:.3f}\n")


    if no_block_points:
        uncovered_plot_file = os.path.join(output_dir, "uncovered_points.png")
        plot_uncovered_points(no_block_points, block_polygons, uncovered_plot_file)


    return lons_full, lats_full, zb_full


def plot_uncovered_points(no_block_points, block_polygons, output_file):
    plt.figure(figsize=(12, 8))
    for block_name, polygon in block_polygons.items():
        poly_lons = [p[0] for p in polygon] + [polygon[0][0]]
        poly_lats = [p[1] for p in polygon] + [polygon[0][1]]
        plt.plot(poly_lons, poly_lats, 'k-', linewidth=1,
                 label=block_name if block_name == list(block_polygons.keys())[0] else "")

    if no_block_points:
        uncovered_lons = [p[0] for p in no_block_points]
        uncovered_lats = [p[1] for p in no_block_points]
        plt.scatter(uncovered_lons, uncovered_lats, c='red', s=20, marker='x', label='no point')

    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


def check_block_coverage(block_polygons, lon_range, lat_range, resolution=0.1):
    lon_min, lon_max = lon_range
    lat_min, lat_max = lat_range
    lons_check = np.arange(lon_min, lon_max + resolution, resolution)
    lats_check = np.arange(lat_min, lat_max + resolution, resolution)

    covered = 0
    not_covered = 0

    for lon in lons_check:
        for lat in lats_check:
            if find_block_for_point((lon, lat), block_polygons) is not None:
                covered += 1
            else:
                not_covered += 1

    coverage_rate = covered / (covered + not_covered) * 100
    print(f"Block coverage inspection: {covered}/{covered + not_covered} The point has been covered.({coverage_rate:.1f}%)")
    return coverage_rate

if __name__ == "__main__":
    b2_list = [1]
    ws_list = [3]
    data_file = "ordos.txt"
    measured_file = "heatflow_0.1.txt"  
    block_dir = "block"

    measured_lons, measured_lats, measured_heat_flows = load_measured_heat_flow(measured_file)
    block_polygons = load_block_polygons(block_dir)

    print(f"101-117°E, 32-42°N")
    print(f"find {len(block_polygons)} block")
    print(f" {len(measured_heat_flows)}heatflow")
    print(f"heatflowrange:{measured_heat_flows.min():.1f} - {measured_heat_flows.max():.1f} mW/m²")

    region_mask = (measured_lons >= 101) & (measured_lons <= 117) & (measured_lats >= 32) & (measured_lats <= 42)
    region_lons = measured_lons[region_mask]
    region_lats = measured_lats[region_mask]
    region_heat_flows = measured_heat_flows[region_mask]

    coverage_rate = check_block_coverage(block_polygons, (101, 117), (32, 42), resolution=0.1)

    successful_combinations = []

    for b2 in b2_list:
        for ws in ws_list:
            print(f"\n▶️ progree: b2={b2}, window_size={ws}")
            all_block_params = {}
            optimization_failed = False

            for block_name, block_polygon in block_polygons.items():
                print(f"  optimization {block_name}...")
                block_lons, block_lats, block_heat_flows = get_points_in_block(
                    block_polygon, region_lons, region_lats, region_heat_flows
                )


                if len(block_lons) == 0:
                
                    default_params = {
                        'myrange1': [3, 7],
                        'myrange2': [1, 4],
                        'rmse': 1000.0
                    }
                    all_block_params[block_name] = default_params
                    continue

                try:
                    best_params, best_loss = optimize_for_b2_ws(
                        b2, ws, data_file, block_lons, block_lats, block_heat_flows, n_trials=50
                    )

                    if best_params is None:
                        default_params = {
                            'myrange1': [3, 7],
                            'myrange2': [1, 4],
                            'rmse': 1000.0
                        }
                        all_block_params[block_name] = default_params
                    else:
                        block_best_params = {
                            'myrange1': [best_params['myrange1_start'], best_params['myrange1_end']],
                            'myrange2': [best_params['myrange2_start'], best_params['myrange2_end']],
                            'rmse': best_loss
                        }
                        all_block_params[block_name] = block_best_params
                        print(f"    ✅ : {best_loss:.4f} mW/m²")
                        print(f"        range1: {block_best_params['myrange1']}")
                        print(f"        range2: {block_best_params['myrange2']}")

                except Exception as e:
                    print(f"    ❌ defeat")
                    default_params = {
                        'myrange1': [3, 7],
                        'myrange2': [1, 4],
                        'rmse': 1000.0
                    }
                    all_block_params[block_name] = default_params

            main_result_dir = f"Curie_Depth_b2={b2}_ws={ws}"
            os.makedirs(main_result_dir, exist_ok=True)

            with open(os.path.join(main_result_dir, "block_parameters.txt"), "w", encoding="utf-8") as f:

            try:
                lons_full, lats_full, zb_full = compute_curie_depth_for_full_map(
                    b2, ws, data_file, all_block_params, block_polygons, main_result_dir, output_resolution=0.1
                )
                with open(os.path.join(main_result_dir, "curie_depth_map.txt"), "w") as f:
                    for i in range(len(lats_full)):
                        for j in range(len(lons_full)):
                            if not np.isnan(zb_full[i, j]):
                                f.write(f"{lons_full[j]:.3f} {lats_full[i]:.3f} {zb_full[i, j]:.3f}\n")

                print(f"✅ save curie_depth_map.txt")

                plt.figure(figsize=(12, 8))
                contour = plt.contourf(lons_full, lats_full, zb_full, levels=20, cmap='jet')
                plt.colorbar(contour, label='Curie Depth (km)')
                plt.xlabel('Longitude')
                plt.ylabel('Latitude')
                plt.title(f'Curie Depth Map (b2={b2}, ws={ws})')
                plt.savefig(os.path.join(main_result_dir, 'curie_depth_map.png'), dpi=300, bbox_inches='tight')
                plt.close()

                print(f"✅ ")
                successful_combinations.append((b2, ws, main_result_dir))

            except Exception as e:
                print(f"❌ ")

            print(f"✅ b2={b2}, ws={ws} save {main_result_dir}")

    for b2, ws, dir_path in successful_combinations:
        print(f"  - b2={b2}, ws={ws}: {dir_path}")
