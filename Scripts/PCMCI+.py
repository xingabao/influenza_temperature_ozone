# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from matplotlib import rcParams
from statsmodels.tsa.seasonal import seasonal_decompose
from statsmodels.tsa.stattools import adfuller

# Tigramite imports
from tigramite import data_processing as pp
from tigramite import plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests.parcorr import ParCorr

# ---------------------------------------------------------
# 1. Global Settings (Publication Quality)
# ---------------------------------------------------------
# Set random seed for reproducibility
np.random.seed(42)

FONT_FAMILY = 'Times New Roman'
FONT_SIZE = 14
rcParams['font.family'] = FONT_FAMILY
rcParams['font.size'] = FONT_SIZE
rcParams['axes.titlesize'] = FONT_SIZE + 2
rcParams['axes.labelsize'] = FONT_SIZE
rcParams['xtick.labelsize'] = FONT_SIZE
rcParams['ytick.labelsize'] = FONT_SIZE
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

# Mapping for LaTeX formatting
LATEX_LABELS = {
    'PM2.5': r'$PM_{2.5}$',
    'PM10': r'$PM_{10}$',
    'O3': r'$O_3$',
    'NO2': r'$NO_2$',
    'SO2': r'$SO_2$',
    'CO': r'$CO$'
}

# ---------------------------------------------------------
# 2. Advanced Data Preprocessing
# ---------------------------------------------------------
def check_and_fix_stationarity(data, var_names, missing_flag = 999.):
    """
    Performs ADF test. If non-stationary (p > 0.05), applies 1st differencing.
    Returns the corrected data.
    """
    print("  [Stationarity Check & Correction]")
    corrected_data = data.copy()
    
    for i, name in enumerate(var_names):
        series = corrected_data[:, i]
        # Create a clean series for ADF test (remove missing flags and NaNs)
        clean_mask = (series != missing_flag) & (~np.isnan(series))
        clean_series = series[clean_mask]
        
        if len(clean_series) < 20:
            print(f"    WARNING: {name} too short for ADF test.")
            continue
        
        try:
            p_value = adfuller(clean_series)[1]
            if p_value > 0.05:
                print(f"    -> {name} is Non-Stationary (p={p_value:.3f}). Applying 1st Differencing.")
                # Apply diff, prepend NaN to keep shape
                diffed = np.diff(series, prepend = np.nan)
                # Restore missing flag where it was originally missing (heuristic)
                # Note: Differencing introduces 1 NaN at start, which Tigramite handles.
                corrected_data[:, i] = diffed
            else:
                print(f"    -> {name} is Stationary (p={p_value:.3f}).")
        except Exception as e:
            print(f"    Error in ADF test for {name}: {e}")
    
    return corrected_data

def preprocess_data(df, target_cols, date_col = 'date', pandemic_col = 'pandemic'):
    """
    Standardizes data, handles seasonality (STL), masks pandemic, and enforces stationarity.
    """
    df = df.copy()
    df[date_col] = pd.to_datetime(df[date_col])
    df = df.set_index(date_col).sort_index()

    # 1. Handle Missing Values (Linear Interpolation)
    df[target_cols] = df[target_cols].interpolate(method = 'time', limit = 3)

    # Determine Periodicity
    freq_days = (df.index[1] - df.index[0]).days
    if freq_days >= 6:
        period = 52
    else:
        period = 365
    
    print(f"  [Preprocessing] Deseasonalizing (STL, period = {period}) and Masking Pandemic ...")
    
    processed_data = pd.DataFrame(index = df.index)
    
    # 2. STL Decomposition
    for col in target_cols:
        try:
            decomp = seasonal_decompose(df[col].dropna(), model = 'additive', period = period, extrapolate_trend = 'freq')
            processed_data[col] = decomp.resid
        except Exception as e:
            print(f"    Warning: STL failed for {col} ({e}). Using raw data.")
            processed_data[col] = df[col]

    # 3. Apply Pandemic Masking
    if pandemic_col in df.columns:
        is_pandemic = df[pandemic_col] == 1
        n_masked = is_pandemic.sum()
        if n_masked > 0:
            print(f"    Masking {n_masked} time points due to Pandemic.")
            processed_data.loc[is_pandemic, target_cols] = np.nan
    
    # 4. Convert to numpy and Fix Stationarity
    data_values = processed_data[target_cols].values
    final_values = check_and_fix_stationarity(data_values, target_cols)
    
    return final_values

# ---------------------------------------------------------
# 3. Visualization Functions
# ---------------------------------------------------------
def plot_heatmap(results, var_names, alpha, save_path):
    """
    Generates a heatmap showing the Maximum MCI strength across lags.
    Annotates with the specific Lag.
    """
    val_matrix = results['val_matrix']
    p_matrix = results['p_matrix']
    N = len(var_names)
    
    # Matrices to hold data for heatmap
    strength_matrix = np.zeros((N, N))
    annot_matrix = np.empty((N, N), dtype = object)
    
    for i in range(N): # Cause (Parents)
        for j in range(N): # Effect (Children)
            if i == j: 
                strength_matrix[i, j] = 0 # Ignore auto-correlation for clarity in this plot
                annot_matrix[i, j] = ""
                continue
            
            # Find the lag with the minimum p-value (most significant)
            # We skip lag 0 for cross-links if we assume no instantaneous effects, 
            # or keep it if PCMCI+ found them. Here we look at all lags.
            best_tau_idx = np.argmin(p_matrix[i, j, :]) 
            min_p = p_matrix[i, j, best_tau_idx]
            val = val_matrix[i, j, best_tau_idx]
            
            if min_p < alpha:
                strength_matrix[i, j] = val
                # Annotation: Strength \n (Lag)
                annot_matrix[i, j] = f"{val:.2f}\n({best_tau_idx})"
            else:
                strength_matrix[i, j] = 0
                annot_matrix[i, j] = ""

    plt.figure(figsize = (6, 5))
    ax = sns.heatmap(
        strength_matrix, 
        annot = annot_matrix,
        fmt = "", 
        xticklabels = var_names, 
        yticklabels = var_names,
        cmap = 'coolwarm', 
        vmin = -0.5,
        vmax = 0.5,
        center = 0,
        cbar = True,
        annot_kws = {'size': 18} 
    )
    
    # Control Left Y-Axis (l)
    plt.ylabel("Cause (Drivers)", fontsize = 18)
    plt.yticks(rotation = 0, fontsize = 18)

    # Control Bottom X-Axis (b)
    plt.xlabel("Effect (Targets)", fontsize = 18)
    ax.set_xticklabels(
        var_names, 
        rotation = 45, 
        ha = 'right', 
        rotation_mode = 'anchor',
        fontsize = 16
    )
    
    # Colorbar
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 16)
    cbar.set_label('MCI Strength', fontsize = 18)
    
    
    plt.xticks(rotation = 45, ha = 'right')
    plt.yticks(rotation = 0)
    
    if save_path:
        heatmap_path = save_path.replace('.pdf', '_heatmap.pdf')
        plt.savefig(heatmap_path, format = 'pdf', bbox_inches = 'tight', transparent = True)
        print(f"  [Saved] Heatmap: {heatmap_path}")
        
    plt.close()

def generate_standalone_legend(save_path):
    """
    
    Generates a separate PDF containing ONLY the colorbars.
    """
    print("  [Legend] Generating standalone legend file...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (17, 0.5))
    
    LABEL_SIZE = 20
    TICK_SIZE = 18
    
    # Edge Colorbar
    cmap_edges = plt.get_cmap('coolwarm')
    norm_edges = mpl.colors.Normalize(vmin = -0.5, vmax = 0.5)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap = cmap_edges, norm = norm_edges, orientation = 'horizontal', extend = 'both')
    cb1.set_label('Cross-MCI (Partial Correlation Strength)', fontsize = LABEL_SIZE)
    cb1.ax.tick_params(labelsize = TICK_SIZE) 
    
    # Node Colorbar
    cmap_nodes = plt.get_cmap('OrRd') 
    norm_nodes = mpl.colors.Normalize(vmin = 0, vmax = 1.0)
    cb2 = mpl.colorbar.ColorbarBase(ax2, cmap = cmap_nodes, norm = norm_nodes, orientation = 'horizontal')
    cb2.set_label('Auto-MCI (Auto-Correlation Strength)', fontsize = LABEL_SIZE)
    cb2.ax.tick_params(labelsize = TICK_SIZE) 
    
    plt.subplots_adjust(hspace = 1.5)
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok = True)
        plt.savefig(save_path, format = 'pdf', bbox_inches = 'tight', transparent = True)
    
    plt.close()

# ---------------------------------------------------------
# 4. Causal Discovery Logic
# ---------------------------------------------------------
def get_domain_knowledge_links(var_names, tau_min, tau_max, flu_label):
    """
    
    Defines 'Hard' constraints based on biological plausibility.
    """
    selected_links = {}
    var_dict = {var: i for i, var in enumerate(var_names)}
    
    temp_idx = var_dict.get('T')
    humid_idx = var_dict.get('AH')
    
    for idx, var in enumerate(var_names):
        parents = []
        # CASE 1: Flu (The Outcome)
        if var == flu_label:
            parents = [
                (other_idx, -tau) 
                for other_idx, other_var in enumerate(var_names) 
                for tau in range(tau_min, tau_max + 1)
            ]
        # CASE 2: Meteorology (Exogenous)
        elif var in ['T', 'AH']:
            met_indices = [v for k, v in var_dict.items() if k in ['T', 'AH']]
            parents = [(p_idx, -tau) for p_idx in met_indices for tau in range(tau_min, tau_max + 1)]
        # CASE 3: Pollutants (Intermediate)
        else: 
            parents = [(idx, -tau) for tau in range(tau_min, tau_max + 1)]
            if temp_idx is not None: 
                parents += [(temp_idx, -tau) for tau in range(max(0, tau_min), tau_max + 1)]
            if humid_idx is not None: 
                parents += [(humid_idx, -tau) for tau in range(max(0, tau_min), tau_max + 1)]
                
        selected_links[idx] = parents
    return selected_links

def save_results_to_csv(results, var_names, alpha, save_path):
    """
    
    Exports significant links to CSV for 'Source Data'.
    """
    val_matrix = results['val_matrix']
    p_matrix = results['p_matrix']
    rows = []
    
    for j, target in enumerate(var_names):
        for i, source in enumerate(var_names):
            if i == j: continue
            for tau in range(val_matrix.shape[2]):
                p_val = p_matrix[i, j, tau]
                if p_val < alpha:
                    rows.append({
                        'Source': source,
                        'Target': target,
                        'Lag': tau,
                        'MCI_Value': val_matrix[i, j, tau],
                        'P_Value': p_val
                    })
    
    if rows:
        df_res = pd.DataFrame(rows)
        csv_path = save_path.replace('.pdf', '_SourceData.csv')
        df_res.to_csv(csv_path, index = False)
        print(f"  [Saved] Source Data: {csv_path}")

def run_pcmci_analysis(data_values, var_names, flu_label, tau_max, pc_alpha, save_path):
    """
    
    Runs PCMCI+ and generates Graph + Heatmap + CSV.
    """
    missing_flag = 999.
    data_values_filled = np.nan_to_num(data_values, nan = missing_flag)
    
    dataframe = pp.DataFrame(data_values_filled, var_names = var_names, missing_flag = missing_flag)
    
    parcorr = ParCorr(significance = 'analytic')
    pcmci = PCMCI(dataframe = dataframe, cond_ind_test = parcorr, verbosity = 0)
    
    selected_links = get_domain_knowledge_links(var_names, tau_min = 0, tau_max = tau_max, flu_label = flu_label)
    
    print(f"  [PCMCI+] Running Causal Discovery (alpha = {pc_alpha})...")
    results = pcmci.run_pcmciplus(tau_min = 0, tau_max = tau_max, pc_alpha = pc_alpha, selected_links = selected_links)
    
    # 1. Network Plot
    plt.figure(figsize = (5, 5)) 
    tp.plot_graph(
        val_matrix = results['val_matrix'],
        graph = results['graph'],
        var_names = var_names,
        show_colorbar = False,
        vmin_edges = -0.5,
        vmax_edges = 0.5,
        cmap_edges = 'coolwarm',
        node_label_size = 22,
        link_label_fontsize = 18,
        arrow_linewidth = 5.0,
        node_size = 0.5,
        node_aspect = 1.4,
        network_lower_bound = 0.1,
    )
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok = True)
        plt.savefig(save_path, format = 'pdf', bbox_inches = 'tight', pad_inches = 0.01, transparent = True)
        print(f"  [Saved] Network Graph: {save_path}")
    plt.close()
    
    # 2. Heatmap Plot (New)
    plot_heatmap(results, var_names, pc_alpha, save_path)
    
    # 3. Export Data (New)
    save_results_to_csv(results, var_names, pc_alpha, save_path)
    
    return results

# ---------------------------------------------------------
# 5. Main Execution Loop
# ---------------------------------------------------------

if __name__ == '__main__':
    
    # Configuration
    # Use raw string for Windows paths to avoid escape character issues
    from pathlib import Path
    BASE_DIR = os.path.dirname(os.path.dirname(Path(__file__).resolve()))
    DATA_DIR = os.path.join(BASE_DIR, 'data')
    OUT_DIR = os.path.join(BASE_DIR, 'Figures/PCMCI+')
    
    MAX_LAG = 4
    
    cities = ['Macau', 'HK', 'Zhuhai']
    pollutants = ['PM10', 'PM2.5', 'SO2', 'NO2', 'O3', 'CO']
    flus = ['FLUA', 'FLUB']
    
    print("Starting Analysis (Publication Ready Mode)...")

    # 1. Generate Shared Legend
    legend_path = os.path.join(OUT_DIR, "Shared_Legend.pdf")
    generate_standalone_legend(legend_path)

    # 2. Run Analysis Loop
    for city in cities:
        
        PC_ALPHA = 0.05

        
        for flu in flus:
            file_path = os.path.join(DATA_DIR, city, 'FLU-CL-AQ.xlsx')
            
            if not os.path.exists(file_path):
                print(f"CRITICAL ERROR: Data file not found: {file_path}")
                continue
                
            raw_df = pd.read_excel(file_path)
            
            # Check for pandemic column
            if 'pandemic' not in raw_df.columns:
                print("WARNING: 'pandemic' column not found. Proceeding without masking.")
            
            for pol in pollutants:
                print(f"\n>>> Processing: {city} | {flu} | {pol}")
                
                # Dynamic column selection
                pol_col = pol
                if pol not in raw_df.columns and f'IAQI_{pol}' in raw_df.columns:
                    pol_col = f'IAQI_{pol}'
                
                required_cols = ['date', flu, 'ave.temp', 'humidity', pol_col]
                if 'pandemic' in raw_df.columns:
                    required_cols.append('pandemic')
                
                if not all(c in raw_df.columns for c in required_cols):
                    print(f"  Skipping: Missing columns in {city}.")
                    continue
                
                # Extract subset
                df_subset = raw_df[required_cols].copy()
                
                # Rename for internal processing
                rename_map = {
                    flu: 'Flu', 
                    'ave.temp': 'T',
                    'humidity': 'AH',
                    pol_col: pol
                }
                df_subset.rename(columns = rename_map, inplace = True)
                
                internal_vars = ['T', 'AH', pol, 'Flu']
                
                # LaTeX Label for Plotting
                flu_tex = r'$Flu_A$' if flu == 'FLUA' else r'$Flu_B$'
                display_names = ['T', 'AH', LATEX_LABELS.get(pol, pol), flu_tex]
                
                # Preprocess (Includes Stationarity Fix)
                data_values = preprocess_data(
                    df_subset, 
                    target_cols = internal_vars, 
                    date_col = 'date', 
                    pandemic_col = 'pandemic' if 'pandemic' in df_subset else None
                )
                
                # Run PCMCI
                out_name = f"{city}_{flu}_{pol}.pdf"
                save_path = os.path.join(OUT_DIR, out_name)
                
                try:
                    run_pcmci_analysis(
                        data_values = data_values,
                        var_names = display_names,
                        flu_label = flu_tex, 
                        tau_max = MAX_LAG,
                        pc_alpha = PC_ALPHA,
                        save_path = save_path
                    )
                except Exception as e:
                    print(f"  ERROR: {e}")
                    import traceback
                    traceback.print_exc()

    print("\nAnalysis Complete.")
    print(f"Outputs saved to: {OUT_DIR}")
    print("Includes: Network Graphs (.pdf), Heatmaps (_heatmap.pdf), and Source Data (_SourceData.csv)")