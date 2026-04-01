# -*- coding: utf-8 -*-
"""
Requires:
- Python packages: pandas, numpy, matplotlib, semopy, plspm, scikit-learn, graphviz (pip install graphviz)
- Graphviz executable installed and on PATH
"""

import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)
import matplotlib.pyplot as plt
from plspm.plspm import Plspm
from plspm.config import Config, Structure, MV
from plspm.mode import Mode
from plspm.scheme import Scheme
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24
plt.rcParams['axes.unicode_minus'] = False

def semplot_plspm(
        plspm_model, 
        filename: str,
        comdata = None, 
        plot_covs = False, 
        plot_exos = True, 
        images = None,
        engine = 'dot', 
        latshape = 'circle', 
        plot_ests = True, 
        show = False, 
        title = 'PLS-PM',
        splines = 'true',
        rankdir = 'TB',
        dpi = 300,
        graph_size = None,
        base_family = 'Times-Roman',
        base_size = 24
    ):
    try:
        import graphviz as gv
    except ImportError:
        raise ModuleNotFoundError("Module `graphviz` not found. Please install system Graphviz first, then pip install graphviz.")
        
    if images is None:
        images = dict()
        
    t = filename.split('.')
    filename_base, ext = '.'.join(t[:-1]), t[-1] if len(t) > 1 else 'png'
    
    try:
        gof_value = plspm_model.goodness_of_fit()
        title_with_gof = f"GOF = {gof_value:.4f}" if title == None else f"{title}\nGOF = {gof_value:.4f}"
    except Exception as e:
        print(f"Error retrieving `GOF` value: {e}")
        title_with_gof = title
    
    g = gv.Digraph('G', format = ext, engine = engine, graph_attr = {'label': title_with_gof, 'labelloc': 'b', 'fontsize': str(base_size), 'fontname': base_family, 'dpi': str(dpi)})
    if graph_size: g.graph_attr['size'] = graph_size; g.graph_attr['ratio'] = 'fill'
    g.attr(overlap = 'scale', splines = splines, rankdir = rankdir)
    g.attr('edge', fontsize = str(base_size), fontname = base_family)

    try:
        inner_model = plspm_model.inner_model()
        latent_vars = inner_model['to'].unique().tolist()
        exogenous_vars = set(inner_model['from'].unique()) - set(inner_model['to'].unique())
        latent_vars.extend(list(exogenous_vars))
        outer_model = plspm_model.outer_model()
        observed_vars = outer_model.index.tolist()
    except Exception as e:
        print(f"Error retrieving model information: {e}")
        return None
    
    g.attr('node', shape = latshape, fillcolor = '#e8f4fd', style = 'filled', fontsize = str(base_size), fontname = base_family, color = '#2e86ab', penwidth = '2')
    for lat in latent_vars:
        if lat in images:
            g.node(lat, label = '', image = images[lat])
        else:
            g.node(lat, label = lat)
    
    g.attr('node', shape = 'box', fillcolor = '#fff2cc', style = 'filled', fontsize =  str(base_size), fontname = base_family, color = '#d6b656', penwidth = '2')
    for obs in observed_vars:
        if obs in images:
            g.node(obs, label='', image=images[obs])
        else:
            g.node(obs, label=obs)
    
    for _, row in inner_model.iterrows():
        from_var, to_var = row['from'], row['to']
        estimate = row.get('estimate', np.nan)
        p_value = row.get('p>|t|', np.nan)
        if plot_ests and pd.notna(estimate):
            label = f'<{estimate:.4f}'
            if pd.notna(p_value):
                if p_value <= 0.001: label += "***"
                elif p_value <= 0.01: label += "**"
                elif p_value <= 0.05: label += "*"
                elif p_value <= 0.1: label += "†"
                label += f'<br/><i>P</i>={p_value:.4f}>'
        else:
            label = ''
        if pd.notna(p_value):
            if p_value <= 0.05:
                edge_color, edge_style, penwidth = '#d32f2f', 'solid', '3'
            elif p_value <= 0.1:
                edge_color, edge_style, penwidth = '#ff9800', 'solid', '2'
            else:
                edge_color, edge_style, penwidth = '#757575', 'dashed', '1'
        else:
            edge_color, edge_style, penwidth = '#757575', 'solid', '1'
        g.edge(from_var, to_var, label = label, color = edge_color, style = edge_style, penwidth = penwidth, fontcolor = 'black')
    
    try:
        if hasattr(plspm_model, 'config') and hasattr(plspm_model.config, 'latent_variables'):
            comdata = {}
            for lv_name, lv_config in plspm_model.config.latent_variables.items():
                if hasattr(lv_config, 'measurement_variables'):
                    comdata[lv_name] = [mv.name for mv in lv_config.measurement_variables]
            # print("Extracting variable relationships from model config:", comdata)
        for lv_name, mv_list in comdata.items():
            for mv_name in mv_list:
                if mv_name in observed_vars:
                    try:
                        loading = outer_model.loc[mv_name, 'loading']
                        label = f'{loading:.4f}' if plot_ests else ''
                    except:
                        label = ''
                    g.edge(mv_name, lv_name, label = label, color = '#4caf50', style = 'solid', penwidth = '2', fontcolor = '#2e7d32')
    except Exception as e:
        print(f"Error plotting loadings: {e}")
    
    try:
        g.render(filename_base, view = show, cleanup = True)
        print(f"Model diagram saved as: {filename_base}.{ext}")
        try:
            print(f"GOF value: {gof_value:.4f}")
        except:
            pass
    except Exception as e:
        print(f"Error rendering graph: {e}")
    return g

def create_effects_plot(df, title = 'Influenza', save_path = None, base_size = 24, base_family = 'Times New Roman'):
    
    '''
    
    ax.set_title(f'Effects on {title}: Direct, Indirect, and Total', fontsize = base_size, fontfamily = base_family, fontweight = 'normal', pad = 20)
    '''
    pivot_df = df.pivot_table(index = 'predictor', columns = 'effect_type', values = 'effect_size', fill_value = 0)
    effect_columns = ['direct_effect', 'indirect_effect', 'total_effect']
    for col in effect_columns:
        if col not in pivot_df.columns:
            pivot_df[col] = 0
    pivot_df = pivot_df[effect_columns]
    categories = pivot_df.index.tolist()
    x = np.arange(len(categories))
    direct = pivot_df['direct_effect'].values
    indirect = pivot_df['indirect_effect'].values
    total = pivot_df['total_effect'].values
    width = 0.25
    fig, ax = plt.subplots(figsize = (10, 10))
    color_direct = '#2E8B8B'
    color_indirect = '#90EE90'
    color_total = '#FFD700'
    ax.bar(x - width, direct, width, label = 'Direct Effect', color = color_direct)
    ax.bar(x, indirect, width, label = 'Indirect Effect', color = color_indirect)
    ax.bar(x + width, total, width, label = 'Total Effect', color = color_total)
    ax.set_ylabel('Effect Size (Standardized)', fontweight = 'normal', fontsize = base_size, fontfamily = base_family)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize = base_size, fontfamily = base_family, fontweight = 'normal', rotation = 90)
    ax.set_ylim(-1.0, 1.0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False) 
    y_ticks = np.arange(-1.0, 1.1, 0.5)
    ax.set_yticks(y_ticks)
    ax.grid(False)
    ax.set_axisbelow(True)
    ax.axhline(0, color = 'black', linewidth = 1.5)
    ax.tick_params(axis = 'both', which = 'major', labelsize = base_size, width = 1.5, length = 6, direction = 'out')
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    legend = ax.legend(loc = 'center', bbox_to_anchor = (0.7, 0.2), frameon = True, fontsize = base_size, fancybox = True, shadow = True)
    for text in legend.get_texts():
        text.set_fontweight('normal')
        text.set_fontfamily(base_family)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_edgecolor('gray')
    legend.get_frame().set_linewidth(1)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi = 300, bbox_inches = 'tight', transparent = True)
        print(f"Figure saved to: {save_path}")
    plt.show()
    return fig, ax

def prepare_effects_data_for_plot(effects, target = 'Influenza'):
    target_effects = effects[effects['to'] == target].copy()
    if target_effects.empty:
        print(f"No effect data found for {target}")
        return pd.DataFrame()
    plot_data = []
    for _, row in target_effects.iterrows():
        predictor = row['from']
        plot_data.append({'predictor': predictor, 'outcome': target, 'effect_type': 'direct_effect', 'effect_size': row['direct']})
        plot_data.append({'predictor': predictor, 'outcome': target, 'effect_type': 'indirect_effect', 'effect_size': row['indirect']})
        plot_data.append({'predictor': predictor, 'outcome': target, 'effect_type': 'total_effect', 'effect_size': row['total']})
    return pd.DataFrame(plot_data)

def extract_path_coefficients(inner_model, inner_summary):
    path_coef_list = []
    for target in inner_model['to'].unique():
        target_data = inner_model[inner_model['to'] == target]
        for _, row in target_data.iterrows():
            path_coef_list.append({
                'Predictor': row['from'],
                'Response': row['to'],
                'Estimate': row.get('estimate', np.nan),
                'StdError': row.get('std error', np.nan),
                'tValue': row.get('t', np.nan),
                'pValue': row.get('p>|t|', np.nan)
            })
    path_coefficients = pd.DataFrame(path_coef_list)
    r_squared_data = inner_summary[inner_summary['type'] == 'Endogenous'][['r_squared']].copy()
    r_squared_data['Response'] = r_squared_data.index
    r_squared_data = r_squared_data.rename(columns={'r_squared': 'R2'})
    path_coefficients = path_coefficients.merge(r_squared_data, on='Response', how='left')
    def get_significance(p_value):
        if pd.isna(p_value): return "NA"
        if p_value < 0.001: return "***"
        elif p_value < 0.01: return "**"
        elif p_value < 0.05: return "*"
        else: return "ns"
    path_coefficients['Significance'] = path_coefficients['pValue'].apply(get_significance)
    path_coefficients = path_coefficients.sort_values(['Response', 'Predictor'])
    return path_coefficients

def print_target_effects(effects, target = 'Influenza'):
    target_effects = effects[effects['to'] == target]
    if target_effects.empty:
        return
    print(f"\nEffects on {target}:")
    for _, row in target_effects.iterrows():
        print(f"{row['from']}: Direct={row['direct']:.4f}, Indirect={row['indirect']:.4f}, Total={row['total']:.4f}")

def add_seasonality_columns(df, date_col = 'date', period = 52, prefix = 'season', data_freq = 'W'):
    if date_col not in df.columns:
        raise ValueError(f"Date column not found: {date_col}")
    df = df.sort_values(date_col).copy()
    days_diff = (df[date_col] - df[date_col].min()).dt.days.values.astype(float)
    if data_freq == 'W':
        t = days_diff / 7.0
    else:
        t = days_diff
    sin_col = f'{prefix} sin'
    df[sin_col] = np.sin(2 * np.pi * t / period)
    return df, [sin_col]

def add_lag_columns(df, cols, lag = 4, group_by = None):
    df = df.copy()
    new_cols = []
    if group_by is None:
        for c in cols:
            new_c = f'{c}_lag{lag}'
            df[new_c] = df[c].shift(lag)
            new_cols.append(new_c)
    else:
        df = df.sort_values(group_by + ['date']).copy()
        for c in cols:
            new_c = f'{c}_lag{lag}'
            df[new_c] = df.groupby(group_by)[c].shift(lag)
            new_cols.append(new_c)
    return df, new_cols

def run_plspm_analysis_optimized(data, comdata, model_formulas, figdir, target, show_plot = True, figname1 = 'PLS-PM.EffectsPlot.png', figname2 = 'PLS-PM.ModelPlot.png'):
    try:
        print(f"Data loaded successfully, shape: {data.shape}")
        
        structure = Structure()
        for target_var, source_vars in model_formulas.items():
            structure.add_path(source_vars, [target_var])
        path_matrix = structure.path()
        
        config = Config(path_matrix, scaled = True)
        for lv_name, mv_list in comdata.items():
            existing_mvs = [mv for mv in mv_list if mv in data.columns]
            if not existing_mvs:
                print(f"Warning: All measurement variables for latent variable {lv_name} are missing from data. Skipped.")
                continue
            mv_objects = [MV(mv_name) for mv_name in existing_mvs]
            config.add_lv(lv_name, Mode.A, *mv_objects)
        
        print("\nRunning PLS-PM model...")
        plspm_model = Plspm(
            data.dropna(), 
            config,
            scheme = Scheme.CENTROID,
            iterations = 1000,
            tolerance = 1e-6,
            bootstrap = True, 
            bootstrap_iterations = 1000
        )
        
        print("\n=== PLS-PM Model Results Summary ===")
        inner_model = plspm_model.inner_model()
        print("\nInner Model (Path Coefficients):")
        print(inner_model)
        
        outer_model = plspm_model.outer_model()
        print("\nOuter Model (Loadings):")
        print(outer_model)
        
        inner_summary = plspm_model.inner_summary()
        print("\nInner Model Summary (R-squared):")
        print(inner_summary)
        
        gof = plspm_model.goodness_of_fit()
        print(f"\nGoodness of Fit (GoF): {gof:.6f}")
        
        effects = plspm_model.effects()
        effects.to_csv(f'{wkdir}/Table_S_PLSPM_Effects_Summary_{place}.csv', index = False)
        print("\nEffects Analysis (Direct, Indirect, Total):")
        print(effects)
        
        path_coefficients = extract_path_coefficients(inner_model, inner_summary)
        path_coefficients.to_csv(f'{wkdir}/Table_S_PLSPM_Path_Coefficients_{place}.csv', index = False)
        
        print("\nCleaned Path Coefficients Table:")
        print(path_coefficients)
        
        scores = plspm_model.scores()
        correlation_matrix = scores.corr()
        correlation_matrix.to_csv(f'{wkdir}/Table_S_Latent_Correlations_{place}.csv', index = False)
        print("\nLatent Variable Correlation Matrix:")
        print(correlation_matrix.round(3))
        
        plot_df = prepare_effects_data_for_plot(effects, target)
        if not plot_df.empty:
            print("\n=== Effects Data Summary ===")
            print(plot_df.pivot_table(index = 'predictor', columns = 'effect_type', values = 'effect_size', fill_value = 0))
            fig, ax = create_effects_plot(plot_df, save_path = f'{figdir}/{figname1}')
        
        print_target_effects(effects)
        
        semplot_plspm(plspm_model, filename = f'{figdir}/{figname2}', comdata = comdata, latshape = 'ellipse', plot_ests = True, show = show_plot, title = None, rankdir = 'TB')
        
        return plspm_model, effects, path_coefficients
    
    except Exception as e:
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()
        return None, None, None

# =============================================================================
# New Functions: Sensitivity Analysis and Plotting
# =============================================================================
def run_sensitivity_analysis(data_clean, comdata_template, model_formulas, lag_cols, lag_range = range(0, 9), target_latent = 'Influenza'):
    """
    Iterate through different lag periods to evaluate model performance.
    """
    results = []
    print(f"\n>>> Starting Sensitivity Analysis (Lag {min(lag_range)} - {max(lag_range)})...")
    
    for lag in lag_range:
        try:
            # 1. Use a copy of the clean data for each loop
            df_curr = data_clean.copy()
            
            # 2. If lag > 0, generate lag columns and replace original columns
            if lag > 0:
                df_curr, _ = add_lag_columns(df_curr, lag_cols, lag=lag)
                
                # Simulate column replacement logic: point original variable names to lagged data columns
                for col in lag_cols:
                    lag_col_name = f"{col}_lag{lag}"
                    if lag_col_name in df_curr.columns:
                        df_curr[col] = df_curr[lag_col_name]
                        # Note: No need to delete here, as we are just overwriting values or letting comdata point to original names
            
            # 3. Build model config (simplified, just to get GoF and R2)
            structure = Structure()
            for target, sources in model_formulas.items():
                structure.add_path(sources, [target])
            
            config = Config(structure.path(), scaled=True)
            for lv, mvs in comdata_template.items():
                # Ensure only existing columns are added
                valid_mvs = [m for m in mvs if m in df_curr.columns]
                if valid_mvs:
                    config.add_lv(lv, Mode.A, *[MV(m) for m in valid_mvs])
            
            # 4. Run PLS-PM (Silent run to reduce output)
            plspm = Plspm(df_curr.dropna(), config, scheme=Scheme.CENTROID, iterations = 300)
            
            gof = plspm.goodness_of_fit()
            r2 = plspm.inner_summary().loc[target_latent, 'r_squared']
            
            results.append({'Lag': lag, 'GoF': gof, 'R2': r2})
            print(f"   Lag {lag}: GoF={gof:.4f}, R2={r2:.4f}")
            
        except Exception as e:
            print(f"   Lag {lag} Failed: {str(e)}")
            
    return pd.DataFrame(results)

def plot_sensitivity_results(df_results, save_path):
    """
    
    Plot sensitivity analysis results
    'Sensitivity Analysis: Model Performance vs. Lag Weeks'
    """
    fig, ax1 = plt.subplots(figsize = (14, 7))
    
    color = 'tab:blue'
    ax1.set_xlabel('Lag Periods (Weeks)', fontsize = 24, fontfamily = 'Times New Roman')
    ax1.set_ylabel('Goodness of Fit (GoF)', color = color, fontsize = 24, fontfamily = 'Times New Roman')
    ax1.plot(df_results['Lag'], df_results['GoF'], color = color, marker = 'o', linewidth = 2, label = 'GoF')
    ax1.tick_params(axis = 'y', labelcolor = color)
    ax1.grid(True, alpha = 0.3)
    
    ax2 = ax1.twinx()  
    color = 'tab:red'
    ax2.set_ylabel('R² (Influenza)', color=color, fontsize = 24, fontfamily = 'Times New Roman')
    ax2.plot(df_results['Lag'], df_results['R2'], color = color, marker = 's', linestyle = '--', linewidth = 2, label = 'R²')
    ax2.tick_params(axis = 'y', labelcolor = color)

    plt.tight_layout()
    plt.savefig(save_path, dpi = 300, transparent = True, bbox_inches = 'tight')
    print(f"Sensitivity analysis plot saved to: {save_path}")

if __name__ == '__main__':

    # Data Path
    import os
    from pathlib import Path
    BASE_DIR = os.path.dirname(os.path.dirname(Path(__file__).resolve()))
    place = 'Macau'   # Options: 'HK' or 'Macau'
    wkdir = f'{BASE_DIR}/Tables'
    file_path = f'{BASE_DIR}/data/{place}/FLU-CL-AQ.xlsx'
    fig_dir = f'{BASE_DIR}/Figures'
    date_col = 'date'
    
    period = 52         
    
    comdata_original = {
        'Atmosphere': ['pressure', 'wind.speed'],
        'Temperature': ['max.temp', 'min.temp', 'ave.temp'],
        'Moisture': ['humidity', 'precipitation'],
        'Particulates': ['PM10', 'PM2.5'],
        'Primary Pollutants': ['SO2', 'NO2', 'CO'],
        'Ozone': ['O3'],
        'FLU': ['FLUA', 'FLUB'], 
        'Policy': ['pandemic']
    }

    # 2. Structural Model
    # Define causal paths between latent variables (Left is affected by Right)
    model_formulas_original = {
        
        "Particulates": ["Atmosphere", "Moisture", "Temperature"],
        "Primary Pollutants": ["Atmosphere", "Temperature", 'Moisture'],
        "Ozone": ["Atmosphere", "Temperature", 'Primary Pollutants'],
        "FLU": [
            "Temperature", 
            "Moisture", 
            "Particulates", 
            "Primary Pollutants", 
            "Ozone", 
            "Policy"
        ]
    }
    
    # Optional: Rename Mapping
    latent_rename_map = {
        'FLU': 'Influenza'
    }
    observed_rename_map = {
        'Pandemic': 'pandemic',
        'precipitation': 'Precipitation',
        'humidity': 'Humidity',
        'pressure': 'Pressure',
        'wind.speed': 'Wind Speed',
        'ave.temp': 'Average Temperature',
        'max.temp': 'Maximum Temperature',
        'min.temp': 'Minimus Temperature',
        'PM10': '<PM<sub>10</sub>>',
        'PM2.5': '<PM<sub>2.5</sub>>',
        'O3': '<O<sub>3</sub>>',
        'SO2': '<SO<sub>2</sub>>',
        'NO2': '<NO<sub>2</sub>>',
        'FLUA': 'Influenza A',
        'FLUB': 'Influenza B'
    }
    
    # Read Data and Preprocess
    data = pd.read_excel(file_path)
    
    # Smooth weekly influenza positivity rates using a 3-week trailing moving average
    # to reduce short-term stochastic fluctuation without using future observations.
    data['FLUA'] = data['FLUA'].rolling(window = 3, min_periods = 1, center = False).mean()
    data['FLUB'] = data['FLUB'].rolling(window = 3, min_periods = 1, center = False).mean()
    
    # Use IAQI
    data['PM10'] = data['IAQI_PM10']
    data['PM2.5'] = data['IAQI_PM2.5']
    data['SO2'] = data['IAQI_SO2']
    data['NO2'] = data['IAQI_NO2']
    data['O3'] = data['IAQI_O3']
    data['CO'] = data['IAQI_CO']

    if date_col not in data.columns: raise ValueError(f"Date column {date_col} not found in Excel")
    data[date_col] = pd.to_datetime(data[date_col])
    data = data[(data[date_col] > '2010-12-31') & (data[date_col] < '2025-01-01')].copy()
    
    # Pandemic Dummy Variable
    pandemic_start = pd.Timestamp('2020-01-20')
    pandemic_end = pd.Timestamp('2022-12-31')
    
    data['pandemic'] = data[date_col].apply(
        lambda x: 1 if pandemic_start <= x <= pandemic_end else 0
    )
    if 'holiday' not in data.columns:
        print("Warning: 'holiday' column not found, automatically filled with 0. Please check data source.")
        data['holiday'] = 0

    # Apply Observed Variable Renaming
    data.rename(columns = observed_rename_map, inplace = True)
    
    # Add Seasonality
    data, season_cols = add_seasonality_columns(
        data, 
        date_col = date_col, 
        period = period, 
        prefix = 'season', 
        data_freq = 'W' 
    )
    
    # Prepare Renamed comdata and model_formulas (For Analysis)
    comdata_renamed = {}
    for old_latent, old_obs_list in comdata_original.items():
        new_latent = latent_rename_map.get(old_latent, old_latent)
        new_obs_list = [observed_rename_map.get(obs, obs) for obs in old_obs_list]
        comdata_renamed[new_latent] = list(new_obs_list)
        
    model_formulas_renamed = {}
    for old_dep, old_srcs in model_formulas_original.items():
        new_dep = latent_rename_map.get(old_dep, old_dep)
        new_srcs = [latent_rename_map.get(s, s) for s in old_srcs]
        model_formulas_renamed[new_dep] = new_srcs

    # Define Columns to Lag (Use Renamed Column Names)
    lag_base_cols = [
        'Pressure', 'Wind Speed',
        'Maximum Temperature', 'Minimus Temperature',
        'Humidity', 'Precipitation',
        '<PM<sub>10</sub>>', '<PM<sub>2.5</sub>>', '<SO<sub>2</sub>>', '<NO<sub>2</sub>>', '<O<sub>3</sub>>', 'CO',
        'Pandemic'
    ]
    lag_base_cols = [c for c in lag_base_cols if c in data.columns]

    # =========================================================================
    # Step 1: Run Sensitivity Analysis (New Part)
    # =========================================================================
    # We do not modify data here, but pass a copy of data for testing
    sensitivity_df = run_sensitivity_analysis(
        data_clean = data,
        comdata_template = comdata_renamed,
        model_formulas = model_formulas_renamed,
        lag_cols = lag_base_cols,
        lag_range = range(0, 9),
        target_latent = 'Influenza'
    )
    
    # Save Sensitivity Analysis Plot
    plot_sensitivity_results(sensitivity_df, save_path = f'{fig_dir}/Sensitivity_Analysis_{place}.pdf')
    
    # Automatically Select Best Lag (Based on R2)
    best_lag = sensitivity_df.loc[sensitivity_df['R2'].idxmax(), 'Lag']
    best_r2 = sensitivity_df['R2'].max()
    print(f"\n>>> Sensitivity analysis complete. Best lag period: {best_lag} weeks (Selected based on R2 = {best_r2:.4f}).")

    # You can choose to force use 4, or use the analyzed best_lag
    # lag_periods = 4  # Uncomment this line if you want to force use 4
    lag_periods = int(best_lag) # Default to using the analyzed best value
    print(f">>> Final model will use Lag = {lag_periods}")

    # =========================================================================
    # Step 2: Run Final Model with Selected Lag (Original Logic)
    # =========================================================================
    
    # Call Lag Function
    if lag_periods > 0:
        data, lagged_cols = add_lag_columns(data, lag_base_cols, lag = lag_periods, group_by = None)
        
        # Replace Data Columns: Replace values of original column names with lagged column data
        # This way variable names in comdata_renamed don't need to change, but data is already lagged
        col_mv_rename = {}
        for latent, mv_list in list(comdata_renamed.items()):
            extended = list(mv_list)
            for mv in mv_list:
                lag_mv = f"{mv}_lag{lag_periods}"
                if lag_mv in data.columns:
                    # Record columns to replace
                    col_mv_rename[mv] = lag_mv
                    # Ensure mv is in list (Actually no need to change, as we are replacing data values)
            comdata_renamed[latent] = sorted(set(extended))
            
        # Execute Data Replacement
        for old_col, new_col in col_mv_rename.items():
            data[old_col] = data[new_col]
            del data[new_col] 
    
    # Run PLS-PM
    model, effects, path_coef = run_plspm_analysis_optimized(
        data = data,
        comdata = comdata_renamed,
        model_formulas = model_formulas_renamed,
        target = 'Influenza',
        show_plot = False,
        figdir = fig_dir,
        figname1 = f'PLS-PM_{place}_effects.pdf',
        figname2 = f'PLS-PM_{place}_plspm.pdf'
    )
    
    if effects is not None and False: effects.to_csv(f'{wkdir}/effects.csv', index = False)
    if path_coef is not None and False: path_coef.to_csv(f'{wkdir}/path_coef.csv', index = False)

    print("\nExecution complete.")