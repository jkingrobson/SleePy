 #!/usr/bin/env python
# coding: utf-8



"""
SleePy v2.2
A python program for the analysis of sleep and circadian metrics with automated management of missing data
Author: Josh King-Robson (j.king-robson@ucl.ac.uk)
"""


import pandas as pd 
import numpy as np
import statistics
import glob
import os
import csv
import datetime, time
import math
from math import isnan
import scipy
from scipy import stats
from collections import defaultdict




def circadian_metrics(raw_df, csv_path,
                      l5_imputation='mean', l5_fallback=None,
                      m10_imputation='mean', m10_fallback=None,
                      is_imputation=None, is_fallback=None,
                      iv_imputation='median', iv_fallback=None,
                      m16l8_time_imputation='mean', m16l8_time_fallback=None,
                     ):
    
    #Calculator for circadian metrics, includes imputation for missing data as per above defaults. 

                    
    # Create datetime; IF using non actiware (?UK_) format, please examine output carefully to check this works. 
    # The Actiware date input is inconsistent; this manages that, and is hopefully robust enough to manage other data. 
    raw = raw_df.copy()
    raw['Activity'] = pd.to_numeric(raw['Activity'], errors='coerce')
    raw['Datetime'] = pd.to_datetime(
        raw['Date'].astype(str) + ' ' + raw['Time'].astype(str), 
        dayfirst=True, 
        errors='coerce',
        format='mixed' # This works ok with the Actiware data, please check if using other data
    )
    
    if 'Interval Status' in raw.columns:
        raw['Interval Status'] = raw['Interval Status'].astype(str).str.upper()


                      
    # Impuation steps. 
                      
    def impute_missing_data(df, imputation_method, fallback_value):
        df_processed = df.copy()
        if imputation_method and 'Interval Status' in df.columns:
            if imputation_method == 'mean':
                impute_values = df[df['Interval Status'] != 'EXCLUDED'].groupby('Time')['Activity'].mean()
            elif imputation_method == 'median':
                impute_values = df[df['Interval Status'] != 'EXCLUDED'].groupby('Time')['Activity'].median()
            else:
                return df_processed
            
            excluded_rows = df_processed['Interval Status'] == 'EXCLUDED'
            df_processed.loc[excluded_rows, 'Activity'] = df_processed.loc[excluded_rows, 'Time'].map(impute_values)
            df_processed.loc[excluded_rows, 'Interval Status'] = 'IMPUTED'
            
        if fallback_value is not None:
            df_processed['Activity'] = df_processed['Activity'].fillna(fallback_value)
        return df_processed
    
    raw_for_L5 = impute_missing_data(raw, l5_imputation, l5_fallback)
    raw_for_M10 = impute_missing_data(raw, m10_imputation, m10_fallback)
    raw_for_IS = impute_missing_data(raw, is_imputation, is_fallback)
    raw_for_IV = impute_missing_data(raw, iv_imputation, iv_fallback)
    raw_for_M16L8 = impute_missing_data(raw, m16l8_time_imputation, m16l8_time_fallback)


                      
    # L5 and M10
                      
    def compute_mean_profile(df):
        df_proc = df.dropna(subset=["Datetime", "Activity"]).copy()
        if df_proc.empty:
            return pd.Series(dtype=float)
        mins = (df_proc["Datetime"].dt.hour * 60 + df_proc["Datetime"].dt.minute + df_proc["Datetime"].dt.second / 60.0)
        mean_profile = df_proc.groupby(mins)["Activity"].mean().sort_index()
        
        target_index = np.arange(0, 1440, 0.5)
        mean_profile_interp = mean_profile.reindex(target_index).interpolate(method="linear", limit_direction="both")
        mean_profile_interp.index.name = "minutes"
        return mean_profile_interp

                      
    average_day_L5 = compute_mean_profile(raw_for_L5)
    average_day_M10 = compute_mean_profile(raw_for_M10)
    base_index = pd.to_datetime("2000-01-01") + pd.to_timedelta(average_day_L5.index, unit="m")
    doubled_L5 = pd.concat([pd.Series(average_day_L5.values, index=base_index), pd.Series(average_day_L5.values, index=base_index + pd.Timedelta(days=1))])
    doubled_M10 = pd.concat([pd.Series(average_day_M10.values, index=base_index), pd.Series(average_day_M10.values, index=base_index + pd.Timedelta(days=1))])
    roll_L5 = doubled_L5.rolling(window=int(round(5*60/0.5)), center=True, min_periods=1).mean()
    roll_M10 = doubled_M10.rolling(window=int(round(10*60/0.5)), center=True, min_periods=1).mean()


                      
    ## This was needed due to rare 'tied' results, where just i.e. dividing by 2 to find average can cause bizarre times. 
    def pick_central_time_of_ties(series_48h, is_max=False):
        if series_48h.empty or not series_48h.notna().any(): return pd.NaT, np.nan
        series_min, series_max = series_48h.min(), series_48h.max()
        if np.isclose(series_min, series_max, atol=1e-8): return pd.NaT, float(series_min) ### In case of error, i.e all 0s or similar
        target_value = series_max if is_max else series_min
        target_mask = np.isclose(series_48h.values, target_value, atol=1e-8)
        times_48h = series_48h.index[target_mask]
        reference_midnight = pd.Timestamp("2000-01-01") ### date is just placeholer, no need to update
        mins_24h = np.array([((time - reference_midnight).total_seconds() / 60.0) % 1440.0 for time in times_48h])
        radians = (mins_24h / 1440.0) * (2 * np.pi)
        mean_minutes = (np.arctan2(np.mean(np.sin(radians)), np.mean(np.cos(radians))) / (2 * np.pi)) * 1440.0
        return reference_midnight + pd.to_timedelta(mean_minutes % 1440.0, unit="m"), float(target_value)
        
    L5_midpoint_time, L5_mean = pick_central_time_of_ties(roll_L5, is_max=False)
    M10_midpoint_time, M10_mean = pick_central_time_of_ties(roll_M10, is_max=True)

                      
    def time_to_minutes_phase_corrected_L(time):
        if pd.isna(time): return np.nan
        m = (time - time.normalize()).total_seconds() / 60.0
        return m - 1440 if m >= 720 else m

                      
    def time_to_minutes_phase_corrected_M(time):
        if pd.isna(time): return np.nan
        m = ((time - time.normalize()).total_seconds() / 60.0) - 720
        return m + 1440 if m <= -720 else m

                      
    L5_midpoint_mins = time_to_minutes_phase_corrected_L(L5_midpoint_time)
    M10_midpoint_mins = time_to_minutes_phase_corrected_M(M10_midpoint_time)
    relative_amplitude = (M10_mean - L5_mean) / (M10_mean + L5_mean) if (M10_mean + L5_mean) > 0 else np.nan


                      
    # Interdaily staubility (IS)
                      
    interdaily_stability_value, interdaily_stability_mean = np.nan, np.nan
    raw_for_IS = raw_for_IS.sort_values('Datetime')
                      
    if not raw_for_IS.empty:
        hourly_data = raw_for_IS.set_index('Datetime')['Activity'].resample('1h').agg(['mean', 'count'])
        valid_hours = hourly_data[hourly_data['count'] > 0].copy()
        if len(valid_hours) > 24:
            interdaily_stability_mean = np.average(valid_hours['mean'], weights=valid_hours['count'])
            valid_hours['hour'] = valid_hours.index.hour
            hourly_means = valid_hours.groupby('hour').apply(lambda group: np.average(group['mean'], weights=group['count']), include_groups=False)
            hourly_weights = valid_hours.groupby('hour')['count'].sum()
            variance_total = np.average((valid_hours['mean'] - interdaily_stability_mean)**2, weights=valid_hours['count'])
            variance_hourly = np.average((hourly_means - interdaily_stability_mean)**2, weights=hourly_weights)
            if variance_total > 0: interdaily_stability_value = variance_hourly / variance_total


                      
    # Intradaily variability (IV)
                      
    hourly_activity_iv = raw_for_IV.groupby(np.arange(len(raw_for_IV)) // 120)['Activity'].mean()
    overall_mean_iv = raw_for_IV['Activity'].mean()
    sum_of_variance = ((hourly_activity_iv - overall_mean_iv) ** 2).sum()
    sum_of_squared_differences = (hourly_activity_iv.diff().dropna() ** 2).sum()
    intradaily_variability_value = (len(hourly_activity_iv) * sum_of_squared_differences) / ((len(hourly_activity_iv) - 1) * sum_of_variance) if sum_of_variance > 0 else np.nan


                      
    # M16/L8 (for the period/window calculations)
                      
    average_m16l8_profile = compute_mean_profile(raw_for_M16L8)
    base_index_m16 = pd.to_datetime("2000-01-01") + pd.to_timedelta(average_m16l8_profile.index, unit="m")
    doubled_m16_profile = pd.concat([pd.Series(average_m16l8_profile.values, index=base_index_m16), pd.Series(average_m16l8_profile.values, index=base_index_m16 + pd.Timedelta(days=1))])
    m16_rolling_mean = doubled_m16_profile.rolling(window=int(16*60/0.5), center=True, min_periods=1).mean().loc[base_index_m16]
    l8_rolling_mean = doubled_m16_profile.rolling(window=int(8*60/0.5), center=True, min_periods=1).mean().loc[base_index_m16]
    m16_score_max, m16_peak_time = (m16_rolling_mean.max(), m16_rolling_mean.idxmax()) if not m16_rolling_mean.empty else (np.nan, pd.NaT)
    l8_score_min, l8_trough_time = (l8_rolling_mean.min(), l8_rolling_mean.idxmin()) if not l8_rolling_mean.empty else (np.nan, pd.NaT)
    
    def comp_is_w(raw_data, start_time, hours, global_mean):
        if pd.isna(start_time) or raw_data.empty or pd.isna(global_mean): return np.nan
        start_time_only, end_time_only = start_time.time(), (start_time + pd.Timedelta(hours=hours)).time()
        time_mask = (raw_data['Datetime'].dt.time >= start_time_only) & (raw_data['Datetime'].dt.time < end_time_only) if start_time_only < end_time_only else (raw_data['Datetime'].dt.time >= start_time_only) | (raw_data['Datetime'].dt.time < end_time_only)
        subset_data = raw_data.loc[time_mask]
        if subset_data.empty: return np.nan
        hourly_stats = subset_data.set_index('Datetime')['Activity'].resample('1h').agg(['mean', 'count'])
        valid_hourly_stats = hourly_stats[hourly_stats['count'] > 0]
        if valid_hourly_stats.empty: return np.nan
        denominator_variance_subset = np.average((valid_hourly_stats['mean'] - global_mean)**2, weights=valid_hourly_stats['count'])
        hourly_means_subset = valid_hourly_stats.groupby(valid_hourly_stats.index.hour).apply(lambda group: np.average(group['mean'], weights=group['count']), include_groups=False)
        return np.average((hourly_means_subset - global_mean)**2, weights=valid_hourly_stats.groupby(valid_hourly_stats.index.hour)['count'].sum()) / denominator_variance_subset if denominator_variance_subset > 0 else np.nan
        
    interdaily_stability_m16 = comp_is_w(raw_for_IS, m16_peak_time - pd.Timedelta(hours=8), 16, interdaily_stability_mean)
    interdaily_stability_l8 = comp_is_w(raw_for_IS, l8_trough_time - pd.Timedelta(hours=4), 8, interdaily_stability_mean)
    
    return pd.DataFrame([[
        os.path.basename(csv_path).replace('.csv', ''), float(L5_mean), L5_midpoint_mins, float(M10_mean), M10_midpoint_mins, relative_amplitude, interdaily_stability_value, intradaily_variability_value, m16_score_max, l8_score_min, interdaily_stability_m16, interdaily_stability_l8
    ]], columns=['ID', 'L5_score_mean', 'L5_midpoint_mins_from_midnight_mean', 'M10_score_mean', 'M10_midpoint_mins_from_midnight_mean', 
                 'relative_amplitude', 'interdaily_stability', 'intradaily_variability', 'M16_score', 'L8_score', 'IS_M16_weighted', 'IS_L8_weighted'])



# 'Sleep' metrics. Also runs the circadian. ?should prpbably separate these to enable easier calculation where rest period not determined?

def sleep_metrics(csv, days_to_remove=0, trim_start=True, trim_end=True, **circadian_params):


 
    ## Load data, handles the Actiware header (for other data this will need adapting)
    try:
        rawest = pd.read_csv(csv, sep='an_unlikely_separator', names=['Column'], engine='python', on_bad_lines='skip')
        line_series = rawest['Column'].str.contains('"Line",', na=False)
        header_row = line_series.idxmax() if line_series.any() else 0
        raw = pd.read_csv(csv, header=header_row, on_bad_lines='skip', engine='python')
    except Exception as e:
        print(f"Failed to load {csv}. Error: {e}")
        return None, None


 
    ### Trim initial data. Essentially sets the start of the first 5 hour uniniterrupted (by missingness) period as
    ### the start of the data. This trims the demonstration period etc (and MRI in our participants) at the beginning of the data. 
 
    first_nonerror_epoch = 0
    win_size = 600
    
    if trim_start and 'Interval Status' in raw.columns and len(raw) > win_size:
        is_excluded = raw['Interval Status'].str.contains('EXCLUDED', na=False)
        rolling_sum = is_excluded.rolling(window=win_size).sum().fillna(0)
        first_valid_window_starts = (rolling_sum == 0)
        if first_valid_window_starts.any():
            first_nonerror_epoch = first_valid_window_starts.idxmax()
    raw = raw.iloc[first_nonerror_epoch:].reset_index(drop=True)
    
    last_valid_epoch = len(raw)
    
    if trim_end and 'Interval Status' in raw.columns and len(raw) > win_size:
        is_excluded_rev = raw['Interval Status'].iloc[::-1].reset_index(drop=True).str.contains('EXCLUDED', na=False)
        rolling_sum_rev = is_excluded_rev.rolling(window=win_size).sum().fillna(0)
        first_valid_rev_starts = (rolling_sum_rev == 0)
        if first_valid_rev_starts.any():
            first_valid_rev = first_valid_rev_starts.idxmax()
            last_valid_epoch = len(raw) - first_valid_rev
    raw = raw.iloc[:last_valid_epoch].reset_index(drop=True)


 
    ## Remove n days (in case of first day effects)
 
    if days_to_remove > 0 and not raw.empty:
        try:
            first_epoch_date = raw.loc[0, 'Date']
            first_epoch_time = raw.loc[0, 'Time']
            first_epoch_date = pd.to_datetime(first_epoch_date, dayfirst=True).date()
            first_epoch_time = pd.to_datetime(first_epoch_time).time()
            first_epoch_datetime = pd.Timestamp.combine(first_epoch_date, first_epoch_time)
            # Trim through midday the next day, regardless of start time. Had to choose something, a little arbitary but works ok. 
            generic_midday = datetime.time(12, 0, 0)
            next_midday_datetime = pd.Timestamp.combine(first_epoch_date, generic_midday) + datetime.timedelta(days=1)
            # Calculate how many epochs that covers (note this is only set for 30s epochs)
            next_midday_epochs = (next_midday_datetime - first_epoch_datetime).total_seconds() / 30
            # Add extra (days_to_remove - 1) full days 
            next_midday_epochs = next_midday_epochs + (2880 * (days_to_remove - 1))
            # Trim
            raw = raw.drop(raw.index[0:int(next_midday_epochs)]).reset_index(drop=True)
        except Exception as e:
            print(f"Could not remove days due to error: {e}") ### Only really a problem if input n strange I think
            pass


 
    # Calculate circadian metrics
 
    circadian_results_df = None
 
    if not raw.empty:
        circadian_results_df = circadian_metrics(raw_df=raw, csv_path=csv, **circadian_params)
    if raw.empty or 'Interval Status' not in raw.columns:
        return circadian_results_df, pd.DataFrame()
    raw['Activity'] = pd.to_numeric(raw['Activity'], errors='coerce').fillna(0)


 
    # Add columns for other sleep scoring algorithms (see refs above). Could remove these, as not very useful withoiut calibration. 
 
    raw["Sadeh_wake"] = 0
    raw["Cole_Kripke_wake"] = 0
    raw["Actiware_wake_20"] = 0
    raw["Actiware_wake_40"] = 0
    raw["Actiware_wake_80"] = 0
    
    for i in range(len(raw)): ### This reworks my previously (much longer way of doing this). 
        epoch_window_dict = {f'p{k}': raw.loc[i+k, 'Activity'] if 0 <= i+k < len(raw) else 0 for k in range(-10, 11)}
        
        mean5minutes = np.mean([epoch_window_dict[f'p{k}'] for k in range(-10, 11)])
        nat_window = [epoch_window_dict[f'p{k}'] for k in range(-10, 11)]
        NAT = sum(1 for x in nat_window if 50 <= x < 100)
        logact = np.log(epoch_window_dict['p0'] + 1)
        
        sadeh_std_window = [epoch_window_dict[f'p{k}'] for k in range(-5, 1)]
        stddev_6min = statistics.stdev(sadeh_std_window) if len(sadeh_std_window) > 1 else 0
        PS = 7.601 - (0.065 * mean5minutes) - (1.08 * NAT) - (0.056 * stddev_6min) - (0.703 * logact)
        if PS < 0:
            raw.loc[i, 'Sadeh_wake'] = 1
        
        D = 0.0001 * ( (50 * epoch_window_dict['p-4']) + (30 * epoch_window_dict['p-3']) + 
                       (14 * epoch_window_dict['p-2']) + (28 * epoch_window_dict['p-1']) + 
                       (121 * epoch_window_dict['p0']) + (8 * epoch_window_dict['p1']) + 
                       (50 * epoch_window_dict['p2']) )
        if D >= 1:
            raw.loc[i, 'Cole_Kripke_wake'] = 1
        actiware_score = ( (epoch_window_dict['p0'] * 2) + 
                           ( (epoch_window_dict['p1'] + epoch_window_dict['p2'] + epoch_window_dict['p-1'] + epoch_window_dict['p-2']) * 0.2 ) + 
                           ( (epoch_window_dict['p3'] + epoch_window_dict['p4'] + epoch_window_dict['p-3'] + epoch_window_dict['p-4']) * 0.04 ) )
        if actiware_score > 20: raw.loc[i, 'Actiware_wake_20'] = 1
        if actiware_score > 40: raw.loc[i, 'Actiware_wake_40'] = 1
        if actiware_score > 80: raw.loc[i, 'Actiware_wake_80'] = 1
        
    is_rest = raw['Interval Status'].str.contains('REST', na=False)
    is_sleep = raw['Interval Status'].str.contains('REST-S', na=False)
    is_wake_input = raw['Sleep/Wake'] == 1
    
    raw['datetime'] = pd.to_datetime(raw['Date'] + ' ' + raw['Time'], dayfirst=True, errors='coerce', format='mixed')
    
    rest_starts = raw.loc[is_rest & ~is_rest.shift(1, fill_value=False), 'datetime']
    rest_ends = raw.loc[is_rest & ~is_rest.shift(-1, fill_value=False), 'datetime']
    sleep_starts = raw.loc[is_sleep & ~is_sleep.shift(1, fill_value=False), 'datetime']
    sleep_ends = raw.loc[is_sleep & ~is_sleep.shift(-1, fill_value=False), 'datetime']
    num_rest_pairs = min(len(rest_starts), len(rest_ends))
    rest_start_times = rest_starts.iloc[:num_rest_pairs].reset_index(drop=True)
    rest_end_times = rest_ends.iloc[:num_rest_pairs].reset_index(drop=True)
    num_sleep_pairs = min(len(sleep_starts), len(sleep_ends))
    sleep_start_times = sleep_starts.iloc[:num_sleep_pairs].reset_index(drop=True)
    sleep_end_times = sleep_ends.iloc[:num_sleep_pairs].reset_index(drop=True)
    rest_periods = len(rest_start_times)
    TIB = is_rest.sum()
    Input_WASO = (is_sleep & is_wake_input).sum()
    Input_TST = (is_sleep & ~is_wake_input).sum()
    Sadeh_TST = (is_sleep & (raw['Sadeh_wake'] == 0)).sum()
    Cole_TST = (is_sleep & (raw['Cole_Kripke_wake'] == 0)).sum()
    Actiware_20_TST = (is_sleep & (raw['Actiware_wake_20'] == 0)).sum()
    Actiware_40_TST = (is_sleep & (raw['Actiware_wake_40'] == 0)).sum()
    Actiware_80_TST = (is_sleep & (raw['Actiware_wake_80'] == 0)).sum()
    
    filename = os.path.basename(csv).replace('.csv', '')


 
    ## Calculate sleep metrics
 
    TIB_avg, Input_SE, Input_TST_avg, Input_WASO_avg = [np.nan] * 4
    Sadeh_SE, Cole_SE, Actiware_20_SE, Actiware_40_SE, Actiware_80_SE = [np.nan] * 5
 
    if rest_periods > 0 and TIB > 0:
        TIB_avg = TIB / 2.0 / rest_periods
        Input_SE = 100.0 / TIB * Input_TST
        Sadeh_SE = 100.0 / TIB * Sadeh_TST
        Cole_SE = 100.0 / TIB * Cole_TST
        Actiware_20_SE = 100.0 / TIB * Actiware_20_TST
        Actiware_40_SE = 100.0 / TIB * Actiware_40_TST
        Actiware_80_SE = 100.0 / TIB * Actiware_80_TST
        Input_TST_avg = Input_TST / 2.0 / rest_periods
        Input_WASO_avg = Input_WASO / 2.0 / rest_periods


 
    #### Sleep timing

    total_sleep_time = (sleep_end_times - sleep_start_times).dt.total_seconds() / 60
    midpoint_sleep = sleep_start_times + pd.to_timedelta(total_sleep_time / 2, unit='m')


 
    # Midpoint_sleep_mins_Mean, uses mignight centering as per L5 etc. 
 
    minutes_from_own_midnight = (midpoint_sleep - midpoint_sleep.dt.normalize()).dt.total_seconds() / 60.0
    
    midpoint_sleep_mins_from_midnight = minutes_from_own_midnight.apply(
        lambda m: m - 1440 if m >= 720 else m
    )
        
    circadian_results = pd.DataFrame({'rest_start_times': rest_start_times, 'Total_sleep_time': total_sleep_time, 'midpoint_sleep_mins_from_midnight': midpoint_sleep_mins_from_midnight})
    
    total_sleep_time_SD = circadian_results['Total_sleep_time'].std()
    Midpoint_sleep_mins_SD = circadian_results['midpoint_sleep_mins_from_midnight'].std()
    total_sleep_time_Mean = circadian_results['Total_sleep_time'].mean()
    Midpoint_sleep_mins_Mean = circadian_results['midpoint_sleep_mins_from_midnight'].mean()
    
    sleep_onset_latency = np.nan
 
    if len(rest_start_times) > 0 and len(rest_start_times) == len(sleep_start_times):
        latencies = (sleep_start_times - rest_start_times).dt.total_seconds() / 60
        sleep_onset_latency = latencies[latencies >= 0].mean()
     
    first_day = pd.to_datetime(raw.loc[0, 'Date'], dayfirst=True).strftime('%A') if not raw.empty else 'N/A'
 
    s1 = pd.DataFrame([[
        filename, rest_periods, TIB_avg, Input_TST_avg, Input_SE, Input_WASO_avg,
        Sadeh_SE, Cole_SE, Actiware_20_SE, Actiware_40_SE, Actiware_80_SE,
        total_sleep_time_Mean, total_sleep_time_SD, 
        Midpoint_sleep_mins_Mean, Midpoint_sleep_mins_SD,
        sleep_onset_latency, first_day
    ]], columns = [
        'ID', 'rest periods', 'TIB (avg)', 'Input TST (avg)', 'Input SE', 'Input WASO (avg)',
        'Sadeh SE', 'Cole-Kripke SE', 'Actiware 20 SE', 'Actiware 40 SE', 'Actiware 80 SE',
        'total_sleep_period_time_Mean', 'total_sleep_time_SD', 
        'Midpoint_sleep_mins_Mean', 'Midpoint_sleep_mins_SD',
        'sleep_onset_latency', 'First_day'
    ])
    if circadian_results_df is not None:
        circadian_results_to_join = circadian_results_df.drop(columns=['ID'])
        s1 = pd.concat([s1.reset_index(drop=True), circadian_results_to_join.reset_index(drop=True)], axis=1)
    return (s1, circadian_results)




#### Function to run it all
   # Runs the sleep and circadian metric calculations on a chosen folder, outputs, enables selection of key variables/missingness management etc. 
def SleePy(
    input_folder, 
    output_path=None, 
    days_to_remove=0,
    trim_start=False, 
    trim_end=False,   
    # Circadian Imputation Defaults
    l5_imputation='mean',
    m10_imputation='mean',
    is_imputation=None,
    iv_imputation='median',
    m16l8_time_imputation='mean'
):
    
    # Output path
 
    if output_path is None:
        output_path = os.path.join(input_folder, "sleep_circadian_metrics_results_minusday_v2point2.csv")
     
    # Find the files
 
    csv_pattern = os.path.join(input_folder, "*.csv")
    csv_files = glob.glob(csv_pattern)
    
    if not csv_files:
        print(f"No CSV files found in {input_folder}")
        return
    all_results = []
    print(f"Starting analysis on {len(csv_files)} files...")
    print("--------------------------------------------------")
    for file_path in csv_files:
        file_name = os.path.basename(file_path)
        print(f"Processing: {file_name}...")
        
        try:
            # Pass the new trim variables down to sleep_metrics
            s1, _ = sleep_metrics(
                file_path, 
                days_to_remove=days_to_remove,
                trim_start=trim_start,
                trim_end=trim_end,
                l5_imputation=l5_imputation,
                m10_imputation=m10_imputation,
                is_imputation=is_imputation,
                iv_imputation=iv_imputation, 
                m16l8_time_imputation=m16l8_time_imputation
            )
            
            if s1 is not None:
                s1 = s1.set_index(['ID'])
                all_results.append(s1)
                
        except Exception as e:
            print(f"  FAILED to process {file_name}. Error: {e}")


 
    # Save them all
 
    if all_results:
        final_df = pd.concat(all_results, sort=False)
        final_df.to_csv(output_path)
        print("--------------------------------------------------")
        print(f"Analysis Complete. Results saved to: {output_path}")
    else:
        print("Error; no results generated, please check input data to confirm all requirements are met, see documentation.")
        print("Likely culprits include date/time formatting, and requirement for data column names as per Actiware output data")
        print("SleePy currently requires data in 30 second epochs")
        print("See documentation and example input data file to understand data requirements.") 






### Example, to run:

# SleePy(input_folder='./example_data/', 
#        days_to_remove=1
#        trim_start=True
#       )

