# *SleePy*: a Python program for the analysis of sleep and circadian metrics with automated management of missing data

**Author:** Josh King-Robson, University College London

**Version:** 2.2

**Description:** Python package for actigraphy data analysis. *SleePy* will calculate a range of parametric and non-parametric sleep and circadian metrics from actigraphy data. It enables automated management of missing data, and automatic removal of the first *n* days, i.e. where there are concerns of 'first day' effects. *SleePy* was designed to work with data exported from Philips Actiware but with minimal adaptation will work with other time series data.

---

## Table of Contents

1. [Input data requirements](#1-input-data-requirements)
2. [Sleep and circadian metrics](#2-sleep-and-circadian-metrics)
3. [SleePy](#3-sleepy)
   - [Quick start](#31-quick-start)
   - [SleePy function](#32-sleepy-function)
4. [References](#4-references)

---

## 1. Input data requirements

*SleePy* was designed to work with a batch of comma-separated value (CSV) files exported from Philips Actiware, but will, with minimal adaptation, work with other time series data from other sources. The results are written to a CSV file within the same folder as the input CSV files.

The calculation of all sleep metrics requires that the rest period (the period during which the study subject is in bed with the intention of sleeping) has already been defined. Automated detection of the rest period is inaccurate, and we recommend that this is performed manually, ideally with the aid of a prospectively collected sleep diary and an event marker. The calculation of circadian metrics does not require that the rest period has been set.

---

## 2. Sleep and circadian metrics

*SleePy* will calculate a range of sleep and circadian metrics, and includes automated imputation of missing data (see section 4) for the circadian metrics.

### A: Sleep metrics

| **Term** | **Explanation** | **Equation** | **Ref.** |
| :--- | :--- | :--- | :--- |
| Rest interval (RI) | The period spent in bed per night (with the intention of sleeping, as opposed to reading or other activities). | | 1-3 |
| Time in bed (TIB) | The duration (generally in minutes) of the rest interval. | | 1, 2 |
| Sleep period (SP) | The time from sleep onset to sleep offset (mean minutes per night). Where there are multiple awakenings during the night, this period represents the time from sleep onset to the final awakening prior to getting out of bed. | | 2 |
| Total sleep time (TST) | The total duration of time identified as sleep during the sleep period (mean minutes per night). | | 4 |
| Sleep efficiency (SE) | The percentage of time in bed (with the intention of sleeping) spent asleep. | $\mathbf{SE} = 100\left(\frac{TST}{TIB}\right)$ | 3 |
| Wake after sleep onset (WASO) | The duration of time awake (generally in minutes) during the sleep period. | | 2 |
| Sleep onset latency (SOL) | The time it takes to fall asleep. Calculated as the interval from the start of the rest interval to the onset of the sleep period. | $\mathbf{SOL} = RI_{\text{StartTime}} - SP_{\text{StartTime}}$ | 2 |
| Midpoint of sleep (MPS) | The time of the midpoint of the sleep period. | $\mathbf{MPS} = SP_{\text{StartTime}} + \frac{SP_{\text{Duration}}}{2}$ | 2 |

### B: Non-parametric circadian metrics

| **Term** | **Explanation** | **Equation** | **Ref.** |
| :--- | :--- | :--- | :--- |
| Interdaily stability (IS) | The ratio of the variance of the average 24-hour pattern around the mean and the overall variance. Evaluates the similarity of 24-hour activity patterns across days. We used a one-hour epoch length for all non-parametric calculations. To ensure that partial or incomplete hours contribute appropriately to the overall IS calculation, we weighted the variance by the number of complete data points contributing to it. | $\mathbf{IS} = \frac{n\sum_{h=1}^{p}\left(\overline{x_{h}} - \overline{x}\right)^{2}}{p\sum_{i=1}^{n}\left(x_{i} - \overline{x}\right)^{2}}$ <br><br> $\mathbf{IS} = \frac{\frac{n\sum_{h=1}^{p}\omega_{h}\left(\overline{x_{h}} - \overline{x}\right)^{2}}{\sum\omega_{h}}}{\frac{p\sum_{i=1}^{n}\omega_{i}\left(x_{i} - \overline{x}\right)^{2}}{\sum\omega_{i}}}$ | 5-7 |
| Daytime interdaily stability (IS<sub>day</sub>) | Interdaily stability calculated within the most active 16 hours for each individual. The calculations for both IS<sub>daytime</sub> and IS<sub>nighttime</sub> use the same overall mean ($\overline{x}$) as the overall IS calculation, from the full dataset. This prevents nocturnal activity from unduly impacting the calculation. As with the IS calculation, this is weighted to account for the contribution of incomplete hours to the calculation. | $IS_{day} = \frac{\frac{n \sum_{h \in window} \omega_{h}(\overline{x_{h}} - \overline{x})^{2}}{\sum \omega_{h}}}{\frac{p \sum_{i \in window} \omega_{i}(x_{i} - \overline{x})^{2}}{\sum \omega_{i}}}$ | |
| Nighttime interdaily stability (IS<sub>night</sub>) | Interdaily stability calculated within the least active 8 hours for each individual. | $IS_{night} = \frac{\frac{n \sum_{h \in window} \omega_{h}(\overline{x_{h}} - \overline{x})^{2}}{\sum \omega_{h}}}{\frac{p \sum_{i \in window} \omega_{i}(x_{i} - \overline{x})^{2}}{\sum \omega_{i}}}$ | |
| Intradaily variability (IV) | Quantifies the frequency and extent of transitions from rest to activity, and vice versa. | $\mathbf{IV} = \frac{n\sum_{i=2}^{n}\left(x_{i} - x_{i-1}\right)^{2}}{(n-1)\sum_{i=1}^{n}\left(x_{i} - \overline{x}\right)^{2}}$ | 5, 8 |
| M10 period | The most active consecutive 10-hour period. From which both the activity (M10$_{\text{activity}}$) and time (M10$_{\text{time}}$) can be derived. | $\mathbf{M10} = \max(x_{10})$ | 5 |
| L5 period | The least active consecutive 5-hour period. From which both the activity (L5$_{\text{activity}}$) and time (L5$_{\text{time}}$) can be derived. | $\mathbf{L5} = \min(x_{5})$ | 5 |
| Relative amplitude (RA) | Amplitude of the rest activity rhythm. | $\mathbf{RA} = \frac{M10_{\text{activity}} - L5_{\text{activity}}}{M10_{\text{activity}} + L5_{\text{activity}}}$ | 5 |


**Key:**

| Symbol | Definition |
| :--- | :--- |
| $n$ | Total number of data points (hours, for non-parametric calculations) |
| $p$ | Total number of data points per day (generally 24, for a full day) |
| $\overline{x_{h}}$ | Hourly mean (for a specific clock-hour, e.g. 2–3 PM, across all days) |
| $\overline{x}$ | Mean activity of all data points |
| $x_{i}$ | Mean activity for individual hours (1-hour epoch length used for all non-parametric calculations) |
| $x_{10}$ | Mean activity over 10 consecutive hours |
| $x_{5}$ | Mean activity over 5 consecutive hours |
| $t$ | Timepoint |
| $\omega_{h}$ | Total weight for a clock-hour (number of epochs for that hour of day, e.g. 2–3 PM, across all days) |
| $\omega_{i}$ | Total weight for an individual hour (number of epochs for that individual hour) |
| $\in \text{window}$ | The calculation window for $IS_{day}$ (most active 16 hours) or $IS_{night}$ (least active 8 hours) |

---

## 3. SleePy

### 3.1 Quick start

1. **Install Dependencies:** Ensure you have Python installed, then run:

```bash
pip install -r requirements.txt
```

2. **Prepare Data:** Place all Actiware CSV files (30s epochs) to be processed into a single folder.

3. **Run Analysis:**

```python
from sleepy import SleePy

SleePy(input_folder='./your_data_path/')
```

---

### 3.2 SleePy function

```
SleePy(input_folder, output_path=None, days_to_remove=0, l5_imputation='mean',
       m10_imputation='mean', is_imputation=None, iv_imputation='median',
       m16l8_time_imputation='mean')
```

Function for the analysis of sleep and circadian metrics with automated management of missing data. This will process sleep and circadian metrics for all files within a selected folder, after imputing missing data using the optimum (or user selected) imputation method for each metric.

**Parameters:**

**`input_folder`** : *str*
> The directory path for the folder containing the Actiware CSV files to be processed. These should be in 30-second epochs.

**`output_path`** : *str, default None*
> Custom path for the results CSV file. If None, results are saved as `sleep_circadian_metrics.csv` within the input folder.

**`days_to_remove`** : *int, default 0*
> Number of days to exclude from the beginning of the recording. The script trims through to midday of the day following the start to account for 'first-day' effects or to trim for other indications.

**`trim_start`** : *bool, default False*
> Trims data from the start of the recording, up until the first 5-hour continuous period of data without any data points marked as 'EXCLUDED'. This is useful where there is a demonstration period or similar at the beginning of the data resulting in unrepresentative/erroneous data or frequent off-wrist periods.

**`trim_end`** : *bool, default False*
> Trims data from the end of the recording, until the final 5-hour continuous period of data without any data points marked as 'EXCLUDED'. Can be useful in cases where there is erratic data at the end of recordings.

**`l5_imputation`** : *str, default 'mean'*
> Method for imputation of missing data for use in calculation of L5 metrics. Options: `'mean'`, `'median'`, or `None`.

**`m10_imputation`** : *str, default 'mean'*
> Method for imputation of missing data for use in calculation of M10 metrics. Options: `'mean'`, `'median'`, or `None`.

**`is_imputation`** : *str, default None*
> Method for imputation of missing data for use in calculation of interdaily stability (IS). Options: `'mean'`, `'median'`, or `None`.

**`iv_imputation`** : *str, default 'median'*
> Method for imputation of missing data for use in calculation of Intradaily Variability (IV).

**`m16l8_time_imputation`** : *str, default 'mean'*
> Method for imputation of missing data for use in M16 and L8 period/window calculations.

**Returns:**

A comma-separated values (CSV) file containing the following metrics:

| Output | Description |
| :--- | :--- |
| `TST` | Total sleep time (mean minutes per night) |
| `TIB` | Time in bed (mean minutes per night) |
| `SE` | Sleep efficiency (%) |
| `WASO` | Time awake after sleep onset (mean minutes per night) |
| `sleep_onset_latency` | Interval from the start of the rest interval to the onset of the sleep period (mean minutes per night) |
| `sleep_period` | Time from sleep onset to sleep offset (mean minutes per night) |
| `TST_SD` | Total sleep time (standard deviation) |
| `midpoint_sleep` | Midpoint of sleep, mean (minutes past midnight) |
| `midpoint_sleep_SD` | Midpoint of sleep, standard deviation (minutes past midnight) |
| `first_day` | First day of data collected (name of day of the week) |
| `L5_score` | Activity score during the L5 period (mean per epoch) |
| `L5_midpoint` | Time, in minutes past midnight, of the L5 midpoint |
| `M10_score` | Activity score during the M10 period (mean per epoch) |
| `M10_midpoint` | Time, in minutes past midday, of the M10 midpoint |
| `relative_amplitude` | Relative amplitude score |
| `interdaily_stability` | Interdaily stability score |
| `intradaily_variability` | Intradaily variability score |
| `M16_score` | Activity score during the M16 period (mean per epoch) |
| `L8_score` | Activity score during the L8 period (mean per epoch) |
| `IS_M16_weighted` | Daytime interdaily stability score (during the most active 16 hours) |
| `IS_L8_weighted` | Nighttime interdaily stability score (during the least active 8 hours) |

The following are additionally output, but are **experimental** and should be interpreted with caution, in particular due to the fact that the thresholds for wake/sleep detection are device-dependent:

| Output | Description |
| :--- | :--- |
| `Sadeh_SE` | Sleep efficiency (%), with sleep/wake determined using the Sadeh method |
| `Cole-Kripke_SE` | Sleep efficiency (%), with sleep/wake determined using the Cole-Kripke method |
| `Actiware_20_SE` | Sleep efficiency (%), Actiware methodology (sensitivity threshold 20) |
| `Actiware_40_SE` | Sleep efficiency (%), Actiware methodology (sensitivity threshold 40) |
| `Actiware_80_SE` | Sleep efficiency (%), Actiware methodology (sensitivity threshold 80) |

**Example usage:**

```python
SleePy(
    input_folder='path_to_folder_containing_csvs',
    output_path='path_to_output_results.csv',
    days_to_remove=1,
    is_imputation='mean'
)
```

---

## 4. References

1. Chow CM, Wong SN, Shin M, *et al*. Defining the rest interval associated with the main sleep period in actigraph scoring. *Nat Sci Sleep*. 2016;8:321-328. doi:10.2147/nss.S114969

2. Fekedulegn D, Andrew ME, Shi M, Violanti JM, Knox S, Innes KE. Actigraphy-Based Assessment of Sleep Parameters. *Ann Work Expo Health*. 2020;64(4):350-367. doi:10.1093/annweh/wxaa007

3. Reed DL, Sacco WP. Measuring Sleep Efficiency: What Should the Denominator Be? *J Clin Sleep Med*. 2016;12(2):263-6. doi:10.5664/jcsm.5498

4. Aili K, Åström-Paulsson S, Stoetzer U, Svartengren M, Hillert L. Reliability of Actigraphy and Subjective Sleep Measurements in Adults: The Design of Sleep Assessments. *Journal of Clinical Sleep Medicine*. 2017;13(01):39-47. doi:10.5664/jcsm.6384

5. Cespedes Feliciano EM, Quante M, Weng J, *et al*. Actigraphy-Derived Daily Rest-Activity Patterns and Body Mass Index in Community-Dwelling Adults. *Sleep*. 2017;40(12). doi:10.1093/sleep/zsx168

6. Cavalcanti-Ferreira P, Berk L, Daher N, *et al*. A nonparametric methodological analysis of rest-activity rhythm in type 2 diabetes. *Sleep Sci*. 2018;11(4):281-289. doi:10.5935/1984-0063.20180044

7. Witting W, Kwa IH, Eikelenboom P, Mirmiran M, Swaab DF. Alterations in the circadian rest-activity rhythm in aging and Alzheimer's disease. *Biol Psychiatry*. 1990;27(6):563-72. doi:10.1016/0006-3223(90)90523-5

8. Gonçalves BSB, Cavalcanti PRA, Tavares GR, Campos TF, Araujo JF. Nonparametric methods in actigraphy: An update. *Sleep Science*. 2014;7(3):158-164. doi:10.1016/j.slsci.2014.09.013

9. Cole RJ, Kripke DF, Gruen W, Mullaney DJ, Gillin JC. Automatic sleep/wake identification from wrist activity. *Sleep*. 1992;15(5):461-9. doi:10.1093/sleep/15.5.461

## 5. Licence
Copyright (c) 2026 Josh King-Robson.  
Source code is licensed under the PolyForm Noncommercial License 1.0.0.  
For commercial or other for-profit licensing enquiries, contact j.king-robson@ucl.ac.uk.
