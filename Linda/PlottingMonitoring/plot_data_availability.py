
#%%
from mats_utils.data_availability.data_availability import load_hourly_time_spans, plot_time_spans, is_data_available  
from datetime import datetime


available_times_l0 = '/Users/lindamegner/MATS/MATS-retrieval/MATS-utility-functions/data/data_availability/available_times_level0-v0.3.txt'
available_times_l1a = '/Users/lindamegner/MATS/MATS-retrieval/MATS-utility-functions/data/data_availability/available_times_level1a-v1.0.txt'
available_times_l1b = '/Users/lindamegner/MATS/MATS-retrieval/MATS-utility-functions/data/data_availability/available_times_level1b-v1.0.1.txt'
# Load time spans
spans_l0 = load_hourly_time_spans(available_times_l0)
spans_l1a = load_hourly_time_spans(available_times_l1a)
spans_l1b = load_hourly_time_spans(available_times_l1b)

#%%
# Plot them
plot_time_spans(spans_l1b, start_date=datetime(2023, 7, 13), end_date=datetime(2024, 7, 17))
#plot_time_spans(time_spans=spans_l0, time_spans_1=spans_l1a, time_spans_2=spans_l1b)

#%%
# Check a specific datetime
dt_to_check = datetime(2025, 1, 24, 17, 0, 0)
print("Data available:", is_data_available(dt_to_check, spans_l1b, nearest_data_prior=True, nearest_data_after=True))

# # %%
# def is_data_available(check_dt, time_spans, nearest_data_prior=False, nearest_data_after=False):
#     """
#     Check if a given datetime falls within any of the provided time spans.

#     Parameters:
#     ----------
#     check_dt : datetime
#         The datetime to check for data availability.
#     time_spans : list of tuples
#         A list of (start_datetime, end_datetime) tuples representing available data intervals.
#     nearest_data_prior : bool, optional
#         If True and no data is available at check_dt, return the closest prior data end time.
#     nearest_data_after : bool, optional
#         If True and no data is available at check_dt, return the closest future data start time.

#     Returns:
#     -------
#     dict
#         A dictionary with:
#         - 'available': True if data is available at check_dt, False otherwise.
#         - 'nearest_prior': datetime of the closest prior data (if requested and not available).
#         - 'nearest_after': datetime of the closest future data (if requested and not available).
#     """
#     found = False
#     nearest_prior = None
#     nearest_after = None

#     for start, end in time_spans:
#         if start <= check_dt <= end:
#             found = True
#             break
#         if end < check_dt:
#             if nearest_prior is None or end > nearest_prior:
#                 nearest_prior = end
#         if start > check_dt:
#             if nearest_after is None or start < nearest_after:
#                 nearest_after = start

#     result = {"available": found}
#     if not found:
#         if nearest_data_prior:
#             result["nearest_prior"] = nearest_prior
#         if nearest_data_after:
#             result["nearest_after"] = nearest_after

#     return result
# # %%

# %%
