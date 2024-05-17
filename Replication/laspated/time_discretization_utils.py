import numpy as np
import pandas as pd
from datetime import timedelta


def apply_time_frequency(ts: pd.Series, window: str, window_size: int):
    '''
    Auxiliar function that applies multiple time frequencies index.

    Parameters
    ----------
    ts: pandas.Series
        Timeseries that will be tranformed into time frequency indexes

    window: str
        Window frequency.
        Y -> year
        M -> month
        W -> week
        D -> day
        H -> hour
        m -> minute
        s -> second

    window_size: int
        How many windows to aggregate into one index.

    Returns
    -------
    w: pd.Series
        Index series.

    Example
    -------
    >>> ts = pd.to_datetime(pd.Series([
            '2022-01-01', '2022-01-04', '2022-01-10', '2022-12-31'
        ]))
    >>> apply_time_frequency(ts, 'W', 1)
        0     0.0
        1     1.0
        2     2.0
        3    52.0
        dtype: float64
    '''
    # apply window unit

    if window == 'Y':
        w = ts.dt.year - ts.dt.year.min()

    elif window == 'M':
        months = ts.dt.month
        years = ts.dt.year
        min_year = years.min()
        w = 12*(years - min_year) + months - 1

    elif window == 'W':
        min_day = ts.min().replace(hour=0, minute=0, second=0)
        min_day = min_day - timedelta(days=min_day.weekday())
        w = (ts - min_day).dt.total_seconds() // (7*24*60*60)

    elif window == 'D':
        min_day = ts.min().replace(hour=0, minute=0, second=0)
        w = (ts - min_day).dt.total_seconds() // (24*60*60)

    elif window == 'H':
        min_hour = ts.min().replace(hour=0, minute=0, second=0)
        w = (ts - min_hour).dt.total_seconds() // (60*60)

    elif window == 'm':
        min_minute = ts.min().replace(hour=0, minute=0, second=0)
        w = (ts - min_minute).dt.total_seconds() // (60)

    elif window == 's':
        min_second = ts.min().replace(second=0)
        w = (ts - min_second).dt.total_seconds()

    # window_size indicates how many consecutive windows should be aggregated
    return w // window_size


def calculate_seasonality(
    ts: pd.Series,
    seasonality_type: str,
    window: int,
    frequency: int
):
    '''
    Calculates seasonality index from multiple time windows and frequencies.

    Parameters
    ----------
    ts: pandas.Series
        Timeseries to be transformed into seasonality indexes.
        MUST BE sorted.

    seasonality_type: str
        Seasonality frequency type
            Y -> year
            M -> month
            W -> week
            D -> day
            H -> hour
            m -> minute
            S -> second

    window: (int, List[int])
        Seasonality window

    frequency: int
        Seasonality frequency

        eg: If seasonality is between 12 months in a year:
            seasonality_type = 'Y'
            window = 1
            frequency = 12

    Returns
    -------
    season_idx: pandas.Series
        Seasonality indexes

    Example
    -------
    >>> ts = pd.to_datetime(pd.Series([
            '2022-01-01', '2022-01-04', '2022-01-10', '2022-12-31'
        ]))
    >>> calculate_seasonality(ts, 'D', 1, 7)
        0     0
        1     3
        2     2
        3     0
        dtype: float64
    '''
    unitary_seasonality = apply_time_frequency(
        ts,
        seasonality_type,
        1
    )
    # if unique frequency, apply basic rule
    if type(window) is int:
        season_idx = (unitary_seasonality // window) % (frequency // window)

    # if variable frequencies, calculate indexes accordingly
    else:
        season_n = unitary_seasonality % frequency
        season_idx = pd.Series(np.nan, index=season_n.index)

        for i_w, w in enumerate(np.cumsum(window)):
            season_idx.loc[(season_n < w) & season_idx.isna()] = i_w

    return season_idx.astype(int)


def apply_custom_time_events(
    ts: pd.Series,
    time_disc_df: pd.DataFrame,
    nan_idx=0
):

    # list of valid years
    years_list = ts.dt.year.unique()

    # validate column type
    for col in ["start", "end"]:
        if time_disc_df[col].dtype == "O":
            time_disc_df[col] = pd.to_datetime(time_disc_df[col])

    # check for time inversion
    for ind, row in time_disc_df.iterrows():
        if row["end"] < row["start"]:
            raise ValueError(f"Event {ind} has 'end time' before 'start time'")

    # add new rows for events with repetition
    new_rows_list = []
    for ind, rep in time_disc_df["repetition"].items():

        # yearly repetition
        if rep == "yearly":
            sta_dt = time_disc_df.loc[ind, "start"]
            end_dt = time_disc_df.loc[ind, "end"]

            # yearly repetition doesn't support start and end dates
            # on different years
            if sta_dt.year != end_dt.year:
                raise ValueError("Yearly repetition is not supported for start and end dates on different years.")

            # add new row for valid year
            for y in years_list:
                y1 = sta_dt.year
                if y1 != y:
                    new_rows_list.append(pd.Series([
                        sta_dt.replace(year=y),
                        end_dt.replace(year=y),
                        time_disc_df["t"].loc[ind],
                        None
                    ], index=time_disc_df.columns))

        elif rep is not None:
            raise ValueError(f"Repetition {rep} is not supported.")

    full_time_disc_df = pd.concat([
        time_disc_df,
        pd.DataFrame(new_rows_list)
    ], axis=0, ignore_index=True)

    # ajust one-day intervals
    for ind, row in full_time_disc_df.iterrows():
        # date formatted start times
            # if full day interval, ajust end date
            if row["start"] == row["end"]:
                full_time_disc_df.loc[ind, "end"] = row["end"] + timedelta(days=1)

    # create new t index
    tdiscr = pd.Series(nan_idx, index=ts.index)
    for t_idx in full_time_disc_df["t"].unique():
        # filter conditions for idx t and apply
        t_df = full_time_disc_df[full_time_disc_df["t"] == t_idx]
        for _, row in t_df.iterrows():
            tdiscr.loc[
                ((ts >= row["start"]) & (ts < row["end"]))
            ] = t_idx

    return tdiscr
