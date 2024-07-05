import pandas as pd
import geopandas as gpd


def addRegressorWeightedAverage(
    df: gpd.GeoDataFrame, regressor_df: gpd.GeoDataFrame, discr_id_col: str = "h3_index"
) -> gpd.GeoDataFrame:
    """
    Apply parameters from regressor database to geographic discretization as a
    weighted average of the areas of intersection.
    Parameters
        df : gpd.GeoDataFrame - geodataframe containing the final geographic
        discretization and a columns 'discr_id' containing unique IDs.
        regressor_df : gpd.GeoDataFrame - geodataframe containing the desired
        regressors and geometries defining geographic boundaries.
        discr_id_col : str - ID column from df
    Returns
        overlay_df : gpd.GeoDataFrame - geodataframe with same structure as df
        with new columns of applied regressors.
    """

    # merge calculating area overlays
    overlay_df = gpd.overlay(df, regressor_df, how="intersection").rename(
        {"geometry": "geo_overlay"}, axis=1
    )
    # merge with df to get original discretization
    overlay_df = pd.merge(overlay_df, df, on=discr_id_col, how="outer").rename(
        {"geometry": "geometry_discr"}, axis=1
    )

    # complete area zero for non overlay dicretizations
    overlay_df["geo_overlay_area"] = overlay_df["geo_overlay"].area.copy()
    overlay_df.loc[overlay_df["geo_overlay"].isna(), "geo_overlay_area"] = 0

    # calculate percentage of areas on each overlay
    overlay_df["overlay_area_percentage"] = (
        overlay_df["geo_overlay_area"] / overlay_df["geometry_discr"].area
    )

    # apply regressor columns multiplying merged values with
    # calculated percentages
    regressor_cols = [
        col for col in regressor_df.columns if not (col in ["geometry", "regr_id"])
    ]
    for regressor_col in regressor_cols:
        overlay_df[regressor_col] *= overlay_df["overlay_area_percentage"]

    # sum weithed partials and recover original columns
    overlay_df = overlay_df.groupby([discr_id_col])[regressor_cols].sum()
    overlay_df = pd.merge(df, overlay_df.reset_index(), on=discr_id_col, how="left")

    return overlay_df


def addRegressorUniformDistribution(
    df: gpd.GeoDataFrame,
    regressor_df: gpd.GeoDataFrame,
    discr_id_col: str = "h3_index",
    type_geo_variable: str = "feature",
) -> gpd.GeoDataFrame:
    """
    Apply parameters from regressor database to geographic discretization as a
    uniform distribution of parameters in regressors area.
    Parameters
        df : gpd.GeoDataFrame - geodataframe containing the final geographic
        discretization and a columns 'discr_id' containing unique IDs.
        regressor_df : gpd.GeoDataFrame - geodataframe containing the desired
        regressors and geometries defining geographic boundaries.
        discr_id_col : str - discretization ID column from df
        type_geo_variable: str - "feature" or "area". specifies if we want to multiply the value of the feature by the
        intersection percentage or just get its area.
    Returns
        overlay_df : gpd.GeoDataFrame - geodataframe with same structure as df
        with new columns of applied regressors.
    """

    # create regressors ID column
    regressor_df["regr_id"] = list(range(len(regressor_df)))

    # merge calculating area overlays
    overlay_df = gpd.overlay(df, regressor_df, how="intersection").rename(
        {"geometry": "geo_overlay"}, axis=1
    )
    # merge with df to get original discretization
    overlay_df = pd.merge(
        overlay_df, regressor_df[["regr_id", "geometry"]], on="regr_id", how="left"
    ).rename({"geometry": "geometry_regr"}, axis=1)

    # calculate percentage of regressor areas on each overlay
    areas = overlay_df["geometry_regr"].area.apply(
        lambda x: 10**-4 if x < 10**-4 else x
    )
    overlay_df["overlay_area_percentage"] = overlay_df["geo_overlay"].area / areas

    # apply regressor columns multiplying merged values with
    # calculated percentages
    regressor_cols = [
        col for col in regressor_df.columns if not (col in ["geometry", "regr_id"])
    ]
    for regressor_col in regressor_cols:
        if type_geo_variable == "feature":
            overlay_df[regressor_col] *= overlay_df["overlay_area_percentage"]
        elif type_geo_variable == "area":
            overlay_df[regressor_col] *= overlay_df["geo_overlay"].area
            overlay_df[regressor_col] /= 10**6

    # sum weithed partials and recover original columns
    overlay_df = overlay_df.groupby([discr_id_col])[regressor_cols].sum()

    overlay_df = pd.merge(df, overlay_df.reset_index(), on=discr_id_col, how="left")

    # If a region in D1 doesn't have intersection with regions in D2, its regressor_cols values may be nan
    # In this case we're replacing them with zero
    for regressor_col in regressor_cols:
        overlay_df[regressor_col] = overlay_df[regressor_col].fillna(0)

    return overlay_df
