import json
import numpy as np

from mlflow.tracking import MlflowClient


def load_artifacts(run_id, client=None):
    client = (
        client or MlflowClient()
    )  # if you dont provide a client, create one,if you do provide a client, use that one
    # dependency injection : the function is not responsible for creating the client, it just uses it. This makes the function more testable and flexible.
    out = {}
    for name in ["bounds", "mahal", "density_band", "legal"]:
        p = client.download_artifacts(
            run_id, f"constraints/{name}.json"
        )  # it downloads the artifact from mlflow and returns the local path to the downloaded file
        with open(p) as f:
            out[name] = json.load(f)  # it loads the json file and returns a dictionary

    return out


def check(x, bounds, mahal, density_band, legal):
    assert (
        list(x.index) == mahal["feature_order"]
    )  # we want them in the same order as the one used to compute the mahalanobis distance, otherwise the distance will be wrong
    v = x.values

    lo = np.array(
        [bounds[f][0] for f in x.index]
    )  # the 1st percentile of the feature f, for each feature in x.index
    hi = np.array(
        [bounds[f][1] for f in x.index]
    )  # the 99th percentile of the feature f, for each feature in x.index
    # location is the mean of each feature in the training set
    # the _ in location_ and precision_ is a convention used by sklearn to indicate that these attributes are estimated from the training data and are not part of the model parameters.

    diff = v - np.array(mahal["location_"])  #  how far from standard
    d2 = (
        diff @ np.array(mahal["precision_"]) @ diff
    )  # the precision of a measurement is the inverse of its variance: τ = 1/σ²

    c = density_band["coef"]
    pred = (
        c["const"]
        + c["alcohol"] * x["alcohol"]
        + c["residual sugar"] * x["residual sugar"]
    )
    band = abs(x["density"] - pred) / density_band["sigma"]

    so2_cap = (
        legal["total sulfur dioxide"]["cap_high_sugar"]
        if x["residual sugar"] >= 5.0
        else legal["total sulfur dioxide"]["cap_low_sugar"]
    )

    return {
        "L1_bounds": bool(((v >= lo) & (v <= hi)).all()),
        "L2_mahal": bool(d2 <= mahal["threshold_d2"]),
        "L3_density": bool(band <= density_band["k"]),
        "L0_legal": bool(
            x["total sulfur dioxide"] <= so2_cap
            and x["volatile acidity"] <= legal["volatile acidity"]["cap"]
        ),
        "d2": float(d2),
        "z_density": float(band),
    }
