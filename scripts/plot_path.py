import json

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


def draw_circle(midpoint, radius, max_height, samples=50):
    if radius > 0:
        xmax = - np.sqrt(radius ** 2 - max_height ** 2) + midpoint
    else:
        xmax = np.sqrt(radius ** 2 - max_height ** 2) + midpoint

    xs = np.linspace(midpoint - radius, xmax, 50)
    ys = np.sqrt(np.power(radius, 2) - np.power(xs - midpoint, 2))

    plt.plot(xs, ys, c="tab:red")
    plt.plot(xs, -ys, c="tab:red")


with open("./boundaries.json", "r") as file:
    for id, boundary in json.load(file).items():
        if boundary["type"] == "line":
            plt.vlines(boundary["midpoint"], -boundary["radius"],
                       boundary["radius"], color="tab:red")
            plt.text(boundary["midpoint"], boundary["radius"] * 0.9,
                     boundary["opt_idx"])
        elif boundary["type"] == "spherical":
            draw_circle(boundary["midpoint"], boundary["radius"],
                        boundary["height"], 150)

            if boundary["radius"] > 0:
                xmax = - np.sqrt(boundary["radius"] ** 2
                                 - boundary["height"] ** 2) \
                    + boundary["midpoint"]
            else:
                xmax = np.sqrt(boundary["radius"] ** 2
                               - boundary["height"] ** 2) \
                    + boundary["midpoint"]

            plt.text(xmax, boundary["height"] * 0.9,
                     boundary["opt_idx"])


cols = ["ray_id", "x", "y", "angle"]
df = pd.read_csv("./path.csv", names=cols)


plt.gca().set_aspect('equal', adjustable='box')

for ray_id in df["ray_id"].unique():
    tmp = df[df["ray_id"] == ray_id]
    # if len(tmp) < 6:
    #     continue
    plt.plot(tmp["x"], tmp["y"], color="tab:blue")

plt.show()
