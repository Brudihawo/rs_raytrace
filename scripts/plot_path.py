import sys
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


def conic(x, radius, conic_param):
    return np.power(x, 2) \
        / (radius * (1 + np.sqrt(1
                                 - (1 + conic_param)
                                 * np.power(x / radius, 2))))


def draw_conic(midpoint, radius, conic_param, max_height, samples=50):
    # input validation
    if (conic_param > -1 and max_height ** 2 > radius**2 / (1 + conic_param)) \
            or radius == 0:
        print("Invalid lens parametrisation")
        return

    xs = np.linspace(-max_height, max_height, samples)
    ys = conic(xs, radius, conic_param)

    plt.plot(ys + midpoint, xs, c="tab:red")


filename = sys.argv[1]

with open(filename, "r") as file:
    sim_config = json.load(file)
    for b in sim_config["ray"]["boundaries"]:
        type = list(b.keys())[0]
        print(f"Btype: {type}")
        b = b[type]

        if type == "Line":
            plt.vlines(b["midpoint"], -b["radius"],
                       b["radius"], color="tab:red")
            plt.text(b["midpoint"], b["radius"] * 0.9,
                     b["opt_idx"])
        elif type == "Spherical":
            draw_circle(b["midpoint"], b["radius"],
                        b["height"], 150)

            if b["radius"] > 0:
                xmax = - np.sqrt(b["radius"] ** 2
                                 - b["height"] ** 2) \
                    + b["midpoint"]
            else:
                xmax = np.sqrt(b["radius"] ** 2
                               - b["height"] ** 2) \
                    + b["midpoint"]

            plt.text(xmax, b["height"] * 0.9,
                     b["opt_idx"])

        elif type == "Conic":
            draw_conic(b["midpoint"], b["radius"],
                       b["conic_param"], b["height"])

            plt.text(conic(b["height"], b["radius"], b["conic_param"]) + b["midpoint"],
                     b["height"] * 0.9,
                     b["opt_idx"])


cols = ["ray_id", "x", "y", "angle"]
df = pd.read_csv("./path.csv", names=cols)


plt.gca().set_aspect('equal', adjustable='box')

for ray_id in df["ray_id"].unique():
    tmp = df[df["ray_id"] == ray_id]
    # if len(tmp) < 6:
    #     continue
    plt.plot(tmp["x"], tmp["y"], color="tab:blue")

plt.show()
