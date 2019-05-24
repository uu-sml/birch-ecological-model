#!/bin/python
"""Processing of the simulation data"""
import json
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
import csv

# Parameters
burnin = 10000
input_file = "input/observation_data.json"
output_file = "output/ecological_model_pmmh.json"
data_folder = "data/"
realization_file = data_folder + "ecological-realization.csv"
histogram_file = data_folder + "ecological-histogram.csv"


with open(output_file, "r") as f:
    data_full = json.load(f)

# Observation data
with open(input_file) as f:
    y = json.load(f)
    y = y["y"]

data_pmmh = data_full[burnin:-1]
N = len(data_pmmh)


c_pmmh = [data['θ']['c'] for data in data_pmmh]

mean_y = np.zeros(len(y))
var_y = np.zeros(len(y))

for i, sample in enumerate(data_pmmh[1:-1]):
    y_prime = np.exp(np.asarray([x['L'] for x in sample["x"]]))
    mean_y += y_prime/N
    var_y += y_prime*y_prime/N

var_y -= mean_y**2


bins, locations = np.histogram(c_pmmh, bins=27, density=True)


def plot():
    fig, (ax1, ax2) = plt.subplots(2, 1)

    ax1.hist(c_pmmh, facecolor='red', bins=27, alpha=0.3, density=True)
    ax1.set_xlabel("c")
    ax1.set_ylabel("Density")

    ax2.plot(y)
    ax2.plot(mean_y)
    ax2.plot(mean_y+3.0*np.sqrt(var_y))
    ax2.plot(mean_y-3.0*np.sqrt(var_y))
    ax2.set_xlabel("t")
    ax2.set_ylabel("n_t")

    plt.show()


def save():
    with open(realization_file, "w") as csv_file:
        fieldnames = ["t", "y", "y_hat", "sigma_y"]
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(fieldnames)
        writer.writerows(zip(range(1975, 1998), y, mean_y, np.sqrt(var_y)))

    with open(histogram_file, "w") as csv_file:
        fieldnames = ["b3", "Pb3"]
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(fieldnames)
        writer.writerows(zip(locations[0:-1], bins))

if __name__ == "__main__":
    plot()
    save()
