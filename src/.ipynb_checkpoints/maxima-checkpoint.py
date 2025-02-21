import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def getDiffs(arr):
    return arr[1:] - arr[:-1]

def derivative(x, y):
    slope = getDiffs(y) / getDiffs(x)
    dx = x[1] - x[0]
    x_alt = x[:-1] + dx/2
    return x_alt, slope

def argClosest(arr, val):
    return np.argmin(abs(arr - val))

def localMaxima(sites, tum, data):
    fig, ax = plt.subplots()
    plot = sns.histplot(ax=ax, x=data['tumor']['beta_values_SELECTION'].loc[sites, tum], kde=True, bins=40)
    fig.set_visible(False)
    plt.close()
    
    x, y = plot.get_lines()[0].get_data()
    x_prime, y_prime = derivative(x, y)
    x_pprime, y_pprime = derivative(x_prime, y_prime)

    max_dy = np.max(np.abs(y[1:] - y[:-1]))
    extrema = x_prime[np.abs(y_prime) < max_dy*10]
    extrema = extrema[(extrema > 0) & (extrema < 1)]

    maxima = []
    last_concavity = 1
    for e in extrema:
        concavity = y_pprime[argClosest(x_pprime, e)]
        if concavity * last_concavity < 0:    # different signed concavity than previous extremum
            if concavity < -1e2:
#                 plt.axvline(e)
                maxima.append(e)
            last_concavity = concavity
    
    return maxima