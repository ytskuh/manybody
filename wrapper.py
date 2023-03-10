from scipy.integrate import quad
import numpy as np

def command_composer(program, args):
    command = program
    for key in args:
        if args[key]:
            command += f" --{key}={args[key]}"
        else:
            command +=f" --{key}"
    return command

def l1divergence(hist, f, low, high, bins):
    h = (high-low) / bins
    hist_expected = np.empty(bins)
    for i in range(bins):
        hist_expected[i] = quad(f, low+h*i, low+h*i+h)[0]
    return np.abs(hist/hist.sum()-hist_expected).sum()