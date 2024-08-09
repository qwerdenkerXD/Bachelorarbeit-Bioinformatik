from matplotlib import pyplot as plt
from scipy.signal import find_peaks
import seaborn as sb
import numpy as np

np.random.seed(0)
# sb.set(font_scale=1.5)
TH_COLOR = "#005B99"


def main():
    plot_method()


def plot_method():
    plt.clf()
    plt.rcParams.update({'font.size': 15, "font.family": "serif"})

    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    set_spines_blue(ax)

    signals = []
    for window in np.random.rand(20, 15):
        signals.append([max(window) * 1.25] + list(window))

    peaks = find_peaks(signals[0])[0]
    for peak in peaks:
        ax.annotate("", xy=(round(peak / 30, 2), signals[0][peak]), xytext=(round(peak / 30, 2), signals[0][peak] + .2), arrowprops=dict(arrowstyle="->", color="red"))
    ax.scatter([round(i / 30, 2) for i in range(len(signals[0]))], signals[0], color=TH_COLOR)
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_xticks([round(i / 30, 2) for i in range(len(signals[0]))])
    ax.set_xlabel("Frequenz")
    ax.set_ylabel("Signalstärke")
    ax.set_title("Frequenzselektion", fontdict={"fontweight": "bold"})
    plt.savefig("../results/plot_frequency_selection.png", bbox_inches='tight')

    ax.clear()

    window_idx = []
    sel_peaks = []
    peak_counts = [0] + list(np.random.randint(0, 5, 19))
    for i, (window, n_peaks) in enumerate(zip(signals, peak_counts), 1):
        p = list(find_peaks(window)[0])
        p = sorted(p, key=lambda x: window[x], reverse=True)[:n_peaks] if n_peaks else p
        sel_peaks.append(p)
        window_idx += len(p) * [i]
    ax.scatter(window_idx, np.concatenate(sel_peaks) / 30, color=TH_COLOR)

    # draw edges
    edges = {}
    for anker_i, peaks in enumerate(sel_peaks[-9:], len(sel_peaks) - 9):
        for anker in peaks:
            for i, other_peaks in enumerate(sel_peaks[anker_i + 1:], anker_i + 2):
                for other in other_peaks:
                    diff = i - anker_i
                    if (prev_anker_i := edges.get((diff, anker, other), False)):
                        ax.plot((prev_anker_i + 1, prev_anker_i + diff), (round(anker / 30, 2), round(other / 30, 2)), linewidth=.5, markersize=0, alpha=1, color="red")
                    edges[diff, anker, other] = anker_i
    window_idx = []
    target_zone = []
    for (diff, anker, other), anker_i in edges.items():
        other_i = anker_i + diff
        window_idx.extend((anker_i + 1, other_i))
        target_zone.extend((round(anker / 30, 2), round(other / 30, 2)))

    ax.plot(window_idx, target_zone, linewidth=.5, markersize=0, alpha=.4)
    y_lo, y_hi = ax.get_ylim()
    ax.set_yticks([round(i / 30, 2) for i in range(int(y_lo * 30) + 1, int(y_hi * 30) + 1)])
    ax.set_xticks(range(1, len(signals) + 1))
    ax.set_xticklabels(["%s—%s" % (i + 1, i + 30) for i in range(0, 300, 15)], rotation=-90)
    ax.set_xlabel("Intervall")
    ax.set_ylabel("Selektierte Frequenz")
    ax.set_title("Constellation-Map und Hashing", fontdict={"fontweight": "bold"})

    plt.savefig("../results/plot_method.png", bbox_inches='tight')


def set_spines_blue(ax):
    for spine in ax.spines.values():
        spine.set_edgecolor(TH_COLOR)


if __name__ == '__main__':
    main()
