from matplotlib.colors import LinearSegmentedColormap as Cmap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from os.path import getsize
import pandas as pd
import numpy as np

plt.rcParams.update({'font.size': 15, "font.family": "serif", "figure.dpi": 200})
np.random.seed(0)
TH_COLOR = "#005B99"
BOXWIDTH = 5


def main():
    plot_method()

    def add_db_sizes(summary_file, db_sizes_file) -> pd.DataFrame:
        summary = pd.read_csv(summary_file)
        summary.rename(columns=dict(n_Peaks="N_Peaks"), inplace=True)
        dbstats = pd.read_csv(db_sizes_file)
        dbstats["DB_Size"] = dbstats["DB_Size"].str.replace("M", "").astype(int)
        dbparams = dbstats[["Winsize", "n_Peaks", "Overlap"]].apply(tuple, axis=1)
        db_sizes = []
        for w, n, o in summary[["Window_Size", "N_Peaks", "Overlap"]].apply(tuple, axis=1):
            db_sizes.append(dbstats["DB_Size"][dbparams == (w, n, o)].iloc[0])
        summary["DB_Size_MB"] = db_sizes
        return summary

    summary = add_db_sizes("summary.csv", "runtimes_createdb.csv")
    plot_single_prot_matching(summary, "plot_previous")


def plot_single_prot_matching(summary: pd.DataFrame, out_file):
    plt.clf()
    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    set_spines_blue(ax)

    ax.axhline(getsize("../materials/protein.fa") / (1024 ** 2))

    cmap = Cmap.from_list(
        'custom_cmap',
        [(0, 'black'), (1/3, 'black'), (2/3, 'red'), (8/10, 'orange'), (1, 'green')]
    )
    norm = Normalize(vmin=0, vmax=summary["Self_Matches"].max())

    def add_gradient_box(ax, position, color_val, ylim):
        gradient = np.linspace(0, 1, 256).reshape(1, -1)
        gradient = np.vstack((gradient, gradient))
        ax.imshow(
            gradient,
            aspect='auto',
            cmap=Cmap.from_list(
                'custom_cmap2',
                [(0, 'white'), (0, cmap(norm(color_val))), (1, cmap(norm(color_val))), (1, 'white')]
            ),
            extent=[position - .5 * BOXWIDTH, position + .5 * BOXWIDTH, *ylim],
            # alpha=0.2
        )

    gradient_boxes = []
    for win_size, data in sorted(summary.groupby("Window_Size"), key=lambda x: x[1]["DB_Size_MB"].max(), reverse=True):
        box = ax.boxplot([data["DB_Size_MB"]], positions=[win_size],
                         widths=BOXWIDTH,
                         **{i: dict(linewidth=3) for i in ["boxprops", "whiskerprops", "capprops"]},
                         medianprops=dict(linewidth=2),
                         showfliers=False
                         )["boxes"][0]
        gradient_boxes.append((box, data["Unique_Self_Matches"].mean()))
    ax.set_xlabel("Intervallgröße")

    ylim = ax.get_ylim()
    for b, matches in gradient_boxes:
        path = b.get_path()
        vertices = path.vertices
        x0, y0 = vertices[0]
        x1, y1 = vertices[2]
        add_gradient_box(ax, (x1 + x0) / 2, matches, (y0, y1))
    ax.autoscale(enable=True, axis='x', tight=True)
    ax.set_ylim(ylim)
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(xmin - .25 * BOXWIDTH, xmax + .25 * BOXWIDTH)
    plt.colorbar(ScalarMappable(cmap=cmap, norm=norm), ax=plt.gca()).set_label("Eindeutige Wiedererkennung", rotation=-90, va="bottom")

    ax.set_ylabel("Datenbankgröße (MB)")
    ax.set_title("Matching Ergebnisse", fontdict={"fontweight": "bold"})
    plt.savefig("%s.png" % out_file, bbox_inches='tight')


def plot_method():
    plt.clf()
    # plt.rcParams.update({'font.size': 15, "font.family": "serif"})

    fig, (peak_sel_plt, const_map_plt) = plt.subplots(2, 1, figsize=(14, 14), gridspec_kw={'height_ratios': [1, 2]})
    set_spines_blue(peak_sel_plt)
    set_spines_blue(const_map_plt)

    signals = []
    for window in np.random.rand(20, 15):
        signals.append([max(window) * 1.25] + list(window))

    peaks = find_peaks(signals[0])[0]
    for peak in peaks:
        peak_sel_plt.annotate("", xy=(round(peak / 30, 2), signals[0][peak]), xytext=(round(peak / 30, 2), signals[0][peak] + .2), arrowprops=dict(arrowstyle="->", color="red"))
    peak_sel_plt.scatter([round(i / 30, 2) for i in range(len(signals[0]))], signals[0], color=TH_COLOR)
    peak_sel_plt.set_yticklabels([])
    peak_sel_plt.set_yticks([])
    peak_sel_plt.set_xticks([round(i / 30, 2) for i in range(len(signals[0]))])
    peak_sel_plt.set_xlabel("Frequenz")
    peak_sel_plt.set_ylabel("Signalstärke")
    peak_sel_plt.set_title("Frequenzselektion", fontdict={"fontweight": "bold"})

    window_idx = []
    sel_peaks = []
    peak_counts = [0] + list(np.random.randint(0, 5, 19))
    for i, (window, n_peaks) in enumerate(zip(signals, peak_counts), 1):
        p = list(find_peaks(window)[0])
        p = sorted(p, key=lambda x: window[x], reverse=True)[:n_peaks] if n_peaks else p
        sel_peaks.append(p)
        window_idx += len(p) * [i]
    const_map_plt.scatter(window_idx, np.concatenate(sel_peaks) / 30, color=TH_COLOR)

    # draw edges
    edges = {}
    for anker_i, peaks in enumerate(sel_peaks[-9:], len(sel_peaks) - 9):
        for anker in peaks:
            for i, other_peaks in enumerate(sel_peaks[anker_i + 1:], anker_i + 2):
                for other in other_peaks:
                    diff = i - anker_i
                    edges[diff, anker, other] = anker_i
    window_idx = []
    target_zone = []
    for (diff, anker, other), anker_i in edges.items():
        other_i = anker_i + diff
        window_idx.extend((anker_i + 1, other_i))
        target_zone.extend((round(anker / 30, 2), round(other / 30, 2)))

    const_map_plt.plot(window_idx, target_zone, linewidth=.5, markersize=0, alpha=.4)
    y_lo, y_hi = const_map_plt.get_ylim()
    const_map_plt.set_yticks([round(i / 30, 2) for i in range(int(y_lo * 30) + 1, int(y_hi * 30) + 1)])
    const_map_plt.set_xticks(range(1, len(signals) + 1))
    const_map_plt.set_xticklabels(["%s—%s" % (i + 1, i + 30) for i in range(0, 300, 15)], rotation=-90)
    const_map_plt.set_xlabel("Intervall")
    const_map_plt.set_ylabel("Selektierte Frequenz")
    const_map_plt.set_title("Constellation-Map und Hashing", fontdict={"fontweight": "bold"})

    plt.savefig("plot_method.png", bbox_inches='tight')


def set_spines_blue(ax):
    for spine in ax.spines.values():
        spine.set_edgecolor(TH_COLOR)


if __name__ == '__main__':
    main()
