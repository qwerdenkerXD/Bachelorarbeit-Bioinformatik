from matplotlib.colors import LinearSegmentedColormap as Cmap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from os.path import getsize
import pandas as pd
import numpy as np

np.random.seed(0)
TH_COLOR = "#005B99"
plt.rcParams.update({'font.size': 15, "font.family": "serif"})


def main():
    plot_method()
    plot_scoring()
    plot_single_prot_matching(pd.read_csv("../material/Target-Zone/summary.csv"), "plot_target_zone")
    plot_single_prot_matching(pd.read_csv("../material/Selection-Method/summary.csv"), "plot_selection_method")

    summary = pd.read_csv("../material/Filter Hashes/summary.csv")
    plot_single_prot_matching(summary[~summary["Hash_Use_First_Appearance"]], "plot_filter_hashes")

    summary = pd.read_csv("../material/UniRef90 Sampling/summary.csv")
    dbstats = pd.read_csv("../material/UniRef90 Sampling/runtimes_createdb.csv")
    dbstats["DB_Size"] = dbstats["DB_Size"].str.replace("M", "").astype(int)
    dbparams = dbstats[["Winsize", "n_Peaks", "Overlap"]].apply(tuple, axis=1)
    db_sizes = []
    for w, n, o in summary[["Window_Size", "N_Peaks", "Overlap"]].apply(tuple, axis=1):
        db_sizes.append(dbstats["DB_Size"][dbparams == (w, n, o)].iloc[0])
    summary["DB_Size_MB"] = db_sizes
    plot_single_prot_matching(summary, "plot_uniref90")

    plot_family_matching(pd.read_csv("../material/Selection-Method/summary_match_family.csv"), "plot_selection_method")

    summary = pd.read_csv("../material/Filter Hashes/summary_match_family.csv")
    plot_family_matching(summary[~summary["Hash_Use_First_Appearance"]], "plot_filter_hashes")


def plot_single_prot_matching(summary: pd.DataFrame, out_file):
    plt.clf()
    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    set_spine_color(ax)

    ax.axhline(getsize("../material/.protein.fa") / (1024 ** 2))

    cmap = Cmap.from_list(
        'custom_cmap',
        [(0, 'black'), (1/3, 'black'), (2/3, 'red'), (1, 'lime')]
    )
    norm = Normalize(vmin=0, vmax=summary["Self_Matches"].max())

    boxwidth = 5

    def add_gradient_box(ax, position, color_val, ylim):
        gradient = np.linspace(0, 1, 256).reshape(1, -1)
        gradient = np.vstack((gradient, gradient))
        ax.imshow(
            gradient,
            aspect='auto',
            cmap=Cmap.from_list(
                'custom_cmap2',
                [(0, 'white'), (.5, cmap(norm(color_val))), (1, 'white')]
            ),
            extent=[position - .5 * boxwidth, position + .5 * boxwidth, *ylim],
            alpha=0.2
        )

    boxes = []

    def boxplot(ax, position, data):
        box = ax.boxplot([data["DB_Size_MB"]], positions=[position],
                         widths=boxwidth,
                         boxprops=dict(linewidth=3),
                         **{i: dict(linewidth=3) for i in ["whiskerprops", "capprops"]},
                         showfliers=False
                         )["boxes"]
        boxes.append((data.get("Mean_Sharpness", np.zeros(1)).mean(), *box))
        gradient_boxes.append((position, data["Unique_Self_Matches"].mean()))

    gradient_boxes = []
    special_params = tuple(i for i in ("Significance", "Quantile", "Skip_First_K_Freqs") if i in summary)
    if len(special_params):
        groups = summary.groupby(["Window_Size", *special_params])
        special_params = ("Fenstergröße",
            *(i.replace("le", "l").replace("cance", "kanz").replace("Skip_First_K_Freqs", "k") for i in special_params)
        )

        plt.clf()
        fig, ax = plt.subplots(1, 1, figsize=(14 + 6 * (len(groups) // 10), 6))
        set_spine_color(ax)

        ax.axhline(getsize("../material/.protein.fa") / (1024 ** 2))

        xticklabels = []
        for pos, (params, data) in enumerate(groups):
            boxplot(ax, pos * 10, data)
            xticklabels.append("\n".join(str(i) for i in params))
        ax.set_xticklabels(xticklabels)
        ax.set_xlabel("\n".join(special_params))
    else:
        for win_size, data in sorted(summary.groupby("Window_Size"), key=lambda x: x[1]["DB_Size_MB"].max(), reverse=True):
            boxplot(ax, win_size, data)
        ax.set_xlabel("Fenstergröße")

    for i, (sharpness, patch) in enumerate(boxes):
        # Get the box's coordinates
        # patch.set_facecolor("white")
        path = patch.get_path()
        vertices = path.vertices
        x0, y0 = vertices[0]
        x1, y1 = vertices[2]
        # Calculate the height of the filled portion
        box_height = y1 - y0
        filled_height = box_height * sharpness
        # Create a filled rectangle patch
        filled_box = plt.Rectangle((x0, y0), x1 - x0, box_height, color="white")
        ax.add_patch(filled_box)
        filled_box = plt.Rectangle((x0, y0), x1 - x0, filled_height, color=TH_COLOR, alpha=.8)
        ax.add_patch(filled_box)

    ylim = ax.get_ylim()
    for b in gradient_boxes:
        add_gradient_box(ax, *b, ylim)
    ax.autoscale(enable=True, axis='x', tight=True)
    # ax.autoscale(enable=True, axis='y', tight=True)
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(xmin - .25 * boxwidth, xmax + .25 * boxwidth)
    plt.colorbar(ScalarMappable(cmap=cmap, norm=norm), ax=plt.gca()).set_label("Unique Self Matches", rotation=-90, va="bottom")

    ax.set_ylabel("Datenbankgröße (MB)")
    ax.set_title("Single-Protein-Matching Ergebnisse")
    plt.savefig("../results/%s.sp.png" % out_file, bbox_inches='tight')


def plot_family_matching(summary: str, out_file):
    plt.clf()
    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    set_spine_color(ax)

    boxwidth = 5

    special_params = tuple(i for i in ("Significance", "Quantile", "Skip_First_K_Freqs") if i in summary)
    if len(special_params):
        groups = summary.groupby(["Window_Size", *special_params])
        special_params = ("Fenstergröße",
            *(i.replace("le", "l").replace("cance", "kanz").replace("Skip_First_K_Freqs", "k") for i in special_params)
        )

        plt.clf()
        fig, ax = plt.subplots(1, 1, figsize=(14 + 5 * (len(groups) // 10), 6))
        set_spine_color(ax)

        xticklabels = []

        for pos, (params, data) in sorted(enumerate(groups), key=lambda x: x[1][1]["Average_F_Score"].max(), reverse=True):
            ax.boxplot([data["Average_F_Score"]], positions=[pos * 10],
                       widths=boxwidth,
                       **{i: dict(linewidth=3) for i in ["boxprops", "whiskerprops", "capprops"]},
                       showfliers=False
                       )
            xticklabels.append("\n".join(str(i) for i in params))
        ax.set_xticklabels(xticklabels)
        ax.set_xlabel("\n".join(special_params))
    else:
        for win_size, data in sorted(summary.groupby("Window_Size"), key=lambda x: x[1]["Average_F_Score"].max(), reverse=True):
            ax.boxplot([data["Average_F_Score"]], positions=[win_size],
                       widths=boxwidth,
                       **{i: dict(linewidth=3) for i in ["boxprops", "whiskerprops", "capprops"]},
                       showfliers=False
                       )
        ax.set_xlabel("Fenstergröße")

    ax.autoscale(enable=True, axis='x', tight=True)
    # ax.autoscale(enable=True, axis='y', tight=True)
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(xmin - .25 * boxwidth, xmax + .25 * boxwidth)

    ax.set_ylabel("F1-Score")
    ax.set_title("Family-Matching Ergebnisse")
    plt.savefig("../results/%s.fam.png" % out_file, bbox_inches='tight')


def plot_scoring():
    plt.clf()
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))

    edges = [(a, b) for a in range(4) for b in range(4)]
    edges_idx = sorted(np.random.choice(len(edges), 14, replace=False), reverse=True)
    edges = np.array(edges)[edges_idx]
    tp = edges[:10]
    qp = np.concatenate((tp[3:], edges[5:6]))

    # draw edges
    def plot_edges(edges, color, alpha=1, start_i=0, linestyle="-"):
        for i, (y1, y2) in enumerate(edges, start_i):
            ax.plot((i//2, i//2+1), (y1, y2), color=color, alpha=alpha, linestyle=linestyle, linewidth=3)
    plot_edges(tp[:6], "red", .3)
    plot_edges(tp[6:], "red", start_i=4)
    plot_edges(tp[6:], "green", start_i=4, linestyle="--")
    plot_edges(edges[10:], "green", .3, 8)

    first_hit_freq = min(tp[6][0], tp[7][0])
    y_min, y_max = ax.get_ylim()
    # ax.annotate("", xy=(2, y_min), xytext=(2, first_hit_freq), arrowprops=dict(arrowstyle="->", linestyle="--", color="black", alpha=.8))

    # draw right constellation map (Query)
    ax.plot((2, len(edges)//2 -.8), (y_max, y_max), color="green", linewidth=3)
    ax.fill_between((2, len(edges)//2 - .8), y_max, y_min, color="green", alpha=.03)
    ax.text(len(edges)//2 - .8, y_max + .05, "Constellation-Map Suchprotein", ha="right")

    # draw right constellation map (Match)
    ax.plot((-.2, len(tp)//2 - 1), (y_max + .05, y_max + .05), color="red", linewidth=3)
    ax.fill_between((-.2, len(tp)//2 - 1), y_max + .05, y_min, color="red", alpha=.03)
    ax.text(-.2, y_max + .1, "Constellation-Map Trefferprotein")

    ax.text(3.5, first_hit_freq / 2, "S1-Score$=4$", va="center", ha="center")

    _, y_max = ax.get_ylim()
    ax.set_ylim(y_min, y_max + .1)

    for s in ax.spines.values():
        s.set_visible(False)
    ax.spines["bottom"].set_visible(True)

    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_xlabel("Offset")
    # ax.set_ylabel("Frequenz")
    ax.set_title("Ermittlung des S1-Score", fontdict={"fontweight": "bold"})
    plt.savefig("../results/plot_scoring.png", bbox_inches='tight')


def plot_method():
    plt.clf()

    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    set_spine_color(ax)

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


def set_spine_color(ax, color=TH_COLOR):
    for spine in ax.spines.values():
        spine.set_edgecolor(color)


if __name__ == '__main__':
    main()
