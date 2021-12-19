import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from utils import get_font


def one_plot(df):
    """
    Plot all dipole moments onto one set of axes - has the consequence of
    making distributions that are spread over a large range appear small, when
    paired with narrow distributions
    """
    p = sns.kdeplot(data=df, x="Dtot", hue="resname", fill=True,)
    p.get_legend().set_title(None)
    p.set_xlabel(r"$\mu$ (D)")
    plt.tight_layout()
    plt.savefig("one_plot.pdf", dpi=300)


def faceted(df):
    """
    Plot distributions of each total dipole moment as an individual facet
    """
    p = sns.displot(
        data=df,
        x="Dtot",
        hue="resname",
        kind="kde",
        fill=True,
        col="resname",
        col_wrap=3,
        facet_kws={"sharey": False},
        common_norm=True,
        height=3,
        aspect=4 / 3,
        legend=False,
    )
    p.set_axis_labels(r"$\mu$ (D)", "")
    p.set_titles("{col_name}")
    p.set(yticks=[])
    sns.despine(left=True, right=True, top=True)
    plt.savefig("faceted.pdf", dpi=300)


def main():
    df = pd.read_csv("dipole_vectors.csv")
    df["resname"] = df["resname"].replace({"ch+": "Ch", "dhp": "DHP", "WAT": "Water"})

    # note that rc params are not used in facetgrid functions (i.e. displot)
    sns.set(
        style="ticks",
        font=get_font(),
        palette="rainbow",
        context="talk",
        rc={"figure.figsize": (5, 5)},
    )
    plt.rcParams["mathtext.default"] = "regular"
    one_plot(df)
    faceted(df)


if __name__ == "__main__":
    main()
