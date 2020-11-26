def remove_legend_title(plot):
    """
    Modifies the global state of the plot so return value
    is needed. 
    Use like:
        plot = sns.lineplot(...)
        remove_legend_title(plot)
    """
    handles, labels = plot.get_legend_handles_labels()
    # plot.legend(handles=handles[1:], labels=labels[1:]) # old version of seaborn
    plot.legend(handles=handles, labels=labels)
