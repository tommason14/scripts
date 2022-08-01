#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import re
from matplotlib.ticker import FuncFormatter


def formatting(x, pos):
    """
    The string formatter directly manipulates the value passed in, so to convert
    1000000 into 1x10^6, you have to divide the value by 10^6, like so:
    f"{x*10**-6}x10$^6$
    """
    if x == 0:
        return 0
    if x < 1e6:
        # like 9.5x10^5
        return f"{x*10**-5:,.1f}x10$^{5}$"
    # expecting lots of timesteps, so report to 0 dp to allow for more space along axis
    return f"{x*10**-6:,.0f}x10$^{6}$"


def read(fname):
    header = None
    data = []
    found = False
    with open(fname) as f:
        for line in f:
            if "Step" in line:
                found = True
                if "," in line:
                    header = [
                        i.replace('"', "").replace("#", "") for i in line.split(",")
                    ]
                else:
                    header = [i.replace('"', "").replace("#", "") for i in line.split()]
                continue
            if found:
                if "," in line:
                    data.append(line.split(","))
                else:
                    data.append(line.split())

    return pd.DataFrame(data, columns=header).apply(lambda x: pd.to_numeric(x))


def column_from_command_line_pipe(df):
    """
    If user pipes in a partial match to any column, return it.
    For example:
    echo dens | openmm_plot.py nvt.log will return 'Density (g/mL)'
    """
    if not sys.stdin.isatty():
        user_input = sys.stdin.readline().split()
        argument = user_input[0]
        for col in df.columns:
            if re.search(argument.lower(), col.lower()):
                return col
    return None


def main():
    if len(sys.argv) < 2 or "-h" in sys.argv:
        print("Usage: openmm_plot.py <file>")
        print("The column to plot can be piped in from the command line:")
        print("$ echo 'density' | openmm_plot.py npt.log")
        print("Otherwise the script will prompt you for the column to plot.")
        sys.exit(1)
    logfile = sys.argv[1]
    df = read(logfile)
    fig, ax = plt.subplots()
    column = column_from_command_line_pipe(df)
    if column:
        df.set_index("Step")[column].plot(ax=ax)
        ax.set(xlabel="Step", ylabel=column)
    else:
        # ask user for column from df - omit Step (the first column)
        choices = [i for i in df.columns[1:]]
        print("Please select a column from the following list:")
        for i, choice in enumerate(choices, 1):
            print(f"{i}: {choice}")
        choice = int(input("Enter number: ")) - 1  # -1 because of Step omission
        df.set_index("Step")[choices[choice]].plot(ax=ax)
        ax.set(xlabel="Step", ylabel=choices[choice])
    # step number scientific formatting
    ax.xaxis.set_major_formatter(FuncFormatter(formatting))
    # if looking at energy, also format y axis
    if "potential" in ax.yaxis.label._text:
        ax.yaxis.set_major_formatter(FuncFormatter(formatting))
    plt.tight_layout()
    plt.savefig("plot.png", dpi=300)


if __name__ == "__main__":
    main()
