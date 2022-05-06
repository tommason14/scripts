#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import re


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
    logfile = sys.argv[1]
    df = read(logfile)
    column = column_from_command_line_pipe(df)
    if column:
        df.set_index("Step")[column].plot().set(xlabel="Step", ylabel=column)
        plt.tight_layout()
        plt.show()
    else:
        # ask user for column from df - omit Step (the first column)
        choices = [i for i in df.columns[1:]]
        print("Please select a column from the following list:")
        for i, choice in enumerate(choices, 1):
            print(f"{i}: {choice}")
        choice = int(input("Enter number: ")) - 1  # -1 because of Step omission
        df.set_index("Step")[choices[choice]].plot().set(
            xlabel="Step", ylabel=choices[choice]
        )
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
