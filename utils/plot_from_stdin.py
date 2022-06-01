#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import sys
import argparse

parser = argparse.ArgumentParser(
    (
        "Expects all data manipulation to be performed before passing it in to this program. This script simply takes a table of data and plots it using matplotlib"
    )
)
parser.add_argument(
    "--header",
    help="Space-separated list of column names to apply, if the data does not include a header line",
    nargs="+",
)
parser.add_argument(
    "-d",
    "--delimeter",
    help="Delimeter is split by, default is to split on any whitespace",
)
parser.add_argument("-o", "--output", help="png filename to save to")

args = parser.parse_args()

if args.delimeter:
    data = [line.strip("\n").split(args.delimeter) for line in sys.stdin]
else:
    data = [line.strip("\n").split() for line in sys.stdin]

# now cannot accept any more input - so asking for columns will fail - so open up fake stdin to
# allow it
sys.stdin = open("/dev/tty")

if args.header:
    data = pd.DataFrame(data, columns=args.header).apply(lambda x: pd.to_numeric(x))
else:
    data = pd.DataFrame(data).apply(lambda x: pd.to_numeric(x))

print("Choose from the following columns:")
for i, col in enumerate(data.columns, 1):
    print(f"{i}: {col}")
# user reads from 1, but df is 0-indexed
x_col = int(input("Enter the number of x-axis column: ")) - 1
y_col = int(input("Enter the number of y-axis column: ")) - 1

xcol = data.columns[x_col]
ycol = data.columns[y_col]
data.set_index(xcol)[ycol].plot().set(xlabel=xcol, ylabel=ycol)

if args.output:
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
else:
    plt.show()
