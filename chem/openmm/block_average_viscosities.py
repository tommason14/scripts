#!/usr/bin/env python3

import argparse
import pandas as pd
import sys

FEMTOSECONDS_TO_NANOSECONDS = 1e-15 / 1e-9


def command_line_args():
    parser = argparse.ArgumentParser(
        description="""\
Average viscosities over various simulation lengths. 
Uses output from the openmm-velocity-verlet plugin that computes viscosity as a
function of applied cosine-acceleration.

Example usage:
block_average_viscosities.py -f viscosity.txt viscosity2.txt -b 0 10 0 20

To take in values from two log files and compute average viscosities from 0-10 ns and 0-20 ns
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-f",
        "--files",
        help="Output files containing cosine acceleration data from OpenMM",
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "-b",
        "--blocks",
        help="Space-separated regions to average over, in ns",
        nargs="+",
        type=float,
        required=True,
    )
    parser.add_argument(
        "-dt", "--timestep", help="Timestep in fs. Default = 1", default=1, type=int
    )
    return parser.parse_args()


def responsive_table(data, strings, min_width=13, decimal_places=5):
    """
    Returns a table that is responsive in size to every column.
    Requires a dictionary to be passed in, with the keys referring to
    the headers of the table.
    Also pass in the number of each column that should be a string, starting
    from 1.

    Usage:
        >>> d = {'col1': [1,2,3,4],
                 'col2': ['One', 'Two', 'Three', 'Four']}
        >>> responsive_table(d, strings = [2])

    Can also give a minimum width, defaults to 13 spaces. A `decimal_places` parameter
    can be passed in to define the number of decimal places of floats.
    """
    num_cols = len(data.keys())
    content = zip(*[data[key] for key in data.keys()])  # dict values into list of lists
    # unknown number of arguments
    max_sizes = {}
    try:
        for k, v in data.items():
            max_sizes[k] = len(max([str(val) for val in v], key=len))
    except ValueError:
        print("Error: No data is passed into responsive_table")
        sys.exit(1)

    # create the thing to pass into .format()- can't have brackets like zip gives
    formatting = []
    index = 0
    all_sizes = []
    for val in zip(data.keys(), max_sizes.values()):
        entry, size = val
        if size < min_width or index + 1 not in strings:
            size = min_width
        # also check dict key length
        if len(entry) > size:
            size = len(entry)
        formatting += [entry, size]
        all_sizes.append(size)
        index += 1
    line_length = sum(all_sizes) + num_cols * 3 - 1  # spaces in header
    print("+" + "-" * line_length + "+")
    output_string = "|" + " {:^{}} |" * len(data.keys())
    print(output_string.format(*formatting))
    print("+" + "-" * line_length + "+")
    for line in content:
        formatting = []
        for val in zip(line, all_sizes):
            entry, size = val
            # if not isinstance(entry, str):
            if isinstance(entry, float):
                size = f"{size}.{decimal_places}f"
            formatting.append(entry)
            formatting.append(size)
        print(output_string.format(*formatting))
    print("+" + "-" * line_length + "+")


def main():
    args = command_line_args()
    blocks = [
        [args.blocks[i], args.blocks[i + 1]] for i in range(0, len(args.blocks), 2)
    ]
    df = pd.concat(pd.read_csv(f, sep="\s+") for f in args.files)
    df = df.rename(columns={'#"Step"': "Step", "1/Viscosity (1/Pa.s)": "Visc"})
    # reset to 0 and convert from steps to ns
    # first step recorded is 1000 steps after the start
    df["Step"] = df["Step"] - df["Step"].min() + 1000
    df["Time"] = df["Step"] * args.timestep * FEMTOSECONDS_TO_NANOSECONDS

    results = {"Start (ns)": [], "End (ns)": [], "Viscosity (mPa.s)": []}
    for block in blocks:
        start, end = block
        subset = df[df["Time"].between(start, end)]
        mean = 1 / subset["Visc"].mean()
        std = 1 / subset["Visc"].std()
        results["Start (ns)"].append(start)
        results["End (ns)"].append(end)
        results["Viscosity (mPa.s)"].append(f"{mean*1e3:.3f} +/- {std*1e3:.3f}")

    responsive_table(results, strings=[3], decimal_places=1)
    print("Results saved in visc_averages.csv")
    pd.DataFrame(results).to_csv("visc_averages.csv", index=False)


if __name__ == "__main__":
    main()
