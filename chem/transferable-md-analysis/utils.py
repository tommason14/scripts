import sys
import platform


def completion(iterable):
    """
    Write percentage completion of a loop to the screen
    """
    total = len(iterable)
    for num, val in enumerate(iterable):
        yield val
        sys.stdout.write(f"\r{num/total*100:.2f} % complete")
    sys.stdout.write("\n")


def get_font():
    """
    Available sans-serif fonts differ according to the operating system used
    """
    _os = platform.platform().split("-")[0].lower()
    _fonts = {"macos": "Helvetica", "linux": "Nimbus Sans"}
    return _fonts.get(_os, "DejaVu Sans")
