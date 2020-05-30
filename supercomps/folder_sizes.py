#!/usr/bin/env python3
#
# import subprocess
# import math
#
# def human_readable_sizes(dirs):
#     """
#     Returns a dictionary with folders as keys and folder sizes as
#     values in a nested dictionary value. Folder sizes found using
#     `du -h -d 1`, or if `-d` is not supported with the version of 
#     du installed, `du -h --max-depth=1` (GNU vs BSD).
#     """
#     try:
#         readable = subprocess.run(['du', '-h', '-d', '1'], encoding = 'utf-8',
#                    stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
#         readable = readable.stdout.split('\n')[:-1]
#
#         for line in readable[:-1]:
#             if 'Operation not permitted' not in line and 'cannot read' not in line:
#                 size, folder = line.split('\t')
#                 folder = folder[2:]
#                 dirs[folder] = {'readable': size.strip()}
#
#     except ValueError:
#         readable = subprocess.run(['du', '-h', '--max-depth=1'], encoding = 'utf-8',
#                    stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
#         readable = readable.stdout.split('\n')[:-1]
#
#         for line in readable[:-1]:
#             if 'Operation not permitted' not in line and 'cannot read' not in line:
#                 size, folder = line.split('\t')
#                 folder = folder[2:]
#                 dirs[folder] = {'readable': size.strip()}
#
#     total_size = readable[-1].split('\t')[0].strip()
#     return dirs, total_size
#
# def readable_to_comparable(dirs):
#     """
#     Saves having to call `du -d 1` to get the sizes in terms of blocks.
#     Instead, parse the output of the `du` command once, then convert size
#     to a floating point value and use that to sort folders.
#     """
#     byte_sizes = {
#         'B': 10 ** 0, 
#         'K': 10 ** 3,  # kilo  
#         'M': 10 ** 6,  # mega  
#         'G': 10 ** 9,  # giga  
#         'T': 10 ** 12, # tera   
#         'P': 10 ** 15, # penta  
#     }
#
#     for k, v in dirs.items():
#         value = v['readable']
#         value = float(value[:-1]) * byte_sizes[value[-1]]
#         v['size'] = value
#     return dirs
#
# def folder_owners(dirs):
#     """
#     Adds the owner of each folder to the output of `folder_sizes` or
#     `human_readable_sizes`.
#     """
#     owners = subprocess.run(['ls', '-la'], encoding = 'utf-8', stdout = subprocess.PIPE)
#     owners = owners.stdout.split('\n')[1:] # remove total = ...
#
#     for line in owners[:-1]:
#         line   = line.split()
#         owner  = line[2]
#         folder = line[-1]
#         for k, v in dirs.items():
#             if folder == k:
#                 v['owner'] = owner
#     return dirs
#
# def sort_dict_by_size(dirs):
#     """
#     Returns a sorted list, largest folder first.
#     """
#     lst = []
#     for k, v in dirs.items():
#         lst.append([k, int(v['size'])])
#
#     lst = sorted(lst, key = lambda kv: kv[1], reverse = True)
#
#     return lst
#
# def print_table(dirs, lst, total):
#     """
#     Print results, using the list of folders sorted
#     by size. The first column of the table is 
#     reponsive to match the longest folder name,
#     avoiding any strange loooking tables with lines of 
#     different lengths.
#     """
#     longest_folder_name = max(dirs.keys(), key = len)
#     folder_col_size = len(longest_folder_name)
#     if folder_col_size < 6:
#         folder_col_size = 6
#     length_dashes = folder_col_size + 31
#     print('+' + '-' * length_dashes + '+')
#     print("| {:^{}} | {:^{}} | {:^{}} |".format(
#         'Folder', folder_col_size,
#         'Owner', 15, 
#         'Size', 8))
#     print('+' + '-' * length_dashes + '+')
#     for d in lst:
#         folder, size = d
#         try:
#             print("| {:^{}} | {:^{}} | {:^{}} |".format(
#                 folder, folder_col_size,
#                 dirs[folder]['owner'], 15, 
#                 dirs[folder]['readable'], 8))
#         except KeyError:
#             continue
#     print('+' + '-' * length_dashes + '+')
#     print(f"\nTotal size: {total}\n")
#
#     per_user = usage_per_user(dirs)
#     print(per_user)
#     per_user = change_size_to_readable(per_user)
#     print('Usage per user:')
#     for user, usage in per_user.items():
#         print(f"{user:>15s}: {usage}")
#     
#
# def print_disclaimer():
#     """
#     The shell command `du` can take a long time
#     to find all folder sizes if directory is large.
#     """
#     print("""\
# This may take a while...
#
# Calling the shell `du` command will 
# take a while if the directory is large.
#
# """)
#
# def change_size_to_readable(data):
#     """
#     Convert the values of a one level dictionary of user: size
#     to human readable sizes. i.e. 12000 bytes-> 12 K
#     """
#
#     def convert_size(size_bytes):
#         if size_bytes == 0:
#             return "0 B"
#         size_name = ("B", "K", "M", "G", "T", "P", "E", "Z", "Y")
#         i = int(math.floor(math.log(size_bytes, 1000)))
#         p = math.pow(1000, i)
#         s = round(size_bytes / p, 2)
#         return "%s %s" % (s, size_name[i])
#     print(data.items())
#     return {k: convert_size(v) for k, v in data.items()}
#     
# def usage_per_user(data):
#     owners = {}
#     for folder, d in data.items():
#         print(folder)
#         try:
#             # here code doesn't work...
#             owner = d['owner']
#         except KeyError:
#             continue
#         if owner not in owners:
#             owners[owner] = 0
#         else:
#             value = d['size']
#             owners[owner] += value
#
#     # sort owners by size
#     tmp = [(k, v) for k, v in owners.items()]
#     tmp = sorted(tmp, key = lambda kv: kv[1], reverse=True)
#     
#     return {k: v for k, v in tmp}
#
# def main():
#     print_disclaimer()
#     dirs = {}
#     dirs, total = human_readable_sizes(dirs)
#     dirs = folder_owners(dirs)
#     dirs = readable_to_comparable(dirs)
#     lst = sort_dict_by_size(dirs)
#     print_table(dirs, lst, total)
#
# if __name__ == '__main__':
#     main()
#
import subprocess
import math
import sys

def exclude_folders(args):    
    """
    Returns all arguments after -e
    """
    for ind, val in enumerate(args):
        if val == '-e':
            return args[ind + 1:]
        

def human_readable_sizes(dirs, exclude):
    """
    Returns a dictionary with folders as keys and folder sizes as
    values in a nested dictionary value. Folder sizes found using
    `du -h -d 1`, or if `-d` is not supported with the version of 
    du installed, `du -h --max-depth=1` (GNU vs BSD).
    """
    def add_excluded_folders(command, exclude):
        if exclude is not None:
            command = command + ['--exclude'] + exclude
        return command
        
    command = ['du', '-h', '-d', '1']
    command = add_excluded_folders(command, exclude)
    readable = subprocess.run(command, encoding = 'utf-8',
               stdout = subprocess.PIPE, stderr = subprocess.STDOUT)

    if any(bad_string in readable.stdout for bad_string in ('invalid', 'not permitted')):     
        command = ['du', '-h', '--max-depth=1']
        command = add_excluded_folders(command, exclude)
        readable = subprocess.run(command, encoding = 'utf-8',
                   stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
   
    readable = readable.stdout.split('\n')[:-1]
    for line in readable[:-1]:
        if 'Operation not permitted' not in line and 'cannot read' not in line:
            size, folder = line.split('\t')
            folder = folder[2:]
            dirs[folder] = {'readable': size.strip()}

    total_size = readable[-1].split('\t')[0].strip()
    return dirs, total_size

def readable_to_comparable(dirs):
    """
    Saves having to call `du -d 1` to get the sizes in terms of blocks.
    Instead, parse the output of the `du` command once, then convert size
    to a floating point value and use that to sort folders.
    """
    byte_sizes = {
        'B': 10 ** 0, 
        'K': 10 ** 3,  # kilo  
        'M': 10 ** 6,  # mega  
        'G': 10 ** 9,  # giga  
        'T': 10 ** 12, # tera   
        'P': 10 ** 15, # penta  
    }

    for k, v in dirs.items():
        value = v['readable']
        value = float(value[:-1]) * byte_sizes[value[-1]]
        v['size'] = value
    return dirs

def folder_owners(dirs):
    """
    Adds the owner of each folder to the output of `folder_sizes` or
    `human_readable_sizes`.
    """
    owners = subprocess.run(['ls', '-la'], encoding = 'utf-8', stdout = subprocess.PIPE)
    owners = owners.stdout.split('\n')[1:] # remove total = ...

    for line in owners[:-1]:
        line   = line.split()
        owner  = line[2]
        folder = line[-1]
        for k, v in dirs.items():
            if folder == k:
                v['owner'] = owner
    return dirs

def sort_dict_by_size(dirs):
    """
    Returns a sorted list, largest folder first.
    """
    lst = []
    for k, v in dirs.items():
        lst.append([k, int(v['size'])])

    lst = sorted(lst, key = lambda kv: kv[1], reverse = True)

    return lst

def print_table(dirs, lst, total):
    """
    Print results, using the list of folders sorted
    by size. The first column of the table is 
    reponsive to match the longest folder name,
    avoiding any strange loooking tables with lines of 
    different lengths.
    """
    longest_folder_name = max(dirs.keys(), key = len)
    folder_col_size = len(longest_folder_name)
    if folder_col_size < 6:
        folder_col_size = 6
    length_dashes = folder_col_size + 31
    print('+' + '-' * length_dashes + '+')
    print("| {:^{}} | {:^{}} | {:^{}} |".format(
        'Folder', folder_col_size,
        'Owner', 15, 
        'Size', 8))
    print('+' + '-' * length_dashes + '+')
    for d in lst:
        folder, size = d
        try:
            print("| {:^{}} | {:^{}} | {:^{}} |".format(
                folder, folder_col_size,
                dirs[folder]['owner'], 15, 
                dirs[folder]['readable'], 8))
        except KeyError:
            continue
    print('+' + '-' * length_dashes + '+')
    print(f"\nTotal size: {total}\n")

    per_user = usage_per_user(dirs)
    per_user = change_size_to_readable(per_user)
    print('Usage per user:')
    for user, usage in per_user.items():
        print(f"{user:>15s}: {usage}")
    

def print_disclaimer():
    """
    The shell command `du` can take a long time
    to find all folder sizes if directory is large.
    """
    print("""\
This may take a while...

Calling the shell `du` command will 
take a while if the directory is large.

""")

def change_size_to_readable(data):
    """
    Convert the values of a one level dictionary of user: size
    to human readable sizes. i.e. 12000 bytes-> 12 K
    """

    def convert_size(size_bytes):
        if size_bytes == 0:
            return "0 B"
        size_name = ("B", "K", "M", "G", "T", "P", "E", "Z", "Y")
        i = int(math.floor(math.log(size_bytes, 1000)))
        p = math.pow(1000, i)
        s = round(size_bytes / p, 2)
        return "%s %s" % (s, size_name[i])

    return {k: convert_size(v) for k, v in data.items()}
    
def usage_per_user(data):
    owners = {}
    for folder, d in data.items():
        try:
            owner = d['owner']
        except KeyError:
            continue
        if owner not in owners:
            owners[owner] = 0
        else:
            value = d['size']
            owners[owner] += value

    return owners

def main():
    if '-h' in sys.argv:
        sys.exit('To exclude folders, run with -e:\n folder_sizes.py -e folder1 folder2')
    excluding = exclude_folders(sys.argv)
    print_disclaimer()
    dirs = {}
    dirs, total = human_readable_sizes(dirs, excluding)
    dirs = folder_owners(dirs)
    dirs = readable_to_comparable(dirs)
    lst = sort_dict_by_size(dirs)
    print_table(dirs, lst, total)

if __name__ == '__main__':
    main()
