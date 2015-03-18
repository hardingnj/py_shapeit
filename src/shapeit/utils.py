__author__ = 'Nicholas Harding'
import re
import sh
import hashlib


def parse_command(parameters):

    cl = re.compile('^-')
    last_was_key = False
    key = None
    parameters = [str(x) for x in parameters]

    command_dict = {}
    for value in parameters:
        value = str(value)
        if cl.match(value):
            key = value
            command_dict[key] = True
            last_was_key = True
        else:
            if last_was_key:
                command_dict[key] = value
                last_was_key = False
            else:
                command_dict[key] = command_dict[key] + ';' + value

    return command_dict


def create_sh_script(filename, commands=None, outfile=None):

    """
    Creates an executable sh script that is prepopulated with useful defaults
    :param filename: name of file to write
    :param commands: list of commands as strings that will be executed
    :param outfile: optional, creates an 'ok' version of this file to validate
        completion.
    :return: None
    """

    # create script file
    if outfile is None:
        outfile = filename

    touch_cmd = "touch {FILE}"

    script = open(filename, 'w')
    script.write("#! /bin/bash" + "\n")
    script.write("set -e" + "\n")
    script.write("set -o pipefail" + "\n")
    for cmd in commands:
        script.write(cmd + "\n")
    script.write(touch_cmd.format(FILE=outfile + '.ok') + "\n")
    sh.chmod('+x', filename)
    script.close()


def calc_regions(size, nbins=20, overlap=0):

    """
    A utility function to assist in working out a regions that overlap. Useful
    for ligateHaplotypes etc.
    :param size: size of total region
    :param nbins: number of regions
    :param overlap: stagger
    :return: a list of tuples each representing a start and an end.
    """
    if overlap is None:
        overlap = size/(nbins*20)
    approx_size = 1 + size/nbins + overlap - (overlap/nbins)
    print approx_size, overlap
    regions = []

    for i in range(nbins):
        start = i*(approx_size-overlap) + 1
        stop = start + approx_size
        regions.append((start, stop))
    return regions


def md5_for_file(f, block_size=2**20):
    fh = open(f, 'rb')
    md5 = hashlib.md5()
    while True:
        data = fh.read(block_size)
        if not data:
            break
        md5.update(data)
    return md5.hexdigest()