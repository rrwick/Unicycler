"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains a Unicycler class for writing output to both stdout and a log file.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import sys
import os
import datetime
import re
import shutil
import textwrap
import subprocess


class Log(object):

    def __init__(self, log_filename=None, stdout_verbosity_level=1, log_file_verbosity_level=None):
        self.log_filename = log_filename

        # Determine if the terminal supports colours or not.
        try:
            self.colours = int(subprocess.check_output(['tput', 'colors']).decode().strip())
        except (ValueError, subprocess.CalledProcessError, FileNotFoundError, AttributeError):
            self.colours = 1

        # There are two verbosity levels: one for stdout and one for the log file. They are the
        # same, except that the log file verbosity level is never 0.
        self.stdout_verbosity_level = stdout_verbosity_level
        if not log_file_verbosity_level:
            self.log_file_verbosity_level = stdout_verbosity_level
        else:
            self.log_file_verbosity_level = log_file_verbosity_level
        self.log_file_verbosity_level = max(1, self.log_file_verbosity_level)

        if self.log_filename:
            log_file_exists = os.path.isfile(self.log_filename)
            self.log_file = open(self.log_filename, 'at', 1, encoding='utf8')  # line buffering

            # If the log file already exists, we pad out a bit of space before appending to it.
            if log_file_exists:
                self.log_file.write('\n\n\n\n\n\n\n\n\n\n\n\n')
        else:
            self.log_file = None

    def __del__(self):
        if self.log_file and not self.log_file.closed:
            self.log_file.close()


# This is the one and only instance of the Log class.
logger = Log()


def log(text, verbosity=1, stderr=False, end='\n', print_to_screen=True, write_to_log_file=True):
    text_no_formatting = remove_formatting(text)

    # The text is printed to the screen with ANSI formatting, if supported. If there are only 8
    # colours available, then remove the 'dim' format which doesn't work.
    if stderr or (verbosity <= logger.stdout_verbosity_level and print_to_screen):
        if logger.colours <= 1:
            text = text_no_formatting
        elif logger.colours <= 8:
            text = remove_dim_formatting(text)
        if stderr:
            print(text, file=sys.stderr, end=end, flush=True)
        else:
            print(text, end=end, flush=True)

    # The text is written to file without ANSI formatting.
    if logger.log_file and verbosity <= logger.log_file_verbosity_level and write_to_log_file:
        logger.log_file.write(text_no_formatting)
        logger.log_file.write('\n')


def log_section_header(message, verbosity=1, single_newline=False):
    """
    Logs a section header. Also underlines the header using a row of dashes to the log file
    (because log files don't have ANSI formatting).
    """
    if single_newline:
        log('', verbosity)
    else:
        log('\n', verbosity)

    time = get_timestamp()
    time_str = '(' + time + ')'
    if logger.colours > 8:
        time_str = dim(time_str)
    log(bold_yellow_underline(message) + ' ' + time_str, verbosity)
    log('-' * (len(message) + 3 + len(time)), verbosity, print_to_screen=False)


def log_progress_line(completed, total, base_pairs=None, end_newline=False):
    """
    Logs a progress line using a carriage return to overwrite the previous progress line. Only the
    final progress line will be written to the log file.
    """
    progress_str = int_to_str(completed) + ' / ' + int_to_str(total)
    if total > 0:
        percent = 100.0 * completed / total
    else:
        percent = 0.0
    progress_str += ' (' + '%.1f' % percent + '%)'
    if base_pairs is not None:
        progress_str += ' - ' + int_to_str(base_pairs) + ' bp'

    end_char = '\n' if end_newline else ''
    log('\r' + progress_str, end=end_char, write_to_log_file=False)
    if end_newline:
        log(progress_str, print_to_screen=False)


def log_explanation(text, verbosity=1, print_to_screen=True, write_to_log_file=True,
                    extra_empty_lines_after=1, indent_size=4):
    """
    This function writes explanatory text to the screen. It is wrapped to the terminal width for
    stdout but not wrapped for the log file.
    """
    text = ' ' * indent_size + text
    if print_to_screen:
        terminal_width = shutil.get_terminal_size().columns
        for line in textwrap.wrap(text, width=terminal_width - 1):
            if logger.colours > 8:
                formatted_text = dim(line)
            else:
                formatted_text = line
            log(formatted_text, verbosity=verbosity, print_to_screen=True, write_to_log_file=False)
    if write_to_log_file:
        log(text, verbosity=verbosity, print_to_screen=False, write_to_log_file=True)

    for _ in range(extra_empty_lines_after):
        log('', verbosity=verbosity, print_to_screen=print_to_screen,
            write_to_log_file=write_to_log_file)


def log_number_list(numbers, verbosity=1, print_to_screen=True, write_to_log_file=True,
                    indent_size=4):
    """
    Some lists of numbers are long, so this function makes them wrap nicely when displayed on
    stdout.
    """
    text = ' ' * indent_size + ', '.join(str(x) for x in numbers)
    if print_to_screen:
        terminal_width = shutil.get_terminal_size().columns
        for line in textwrap.wrap(text, width=terminal_width - 1):
            log(line, verbosity=verbosity, print_to_screen=True, write_to_log_file=False)
    if write_to_log_file:
        log(text, verbosity=verbosity, print_to_screen=False, write_to_log_file=True)


def int_to_str(num, max_num=0):
    if num is None:
        num_str = 'n/a'
    else:
        num_str = '{:,}'.format(num)
    max_str = '{:,}'.format(int(max_num))
    return num_str.rjust(len(max_str))


def get_timestamp():
    return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
YELLOW = '\033[93m'
DIM = '\033[2m'


def bold_yellow_underline(text):
    return YELLOW + BOLD + UNDERLINE + text + END_FORMATTING


def dim(text):
    return DIM + text + END_FORMATTING


def remove_formatting(text):
    return re.sub('\033.*?m', '', text)


def remove_dim_formatting(text):
    return re.sub('\033\\[2m', '', text)
