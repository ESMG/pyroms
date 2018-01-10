#!/usr/bin/env python

import __future__
import json
import re
import copy
#import pdb; pdb.set_trace()

"""
Some tools for parsing the ROMS ocean.in file, saving it to JSON and/or
writing it back out to ocean.in format. The Python format is a list and
a dictionary. The dictionary has the settable variables such as 'Ngrids'
as keys, the values as values in string form. The list is more complex,
containing one element for comment lines and blank lines and up to four
elements for assignment lines: The variable name ('Ngrids'), the equals
string ("=" or "=="), the position of the first equals (13), and
optionally any trailing comment.

There are lingering issues with comments on continuation lines...
"""

class Line:
    """
    Helper class containing information about lines in an ocean.in file
    so that the file can be reconstructed in the same order, with the
    same comments and white space.
    """
    def __init__(self, text, ngrids=False, pos=0, comment="", extra=[]):
        """(Line, str, bool, int, str) -> NoneType
        Constructor for line with ROMS variable name.

        text -    String with ROMS variable name.
        ngrids -  One value per grid (True) or one for all grids (False)
        pos -     Position of first = in file.
        comment - String containing a final comment (if any)
        """
        self.text = text
        self.ngrids = ngrids
        self.pos = pos
        self.comment = comment
        self.extra = extra


class ocean_in:
    """
    Data structures needed for holding the information contained in a
    ROMS ocean.in file. It requires a dictionary of the settable values
    plus other things to keep the order and whitespace correct.
    """

    def do_lines(self, line_list):
# Private helper function for dealing with continuation lines
        first = line_list.pop(0)
        extra = []
        if len(first) > 1:
            extra.append(first[1])
        else:
            extra.append("")
        p = re.compile('([^=]+)(=+)([^=!]+)')
        try:
            eq_pos = first[0].find('=')
            m = p.match(first[0])
            key = m.group(1)
            eq = m.group(2)
            value = m.group(3)
            key = key.strip()
            value = value.lstrip()
        except:
            print(("trouble!", first))
            exit()
        key = key.strip()
        value_list = [value]
        for item in line_list:
            line = item[0]
            line = line.strip()
            value_list.append(line)
            if len(item) > 1:
                extra.append(item[1])
            else:
                extra.append("")
        ngrids = False
        if eq == "==": ngrids = True
        my_list = Line(key, ngrids=ngrids, pos=eq_pos, extra = extra)
        self.var_list.append(my_list)
        self.var_dict[key] = value_list

    def __init__(self, fname):
        """ (ocean_in, str) -> NoneType

        fname is a filename in the ROMS ocean.in format.
        """

        self.var_list = []
        self.var_dict = {}

        continuation = False
        temp_list = []

        fh = open(fname)
        for line in fh:
            # Handle blank lines and lines starting with !
            line = line.rstrip()
            if len(line) == 0 or line[0] == '!':
                line_list = Line(line)
                self.var_list.append(line_list)
                continue

            # Look for comments after the assignment and stash them away
            p = re.compile('(.*\S)( *!.*)')
            m = p.search(line)
            if m:
                line = m.group(1)
                after = m.group(2)
            else:
                after = None

            # Lots of extra work for continuations
            if line[-1] == '|' or line[-1] == '\\':
                continuation = True
                temp_list.append([line, after])
            elif continuation:
                # last in continue cluster
                continuation = False
                temp_list.append([line, after])
                self.do_lines(temp_list)
                temp_list = []
            else:
                # Solo lines
                eq_pos = line.find('=')
                line = line.strip()
                try:
                    p = re.compile('([^=]+)(=+)([^=!]+)')
                    m = p.match(line)
                    key = m.group(1)
                    eq = m.group(2)
                    value = m.group(3)
                except:
                    print(("trouble with line:", line))
                    exit()
                key = key.strip()
                value = [value.strip()]
                ngrids = False
                if eq == "==": ngrids = True
                line_list = Line(key, ngrids=ngrids, pos=eq_pos, comment=after)
                self.var_list.append(line_list)
                self.var_dict[key] = value

    def write_ocean_in(self, fname):
        """ (self, str)

        Given a filename and ROMS input data, create the ocean.in file.

        fname - name of file to write
        """
        try:
            fh = open(fname, 'w')
        except:
            print(('trouble opening file', fname))
            exit()

        for item in self.var_list:
            # Deal with non-var lines first
            if item.pos == 0:
                fh.write(item.text)
                fh.write("\n")
                continue

            varname = item.text
            eq = "="
            if item.ngrids: eq = "=="
            value = self.var_dict[varname]
            string = ' '.join(value)
            string = item.text + ' ' + eq + ' ' + string
            if item.comment:
                string += item.comment

            index = string.find('=')
            # We might need pos below.
            pos = item.pos
            if pos > index:
                string = ' '*(item.pos-index) + string

            # Look for continuation lines and split them up
            first = True
            p = re.compile('(\\\\|\|)')
            m = p.findall(string)
            if m:
                # All continuation chunks should have an "extra"
                # except after a merge
                extra = item.extra
                count = 0
                for match in m:
                    index = string.find(match)
                    part = string[0:index+len(match)]
                    if len(extra) > count and extra[count]:
                        part += extra[count]
                    count += 1
                    string = string[index+len(match):]

                    # Add space to front of trailing lines
                    if first:
                        first = False
                    else:
                        p = re.compile('([^ ])')
                        m = p.search(part)
                        index = part.find(m.group(1))
                        part = ' '*(pos-index + 3) + part

                    fh.write(part)
                    fh.write("\n")
                # Get the last continued line and pad it too
                p = re.compile('([^ ])')
                m = p.search(string)
                index = string.find(m.group(1))
                string = ' '*(pos-index + 3) + string
                if len(extra) > count and extra[count]:
                    string += extra[count]
                fh.write(string)
                fh.write("\n")
            else:
                fh.write(string)
                fh.write("\n")

    def write_json_dict(self, fname, indent=0):
        """(ocean_in, str, int) -> NoneType

        fname - filename to write to
        indent - set to non-zero for pretty-printed JSON file

        Save the ROMS.in dictionary to a JSON file. The default is compact while
        setting indent to a positive integer will create a more readable file.
        """
        fp = open(fname, 'w')
        if indent != 0:
            json.dump([self.var_dict], fp, indent=indent)
        else:
            json.dump([self.var_dict], fp)


    def merge_dicts(self, my_ocn_list):
        """(ocean_in, my_ocn_list)

        my_ocn_list - a list of dictionaries to add to self.

        Given a ROMS ocean_in object, merge the dictionaries from one or
        more other grids to create a multi-grid ocean_in.
        """
        # Copy here in case we're adding this object to itself
        t_dict = copy.deepcopy(self.var_dict)
        num_lists = len(my_ocn_list)
        for i in range(num_lists):
            new_dict = my_ocn_list[i].var_dict
            for it in self.var_list:
                var = it.text
                if it.ngrids:
                    try:
# If the first string ends in "\" or "|", it is likely to be
# a set of long filenames for which we need to add "\" before appending
                        first = t_dict[var][0]
                        if first[-1] == '|' or first[-1] == '\\':
                            t_dict[var][-1] += ' \\'
                        t_dict[var] = t_dict[var] + new_dict[var]
                    except:
                        print(("List", i, "is missing variable", var))

            t_dict['Ngrids'][0] = str(int(t_dict['Ngrids'][0]) + \
                          int(new_dict['Ngrids'][0]))
        self.var_dict = t_dict

def main():
#    fname = raw_input("Name of ocean.in file:")
    ocean_1 = ocean_in("ocean_foo.in")
#    jname = raw_input("Name of json file:")
    ocean_1.write_json_dict("foo.json")
    ocean_1.write_json_dict("foo2.json", indent=2)

    ocean_1.merge_dicts([ocean_1, ocean_1])
#    ocean_1.merge_dicts([ocean_1])
    ocean_1.write_ocean_in('ocean_foo2.in')

if __name__ == "__main__":
    main()

def test_one():
    ocean_2 = ocean_in("ocean_foo.in")
    assert ocean_2.var_dict['Lm'] == ['688']
    assert ocean_2.var_dict['Ngrids'] == ['1']
    assert len(ocean_2.var_dict['LBC(isTvar)']) == 6
    assert len(ocean_2.var_dict['LBC(isAice)']) == 1
    assert len(ocean_2.var_dict['FRCNAME']) == 24
    assert len(ocean_2.var_dict) == 17
    assert len(ocean_2.var_list) == 35

def test_two():
    import subprocess
    ocean_1 = ocean_in('ocean_foo.in')
    ocean_1.write_ocean_in('ocean_foo1.in')
    ocean_2 = ocean_in('ocean_foo1.in')
    ocean_2.write_ocean_in('ocean_foo2.in')
    foo = subprocess.check_output(['diff', 'ocean_foo1.in', 'ocean_foo2.in'])
    assert foo == ''
    subprocess.call('rm ocean_foo1.in ocean_foo2.in', shell=True)

    ocean_1.merge_dicts([ocean_2, ocean_1])
    assert ocean_1.var_dict['Lm'] == ['688', '688', '688']
    assert ocean_1.var_dict['Ngrids'] == ['3']
    assert len(ocean_1.var_dict['LBC(isTvar)']) == 18
    assert len(ocean_1.var_dict['LBC(isAice)']) == 3
    assert len(ocean_1.var_dict['FRCNAME']) == 72
    assert len(ocean_1.var_dict) == 17
    assert len(ocean_2.var_list) == 35
    ocean_1.write_ocean_in('ocean_foo3.in')
