#!/usr/bin/env python

from pyroms_toolbox import ocean_in

def main():
#    fname = raw_input("Name of ocean.in file:")
    ocean_1 = ocean_in("ocean_foo.in")
#    jname = raw_input("Name of json file:")
    ocean_1.write_json_dict("foo.json")
    ocean_1.write_json_dict("foo2.json", indent=2)

    ocean_1.merge_dicts([ocean_1, ocean_1])
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
