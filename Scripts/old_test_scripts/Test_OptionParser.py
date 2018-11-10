#! /usr/bin/env python
"""TESt
"""
#python test.py one two d cual 546
#python test.py what -s 0


#import sys
from optparse import OptionParser

def foo_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

def make_list(option, opt_str, value, parser):
    setattr(parser.values, option.dest, map(int,value.replace('[','').replace(']','').split(',')))

parser = OptionParser()
parser.add_option("-s","--substation",type="int",default=1,help="If using HBA substations (HBA0 and HBA1) then 1. Otherwise use 0.")
parser.add_option("-f","--flag_antenna",type="str",action='callback',callback=make_list,help="Bad antennas to flag.")
parser.add_option('-d', '--foo',
                  type='string',
                  action='callback',
                  callback=foo_callback)

(options, args) = parser.parse_args()

print 'options',options
print 'args',args


if not parser.get_prog_name()=="test_OptionParser.py":
    #   Program was run from within python
    substation = 1
    flag_antenna = None
else:
    substation = options.substation
    flag_antenna = options.flag_antenna


print 'substation', substation
print 'flag_antenna', flag_antenna

#for item in sys.argv:
#    print item, type(item)



# def record_foo_seen(option, opt_str, value, parser):
#     parser.saw_foo = True
#
# parser.add_option("--foo", action="callback", callback=record_foo_seen)
#
#
# def check_moon(option, opt_str, value, parser):
#     if is_moon_full():
#         raise OptionValueError("%s option invalid when moon is full"
#                                % opt_str)
#     setattr(parser.values, option.dest, 1)
# [...]
# parser.add_option("--foo",
#                   action="callback", callback=check_moon, dest="foo")
#
#
# def store_value(option, opt_str, value, parser):
#     setattr(parser.values, option.dest, value)
# [...]
# parser.add_option("--foo",
#                   action="callback", callback=store_value,
#                   type="int", nargs=3, dest="foo")


