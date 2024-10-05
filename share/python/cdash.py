#!/usr/bin/env python

# this tool searches the openCFS server for reports on a given test case

import sys
import pandas
import datetime
import argparse

today = datetime.datetime.today()
#url = 'http://movm.mi.uni-erlangen.de/cdash/api/v1/testSummary.php?date=2024-09-16&name=Coupledfield_LinFlowMechSmooth_CompressedChamber2D&project=1&export=csv'
server = 'http://movm.mi.uni-erlangen.de/cdash/api/v1/testSummary.php?'
url = server + 'project=1'
      
txt =  'searches the openCFS cdash server for a specific test case and reports when it failed'
parser = argparse.ArgumentParser(description=txt)
parser.add_argument('test', help="e.g. 'Coupledfield_LinFlowMechSmooth_CompressedChamber2D'")
parser.add_argument('--end', help="end date YYYY-MM-DD, default is today",default=today.strftime('%Y-%m-%d'))
parser.add_argument('--days', help='how many days to go back. Default is 60',type=int, default=60)
parser.add_argument('--start', help="alternative to --days in format YYYY-MM-DD")
args = parser.parse_args()

end   = datetime.datetime.fromisoformat(args.end) if args.end else today
start = datetime.datetime.fromisoformat(args.start) if args.start else end - datetime.timedelta(days=args.days)

url += '&name=' + args.test


errors = 0 # we count the errors (days we cannot get something from cdash)
good = 0 # count 'Passed'
bad  = 0 # count 'Failed'

days = (end-start).days
for d in range(1,days):
  date = (end-datetime.timedelta(days=d)).strftime('%Y-%m-%d')
  du = '&date=' + date 
  df = pandas.read_csv(url + du + '&export=csv',delimiter=',')
  # df has columns without data for wrong test or date
  status = list(df['Status'])
  if len(status) == 0:
    errors += 1
    if errors > 4:
      print('we got',errrors,'from cdash - skip the remeaining',days-d,'requests')
      break
  else:
    good += status.count('Passed')
    bad  += status.count('Failed')
    
    if status.count('Failed') > 0:
      print('-' + str(d),'days:',date,':',url + du)
      site = list(df['Site'])
      for t in range(len(status)):
        if status[t] != 'Passed':
           print(' -> ',site[t])  
      
print("checked test '",args.test,'" for',days,'days:')
print(bad,'failed tests')
print(good,'passed tests')
print(errors,'errors reported from',server +  'project=1')

  