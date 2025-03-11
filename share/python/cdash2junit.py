#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="Get information from CDash and write a JUnit file")
parser.add_argument('buildid', type=int, nargs=1,
                    help='the build id of the build for which to collext test details')
parser.add_argument('--cdash-url', default="http://movm.mi.uni-erlangen.de/cdash/",
                    help='the url of the cdash server')

args = parser.parse_args()
cdash_url = args.cdash_url
if not cdash_url[-1] == '/': # add one
     cdash_url += '/'
buildid = args.buildid[0]

# get json from CDASH

# build the URL as seen on the JSON link at the bottom of the CDASH page
url = cdash_url + "api/v1/viewTest.php?buildid=%i"%(buildid)

print("getting JSON from CDash url:",url)
from urllib.request import urlopen
import json
# store the response of URL
response = urlopen(url)
# storing the JSON response
cdash_json = json.loads(response.read())

# check if there were tests found
if len(cdash_json['tests']) == 0:
  raise Exception("No tests found in JSON response: cdash down or wrong buildid?")

# write the JUnit file
junit_out = '%i.xml'%buildid
print("writing JUnit file:",junit_out)
f = open(junit_out, 'w')
def write(x,f=f):
    #print(x.strip())
    f.write(x)
buildname = cdash_json['build']['name']
write('<testsuite name="%i">\n'%(buildid)) # name probably not used
for test in cdash_json['tests']:
    # file = linked in gitlab against the source (should be the file tested)
    write('  <testcase id="%i" name="%s" classname="%s" file="%s" time="%.3f">\n'%(test['id'],test['name'],buildname,test['name'],test['execTimeFull']))
    if test['status']=='Failed':
        write('    <failure>'+cdash_url+test['detailsLink'].replace('&','&#38;')+'</failure>\n')
        test_url = cdash_url + "api/v1/testDetails.php?test=%i&build=%i"%(test['id'],buildid)
        print("  getting failed test details JSON from CDash url:",test_url)
        test_response = urlopen( test_url )
        test_json = json.loads(test_response.read())
        output = test_json['test']['command']+'\n\n'+test_json['test']['output'].strip()
    else:
        output = 'gitlab only displays output for failed tests, so we do not need to query in other cases'
    write('    <system-out>%s</system-out>\n'%(output))
    write('  </testcase>\n')
write('</testsuite>')
f.close()
