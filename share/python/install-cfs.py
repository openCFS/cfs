#!/bin/env python

# the python-gitlab library: https://python-gitlab.readthedocs.io/en/stable/index.html
import gitlab

import os, sys

tipp = """
can be used in a crontab entry to regularly download the latest master build for Linux
e.g. to update every hour set (via 'croantab -e')
* */1 * * * /path/to/python /path/to/install-cfs.py -s /path/to/install/dir > /path/to/log/file.log
"""

def get_package_ids(version):
    ids = []
    packages = cfs.packages.list(get_all=True)
    for pkg in packages:
        if version in pkg.version:
            ids.append( pkg.id )
    return ids

# from https://stackoverflow.com/questions/22058048/hashing-a-file-in-python
import hashlib
def sha256sum(filename):
    with open(filename, 'rb', buffering=0) as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def get_package_files(pkg,file_name_component='Linux',latest=True):
    result = {}
    # list the files in the package
    files = pkg.package_files.list(get_all=True)
    for f in files:
        if file_name_component in f.file_name:
            result[f.id] = f
    if len(result.keys()) > 1 and latest: # try to determine the latest copy
        from datetime import datetime
        latest = datetime.strptime("2010-01-01T00:00:00.00Z","%Y-%m-%dT%H:%M:%S.%fZ") # an early date
        for key,f in result.items():
            created_at = datetime.strptime(f.created_at,"%Y-%m-%dT%H:%M:%S.%fZ") 
            if created_at > latest:
                latest = created_at
                latest_key = key
        result = { result[latest_key].file_name: result[latest_key] }
    elif len(result.keys()) == 1:
        for _,f in result.items():
            result = { f.file_name: f }
    # try to add the commit sha
    for key, f in result.items():
        try:
            pipelines = f.pipelines
            if not len(pipelines) == 1:
                print("Why do we have more than one pipeline for this file?",pipelines)
            sha = f.pipelines[0]['sha']
            if pkg.pipeline['sha'] == sha:
                print("setting file commit-sha =",sha)
                result[key].sha = sha
        except:
            print("could not set sha on file "+str(f))
    return result


def download_file(package_file,destination,attempts=3):
    if attempts >= 0:
        if os.path.exists(destination):
            print("file",destination,"exists")
            expected_sha256 = package_file.file_sha256
            if not expected_sha256 == sha256sum(destination):
                print("  but does not have expected sha256sum -> retry")
                download_file(package_file,destination,attempts=attempts-1)
            else:
                print("  sha256sum = %s -> OK"%(expected_sha256))
                return
        else:
            url = "https://gitlab.com/openCFS/cfs/-/package_files/%i/download"%(package_file.id)
            print("downloading",destination,"from", url)
            if attempts > 0 :
                import urllib.request
                urllib.request.urlretrieve(url, destination)
                download_file(package_file,destination,attempts=attempts-1)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog='install-cfs.py',
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    description='Installs CFS from GitLab package repository (from GitLab.com/openCFS/cfs)',
                    epilog=tipp)
    parser.add_argument('installdir',help='the parent directory to install the CFS package to') # positional argument
    parser.add_argument('-s', '--set-latest-symlink', action='store_true', help='set a symlink pointing to the latest release')
    parser.add_argument('-e', '--no-extract', action='store_true', default=False, help='do not extract the download')
    parser.add_argument('-a', '--install-all-files', action='store_true', default=False, help='install all files for a package version')
    parser.add_argument('-d', '--download-attempts', default=3, type=int, help='number of download attempts, default=3, 0 will only verify sha and print links')
    parser.add_argument('-p', '--package-version-pattern', default='master', type=str, help='selects a package version containing the given string, default=master')
    parser.add_argument('-f', '--file-name-pattern', default='Linux.tar.gz', type=str, help='selects files with names containing the given string, default=Linux.tar.gz')
    args = parser.parse_args()
    #print(args.installdir, args.set_latest_symlink, args.no_extract, args.package_version_pattern)

    install_dir = args.installdir
    if not os.path.isdir( install_dir ):
        sys.exit("given installdir %s does not exist."%install_dir)

    # anonymous read-only access for public resources (GitLab.com)
    gl = gitlab.Gitlab()

    # get the CFS project (get the ID from GitLab.com)
    cfs = gl.projects.get(12930334)

    # get all packages for the project
    packages = cfs.packages.list(get_all=True)

    # select the packages based on pattern
    version = args.package_version_pattern
    package_ids = get_package_ids(version)
    if len(package_ids) == 0 :
        sys.exit("Did not find a package version matching '"+version+"' found: "+str([pkg.version for pkg in packages]))

    for package_id in package_ids:
        pkg = cfs.packages.get(package_id) # 21499193 = master
        print("Found package %s:%s"%(pkg.name,pkg.version))

        # list the files in the package
        files = pkg.package_files.list(get_all=True)
        
        file_dict = get_package_files(pkg,args.file_name_pattern,not args.install_all_files)

        for key,f in file_dict.items():
            # if we have the master build, rename to a unique name
            if type(key) == int or "master" in key:
                name = f.file_name
                try:
                    newname = name.replace("master",f.sha)
                except:
                    newname = name.replace("master",str(f.id))
                destination = os.path.join(install_dir,newname)
            else:
                destination = os.path.join(install_dir,key)
            download_file(file_dict[key],destination,args.download_attempts)

            if (not args.no_extract):
                import tarfile
                with tarfile.open(destination) as archive:
                    print("Attempting to extract",destination,"to",install_dir)
                    for member in archive.getmembers():
                        if not os.path.exists( os.path.join(install_dir,member.name) ):
                            print(" ",member.name)
                            archive.extract(member)
                        else:
                            print(" ",member.name,"already exists")
    
    if args.set_latest_symlink:
        if len(package_ids) == 1 and len(file_dict.keys()) == 1:
            target = os.path.basename(destination).split('.')[0]
            if os.path.islink( os.path.join(install_dir,"latest") ) and os.readlink( os.path.join(install_dir,"latest") ) == target :
                print("Symlink is already correct")
            else :
                print("Setting symlink '%s/latest' -> '%s'"%(install_dir,target))
                os.symlink( target,  os.path.join(install_dir,"new-link") )
                os.rename( os.path.join(install_dir,"new-link"), os.path.join(install_dir,"latest")  )
        else :
            sys.exit("Cannot set latest symlink if there are more than one package or file")
