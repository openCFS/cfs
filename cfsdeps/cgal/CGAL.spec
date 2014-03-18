%define boost_version 1.32

Name:           cgal-lse-devel
Version:        3.3
Release:        5%{?dist}
Summary:        Computational Geometry Algorithms Library

Group:          System Environment/Libraries
License:        QPL/LGPL
URL:            http://www.cgal.org/
Packager:	Simon Triebenbacher <simon.triebenbacher@gmx.net>

Source0:        ftp://ftp.mpi-sb.mpg.de/pub/outgoing/CGAL/CGAL-%{version}.tar.gz
Source1:	distro_sh.tgz
Source2:        patched_headers.tgz
Source10:       CGAL-README.Fedora

Patch1:         CGAL-install_cgal-SUPPORT_REQUIRED.patch
Patch2:         CGAL-3.3-build-library.patch
Patch4:         CGAL-3.2.1-install_cgal-no_versions_in_compiler_config.h.patch

BuildRoot:      %{_tmppath}/%{name}-%{version}-%{release}-buildroot
Prefix:		/opt/lse


# Required packages.
BuildRequires: gmp-devel
BuildRequires: mpfr-devel
BuildRequires: boost-devel >= %boost_version
# BuildRequires: qt3-devel >= 3.0
BuildRequires: zlib-devel
BuildRequires: blas lapack

# for chmod in prep:
BuildRequires: coreutils

# CGAL-libs is renamed to CGAL.
Obsoletes:     %{name}-libs < %{version}-%{release}
Provides:      %{name}-libs = %{version}-%{release}

%description
Libraries for CGAL applications.
CGAL is a collaborative effort of several sites in Europe and
Israel. The goal is to make the most important of the solutions and
methods developed in computational geometry available to users in
industry and academia in a C++ library. The goal is to provide easy
access to useful, reliable geometric algorithms.


%package demos-source
Group:          Documentation
Summary:        Examples and demos of CGAL algorithms
Requires:       %{name} = %{version}-%{release}
Obsoletes:      %{name}-demo < %{version}-%{release}
Provides:       %{name}-demo = %{version}-%{release}
%description demos-source
The %{name}-demos-source package provides the sources of examples and demos of
CGAL algorithms.


%prep
rm -rf $RPM_BUILD_ROOT

%setup -q -n CGAL-%{version}
tar xvzf ../../SOURCES/distro_sh.tgz
tar xvzf ../../SOURCES/patched_headers.tgz

%patch1 -p0 -b .support-required.bak
%patch2 -p1 -b .build-library.bak

%patch4 -p1 -b .no_versions.bak

chmod a-x examples/Nef_3/handling_double_coordinates.cin



# fix end-of-lines of several files
sed -i 's/\r//' \
    examples/Surface_mesh_parameterization/data/mask_cone.off \
    examples/Boolean_set_operations_2/test.dxf

for f in demo/Straight_skeleton_2/data/vertex_event_9.poly \
         demo/Straight_skeleton_2/data/vertex_event_0.poly;
do
  [ -r $f ] && sed -i 's/\r//' $f;
done

# Install README.Fedora here, to include it in %doc
install -m 644 %{SOURCE10} ./README.Fedora

mv *.h include/CGAL


%build
# Check if the C compiler has been set externally
if [ "xxx$CC" == "xxx" ]; then
  export CC=gcc
fi

# Check if the C++ compiler has been set externally
if [ "xxx$CXX" == "xxx" ]; then
  export CXX=g++
fi

# Check if the Fortran compiler has been set externally
if [ "xxx$FC" == "xxx" ]; then
  export FC=gfortran
fi

echo "Compiled for:" > compiler.txt
echo "---------------------------------------------" >> compiler.txt
echo `./distro.sh -h` >> compiler.txt
echo >> compiler.txt
echo "C compiler:" >> compiler.txt
echo "---------------------------------------------" >> compiler.txt
$CC --version >> compiler.txt
echo >> compiler.txt
echo "C++ compiler:" >> compiler.txt
echo "---------------------------------------------" >> compiler.txt
$CXX --version >> compiler.txt
echo >> compiler.txt
echo "Fortran compiler:" >> compiler.txt
echo "---------------------------------------------" >> compiler.txt
$FC --version >> compiler.txt

if [ "$FC" == "gfortran" ]; then
  FORTRAN_FLAGS="-lgfortran";
fi

# source %{_sysconfdir}/profile.d/qt3.sh

./install_cgal -ni $CXX --CUSTOM_CXXFLAGS "$RPM_OPT_FLAGS" \
               --CUSTOM_LDFLAGS "$FORTRAN_FLAGS" \
               --without-autofind \
               --with-ZLIB \
               --with-BOOST \
               --with-BOOST_PROGRAM_OPTIONS \
               --with-X11 \
               --with-GMP \
               --with-GMPXX \
               --with-MPFR \
               --with-QT3MT \
               --with-REFBLASSHARED \
               --with-LAPACK \
               --with-OPENGL \
               --disable-shared
#               --prefix $RPM_BUILD_ROOT%{prefix}


%install
rm -rf %{buildroot}

# Install headers
mkdir -p %{buildroot}%{prefix}/include
cp -a include/* %{buildroot}%{prefix}/include
rm -rf %{buildroot}%{prefix}/include/CGAL/config/msvc*
mv %{buildroot}%{prefix}/include/CGAL/config/*/CGAL/compiler_config.h %{buildroot}%{prefix}/include/CGAL/
rm -rf %{buildroot}%{prefix}/include/CGAL/config

# Install scripts (only those prefixed with "cgal_").
mkdir -p %{buildroot}%{prefix}/bin
cp -a scripts/cgal_* %{buildroot}%{prefix}/bin

# Install libraries
mkdir -p %{buildroot}%{prefix}/%{_lib}
cp -a lib/*/lib* %{buildroot}%{prefix}/%{_lib}

# Install makefile:
mkdir -p %{buildroot}%{prefix}/share/CGAL
touch -r make %{buildroot}%{prefix}/share/CGAL
cp -p make/makefile_* %{buildroot}%{prefix}/share/CGAL/cgal.mk

# Install demos and examples
mkdir -p %{buildroot}%{prefix}/share/CGAL/
touch -r demo %{buildroot}%{prefix}/share/CGAL/
cp -a demo %{buildroot}%{prefix}/share/CGAL/demo
cp -a examples %{buildroot}%{prefix}/share/CGAL/examples

# Modify makefile
cat > makefile.sed <<'EOF'
s,CGAL_INCL_DIR *=.*,CGAL_INCL_DIR = %{prefix}/include,;
s,CGAL_LIB_DIR *=.*,CGAL_LIB_DIR = %{prefix}/%{_lib},;
/CUSTOM_CXXFLAGS/ s/-O2 //;
/CUSTOM_CXXFLAGS/ s/-g //;
/CGAL_INCL_DIR/ s,/CGAL/config/.*,,;
s,/$(CGAL_OS_COMPILER),,g;
/-I.*CGAL_INCL_CONF_DIR/ d
EOF

sed -i -f makefile.sed %{buildroot}%{prefix}/share/CGAL/cgal.mk

# check if the sed script above has worked:
grep -q %{_builddir} %{buildroot}%{prefix}/share/CGAL/cgal.mk && false
grep -q %{buildroot} %{buildroot}%{prefix}/share/CGAL/cgal.mk && false
grep -q CGAL/config %{buildroot}%{prefix}/share/CGAL/cgal.mk && false
grep -q -E 'CUSTOM_CXXFLAGS.*(-O2|-g)' %{buildroot}%{prefix}/share/CGAL/cgal.mk && false

# Remove -L and -R flags from the makefile
cat > makefile-noprefix.sed <<'EOF'
/'-L$(CGAL_LIB_DIR)'/ d;
/-R$(CGAL_LIB_DIR)/ d;
/'-I$(CGAL_INCL_DIR)'/ d;
EOF

sed -i -f makefile-noprefix.sed  %{buildroot}%{prefix}/share/CGAL/cgal.mk

# check that the sed script has worked
grep -q -E -- '-[LI]\$' %{buildroot}%{prefix}/share/CGAL/cgal.mk && false
grep -q -E -- '-R' %{buildroot}%{prefix}/share/CGAL/cgal.mk && false


# Create %{_sysconfdir}/profile.d/ scripts
# cd %{buildroot}
# mkdir -p %{prefix}.%{_sysconfdir}/profile.d
# cat > %{prefix}.%{_sysconfdir}/profile.d/cgal.sh <<EOF
# if [ -z "\$CGAL_MAKEFILE" ] ; then
#   CGAL_MAKEFILE="%{prefix}%{_datadir}/CGAL/cgal.mk"
#   export CGAL_MAKEFILE
# fi
# EOF

# cat > %{prefix}.%{_sysconfdir}/profile.d/cgal.csh <<EOF
# if ( ! \$?CGAL_MAKEFILE ) then
#  setenv CGAL_MAKEFILE "%{prefix}%{_datadir}/CGAL/cgal.mk"
# endif
# EOF
# chmod 644 %{prefix}.%{_sysconfdir}/profile.d/cgal.*sh

mkdir -p %{buildroot}/%{prefix}/share/doc/packages/%{name}
cp compiler.txt %{buildroot}/%{prefix}/share/doc/packages/%{name}
cp AUTHORS %{buildroot}/%{prefix}/share/doc/packages/%{name}
cp LICENSE %{buildroot}/%{prefix}/share/doc/packages/%{name}
cp LICENSE.FREE_USE %{buildroot}/%{prefix}/share/doc/packages/%{name}
cp LICENSE.LGPL %{buildroot}/%{prefix}/share/doc/packages/%{name}
cp LICENSE.QPL %{buildroot}/%{prefix}/share/doc/packages/%{name}
cp CHANGES %{buildroot}/%{prefix}/share/doc/packages/%{name}
cp README.Fedora %{buildroot}/%{prefix}/share/doc/packages/%{name}

%clean
rm -rf %{buildroot}


%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig


%files
%defattr(-,root,root,-)
%{prefix}/%{_lib}/libCGAL*
%{prefix}/share/doc/packages/%{name}
%{prefix}/include/CGAL
%dir %{prefix}/share/CGAL
%{prefix}/share/CGAL/cgal.mk
%{prefix}/bin/*
%exclude %{prefix}/bin/cgal_make_macosx_app
# %config(noreplace) %{prefix}%{_sysconfdir}/profile.d/cgal.*


%files demos-source
%defattr(-,root,root,-)
%dir %{prefix}/share/CGAL
%{prefix}/share/CGAL/demo
%{prefix}/share/CGAL/examples
%exclude %{prefix}/share/CGAL/*/*/*.vcproj
%exclude %{prefix}/share/CGAL/*/*/skip_vcproj_auto_generation


%changelog
* Mon Jun 18 2007 Simon Triebenbacher <simon.triebenbacher@gmx.net> - 3.3-5%{?dist}
- CGAL gets installed in /opt/lse.
- CGAL builds on SUSE.

* Thu Jun  7 2007 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-4%{?dist}
- Move the makefile back to %%{_datadir}/CGAL, and rename it cgal.mk (sync
  with Debian package). That file is not a config file, but just an example
  .mk file that can be copied and adapted by users.
- Fix the %%{_sysconfdir}/profile.d/cgal.* files (the csh one was buggy).
- CGAL-devel now requires all its dependancies.

* Sat Jun  2 2007 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-2%{?dist}
- Officiel CGAL-3.3 release
- Skip file named "skip_vcproj_auto_generation"

* Wed May 30 2007 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-0.1.RC1%{?dist}
- New upstream version: 3.3-RC1
- Obsolete patches CGAL-3.2.1-build-libCGALQt-shared.patch,
                   CGAL-3.2.1-build-no-static-lib.patch,
                   CGAL-3.2.1-config.h-endianness_detection.patch.
  These patchs have been merged and adapted by upstream.
- New option --disable-static
- Shipped OpenNL and CORE have been renamed by upstream:
    - %%{_includedir}/OpenNL is now %{_includedir}/CGAL/OpenNL
    - %%{_includedir}/CORE is now %{_includedir}/CGAL/CORE
    - libCORE has been rename libCGALcore++
  Reasons:
    - CGAL/OpenNL is a special version of OpenNL, rewritten for CGAL 
      in C++ by the OpenNL author,
    - CGAL/CORE is a fork of CORE-1.7. CORE-1.7 is no longer maintained by 
      its authors, and CORE-2.0 is awaited since 2004.
  In previous releases of this package, CORE was excluded from the package, 
  because %%{_includedir}/CORE/ was a name too generic (see comment #8 of
  #bug #199168). Since the headers of CORE have been moved to
  %%{_includedir}/CGAL/CORE, CORE is now shipped with CGAL.
- move %%{_datadir}/CGAL/make/makefile to %%{_sysconfdir}/CGAL/makefile
(because it is a config file).

* Thu Nov  2 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2.1-19
- Fix CGAL-3.2.1-build-libCGALQt-shared.patch (bug #213675)

* Fri Sep 15 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2.1-18
- Move LICENSE.OPENNL to %%doc CGAL-devel (bug #206575).

* Mon Sep 11 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2.1-17
- libCGALQt.so needs -lGL
- remove %%{_bindir}/cgal_make_macosx_app

* Sat Sep  11 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2.1-16
- Remove CORE support. Its acceptance in Fedora is controversial (bug #199168).
- Exclude .vcproj files from the -demos-source subpackage.
- Added a patch to build *shared* library libCGALQt.
- Added a patch to avoid building static libraries.
- Fixed the License: tag.

* Thu Aug 17 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2.1-15
- Change the permissions of %%{_sysconfdir}/profile.d/cgal.*sh
- Remove the meta package CGAL. CGAL-libs is renamed CGAL.
- Added two patchs:
  - CGAL-3.2.1-config.h-endianness_detection.patch which is an upstream patch
    to fix the endianness detection, so that is is no longer hard-coded in
    <CGAL/compiler_config.h>,
  - CGAL-3.2.1-install_cgal-no_versions_in_compiler_config.h.patch that
    removes hard-coded versions in <CGAL/compiler_config.h>.

* Wed Aug 16 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2.1-14
- Simplified spec file, for Fedora Extras.

* Mon Jul 17 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2.1-13
- Change CGAL-README.Fedora, now that Installation.pdf is no longer in the
tarball.

* Mon Jul 17 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2.1-12
- Remove unneeded  -R/-L/-I flags from %%{_datadir}/CGAL/make/makefile

* Mon Jul 17 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2.1-11
- Fix the soversion.
- Fix %%{cgal_prefix} stuff!!
- Quote 'EOF', so that the lines are not expanded by the shell.

* Tue Jul  4 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2.1-10
- Fix makefile.sed so that %%{buildroot} does not appear in 
  %%{_datadir}/CGAL/make/makefile.

* Sun Jul  2 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2.1-9
- Remove Obsoletes: in the meta-package CGAL.

* Sun Jul  2 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2.1-8
- Fix the localisation of demo and examples.

* Sun Jul  2 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2.1-6
- Set Requires, in sub-packages.

* Sun Jul  2 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2.1-5
- CGAL-3.2.1
- Sub-package "demo" is now named "demos-source" (Fedora guidelines).
- Fix some rpmlint warnings
- Added README.Fedora, to explain why the documentation is not shipped, and how CGAL is divided in sub-packages.


* Sat Jul  1 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2-4
- Use %%{_datadir}/CGAL instead of %%{_datadir}/%%{name}-%%{version}
- Fix %%{_datadir}/CGAL/makefile, with a sed script.
- Added a new option %%set_prefix (see top of spec file).

* Sat Jul  1 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2-3
- Use less "*" in %%files, to avoid futur surprises.
- Remove %%{_sysconfdir}/profile.d/cgal.* from %%files if %%cgal_prefix is not empty.
- Fix %%build_doc=0 when %%fedora is set. New option macro: %%force_build_doc.

* Fri Jun 30 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2-2
- Fix some end-of-lines in %%prep, to please rpmlint.

* Mon May 22 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.2-1
- Remove README from %%doc file: it describes the tarball layout.
- Updated to CGAL-3.2.
- Added examples in the -demo subpackage.
- Cleaning up, to follow Fedora Guidelines.
- The -doc subpackage cannot be build on Fedora (no license).
- Add ldconfig back.
- No prefix.

* Fri Apr 28 2006 Laurent Rineau <laurent.rineau__fc_extra@normalesup.org> - 3.2-0.447
- Update to CGAL-3.2-447.

* Fri Apr 21 2006 Laurent Rineau <laurent.rineau__fc_extra@normalesup.org> - 3.2-0.440
- Updated to CGAL-3.2-I-440.

* Wed Apr 19 2006 Laurent Rineau <laurent.rineau__fc_extra@normalesup.org> - 3.2-0.438
- Added a patch to install_cgal, to require support for BOOST, BOOST_PROGRAM_OPTIONS, X11, GMP, MPFR, GMPXX, CORE, ZLIB, and QT.
- Move scripts to %%{_bindir}
- %%{_libdir}/CGAL-I now belong to CGAL and CGAL-devel, so that it disappears when the packages are removed.

* Wed Apr 12 2006 Laurent Rineau <laurent.rineau__fc_extra@normalesup.org> - 3.2-0.431
- Updated to CGAL-3.2-I-431.
- Remove the use of ldconfig.
- Changed my email address.
- No longer need for patch0.
- Pass of rpmlint.
- Remove unneeded Requires: tags (rpm find them itself).
- Change the release tag.
- Added comments at the beginning of the file.
- Added custom ld flags, on 64bits archs (so that X11 is detected).

* Tue Apr 11 2006 Laurent Rineau <laurent.rineau__fc_extra@normalesup.org>
- Removed -g and -O2 from CUSTOM_CXXFLAGS, in the makefile only.
  They are kept during the compilation of libraries.
- Added zlib in dependencies.
- Added a patch to test_ZLIB.C, until it is merged upstream.

* Fri Mar 31 2006 Naceur MESKINI <nmeskini@sophia.inria.fr>
- adding a test in the setup section.

* Mon Mar 13 2006 Naceur MESKINI <nmeskini@sophia.inria.fr>
- delete the patch that fixes the perl path.
- add build_doc and build_demo flags.

* Fri Mar 10 2006 Naceur MESKINI <nmeskini@sophia.inria.fr>
- Adding new sub-packages doc(pdf&html) and demo.
- Add internal_release flag. 

* Thu Mar 09 2006 Naceur MESKINI <nmeskini@sophia.inria.fr>
- Cleanup a specfile.

