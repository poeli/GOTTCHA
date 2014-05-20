#!/bin/bash

export PATH=ext/bin/:$PATH;
export PERL5LIB=ext/lib/:$PERL5LIB;
mkdir -p bin;

done_message () {
  if [ $? == 0 ]; then
      if [ "x$1" != "x" ]; then
          echo $1;
      else
          echo "done.";
      fi
  else
      echo "Installation failed." $2
      exit 1;
  fi  
}

download_ext () {
  if [ -e $2 ]; then
      echo "$2 existed. Skiping downloading $1."
      return;
  fi; 
  
  if hash curl 2>/dev/null; then
      curl -L $1 -o $2; 
  else
      wget -O $2 $1; 
  fi; 

  if [ -z $2 ]; then
      echo "ERROR: $1 download failed."
  fi; 
}
 
if ! hash unzip 2>/dev/null; then
    echo "WARNNING: Unzip not found. You might experience error in further installation.";
fi;

echo "Checking BWA ..."

BWA_VER=`bwa 2>&1 | grep -Po "Version: \d+\.\d+" | grep -Po "\d+\.\d+"`;

if ( hash bwa 2>/dev/null ) && ( echo $BWA_VER | awk '{if($_>=0.7) exit 0; else exit 1}' )
then
  echo "BWA >=0.7 found.";
else
  echo "BWA >=0.7 or above not found. Trying to download from https://github.com/lh3/bwa/archive/master.zip ...";
  mkdir -p ext/opt;
  mkdir -p ext/bin;
  download_ext https://github.com/lh3/bwa/archive/master.zip ext/opt/bwa.zip;
  unzip ext/opt/bwa.zip -d ext/opt/;
  cd ext/opt/bwa-master;
  make;
  cd ../../../;
  cp ext/opt/bwa-master/bwa ext/bin/;
fi;
done_message " Done." "";

echo "Checking D Programming Language compiler ..."

if hash dmd 2>/dev/null; then
  echo "D compiler found."
else
  echo "D not found. Trying to download from http://downloads.dlang.org ..."
  mkdir -p ext/opt;
  mkdir -p ext/bin;

  #determine platform
  DMDPATH="ext/opt/dmd2"
  UNAME=`uname`
  UNAMEM=`uname -m`
  if echo $UNAME | grep -i "^linux$" > /dev/null; then
     DMDPATH="$DMDPATH/linux"
  elif echo $UNAME | grep -i "^freebsd$" > /dev/null; then
     DMDPATH="$DMDPATH/freebsd"
  fi; 

  if echo $UNAME | grep -i "^darwin$" > /dev/null; then
     DMDPATH="$DMDPATH/osx/bin/"
  elif echo $UNAMEM | grep -i "^x86_64$" > /dev/null; then
     DMDPATH="$DMDPATH/bin64/"
  elif echo $UNAMEM | grep -i "^i686$" > /dev/null; then
     DMDPATH="$DMDPATH/bin32/"
  elif echo $UNAMEM | grep -i "^i386$" > /dev/null; then
     DMDPATH="$DMDPATH/bin/"
  fi;

  if [ -x ${DMDPATH}dmd ]; then
      echo "DMD existed. Skiping downloading $1."
  else
      echo "Download DMD from http://downloads.dlang.org/..."
      download_ext http://downloads.dlang.org/releases/2014/dmd.2.065.0.zip ext/opt/dmd.zip;
      unzip -o ext/opt/dmd.zip -d ext/opt/;
      chmod a+x ${DMDPATH}dmd;
  fi; 
fi;
( set -xe;
  ${DMDPATH}dmd -O -release -op -inline src/splitrim.d;
  mv splitrim bin/;
)
done_message " Done compiling splitrim." "";

echo "Installing Perl dependencies..."

if hash cpanm 2>/dev/null; then
  echo "cpanm found. Start installing Perl dependencies..."
else
  echo "cpanm not found. Downloading from http://cpanmin.us..."
  download_ext http://cpanmin.us ext/bin/cpanm;
  chmod a+x ext/bin/cpanm;
fi
( set -xe;
  perl -Mthreads -e 1 > /dev/null 2>&1            || cpanm -v --notest -l ext threads;
  perl -Mthreads::shared -e 1 > /dev/null 2>&1    || cpanm -v --notest -l ext threads::shared;
  perl -MGetopt::Long -e 1 > /dev/null 2>&1       || cpanm -v --notest -l ext Getopt::Long;
  perl -MBenchmark -e 1 > /dev/null 2>&1          || cpanm -v --notest -l ext Benchmark;
  perl -MTime::HiRes -e 1 > /dev/null 2>&1        || cpanm -v --notest -l ext Time::HiRes;
  perl -MFile::Path -e 1 > /dev/null 2>&1         || cpanm -v --notest -l ext File::Path;
  perl -MFile::Basename -e 1 > /dev/null 2>&1     || cpanm -v --notest -l ext File::Basename;
  perl -MFile::Copy -e 1 > /dev/null 2>&1         || cpanm -v --notest -l ext File::Copy;
  perl -MIO::Handle -e 1 > /dev/null 2>&1         || cpanm -v --notest -l ext IO::Handle;
  perl -MYAML::XS -e 1 > /dev/null 2>&1           || cpanm -v --notest -l ext YAML::XS;
  perl -MYAML -e 1 > /dev/null 2>&1               || cpanm -v --notest -l ext YAML;
  perl -MXML::Simple -e 1 > /dev/null 2>&1        || cpanm -v --notest -l ext XML::Simple;
  perl -MStorable -e 1 > /dev/null 2>&1           || cpanm -v --notest -l ext Storable;
  perl -MStatistics::Descriptive -e 1 > /dev/null 2>&1 || cpanm -v --notest -l ext Statistics::Descriptive;
  perl -MTie::IxHash -e 1 > /dev/null 2>&1        || cpanm -v --notest -l ext Tie::IxHash;
  perl -MAlgorithm::Combinatorics -e 1 > /dev/null 2>&1 || cpanm -v --notest -l ext Algorithm::Combinatorics;
  perl -MDevel::Size -e 1 > /dev/null 2>&1        || cpanm -v --notest -l ext Devel::Size;
  perl -MSort::Key::Radix -e 1 > /dev/null 2>&1   || cpanm -v --notest -l ext Sort::Key::Radix;
  perl -MSort::Key -e 1 > /dev/null 2>&1          || cpanm -v --notest -l ext Sort::Key;
  perl -MBit::Vector -e 1 > /dev/null 2>&1        || cpanm -v --notest -l ext Bit::Vector;
  #perl -M"feature 'switch'" -e 1 > /dev/null 2>&1 || cpanm -v --notest -l ext feature;
)
done_message "Done installing Perl dependencies." "Failed installing Perl dependencies.";

echo -n "Moving scripts...";
( set -e;
  cp src/*.pl bin/
  chmod a+x bin/*
)
done_message "done." "Failed installing Perl dependencies.";

echo "
================================================================================
                 GOTTCHA installed successfully.
================================================================================

The pre-computed bacterial and viral GOTTCHA databases are available at our SFTP
server. If you want to create your own GOTTCHA signatue database, the detail
instruction can be found in README_FULL.txt file.

    SFTP server: img-gp.lanl.gov 
    Port: 33001 
    username: gottcha 
    password: 9001gottcha

Quick start:
    bin/gottcha.pl -i <FASTQ> -d <DATABASE>

Check our bitbucket site for update:
    https://bitbucket.org/poeli/gottcha

";
