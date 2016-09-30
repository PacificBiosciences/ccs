
Test if correctly stream to a log file

  $ $__PBTEST_CCS_EXE 2> cerr.log
  [1]
  $ if [ -s cerr.log ] ; then echo "dirty"; else echo "clean"; fi
  dirty

  $ $__PBTEST_CCS_EXE --logFile out.log 2> cerr.log
  [1]
  $ if [ -s cerr.log ] ; then echo "dirty"; else echo "clean"; fi
  clean
  $ if [ -s out.log ] ; then echo "dirty"; else echo "clean"; fi
  dirty
