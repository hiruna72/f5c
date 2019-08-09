#! /bin/sh
#
# Copyright by The HDF Group.
# Copyright by the Board of Trustees of the University of Illinois.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the COPYING file, which can be found at the root of the source code
# distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.
# If you do not have access to either file, you may request a copy from
# help@hdfgroup.org.
#
# Test to verify that the assertion/abort failure is fixed when the application
# does not close the file. (See HDFFV-10160)

srcdir=.

nerrors=0

##############################################################################
##############################################################################
###			  T H E   T E S T                                              ###
##############################################################################
##############################################################################

echo "Testing file not closed assertion/abort failure"
TEST_NAME=filenotclosed     # The test name
TEST_BIN=`pwd`/$TEST_NAME 	# The path of the test binary
#
# Run the test
#$RUNSERIAL $TEST_BIN >/dev/null 2>&1 
$RUNSERIAL $TEST_BIN 2>&1 
exitcode=$?
if [ $exitcode -eq 0 ]; then
    echo "Test PASSED"
else
 	nerrors="`expr $nerrors + 1`"
	echo "***Error encountered***"
fi
exit $nerrors
