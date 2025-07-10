#!/usr/bin/env bash
#
# hrdp-variant-load-pipeline
#
. /etc/profile
APPNAME=hrdp-variant-load-pipeline
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

EMAILLIST="mtutaj@mcw.edu llamers@mcw.edu motutaj@mcw.edu"

APPDIR=/home/rgddata/pipelines/$APPNAME
cd $APPDIR

java -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$APPDIR/properties/log4j2.xml \
    -Xmx20g -jar lib/$APPNAME.jar --genicQc "$@" > run.log 2>&1

mailx -s "[$SERVER] HRDP Variant load Pipeline QC Run" llamers@mcw.edu < $APPDIR/logs/qcSummary.log
