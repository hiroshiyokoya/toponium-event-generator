#!/bin/sh
ln -sf thr11.inc thr.inc
touch matrix*.f
../bin/generate_events 2 8 Test_11
ln -sf thr01.inc thr.inc
touch matrix*.f
../bin/generate_events 2 8 Test_01
ln -sf thr10.inc thr.inc
touch matrix*.f
../bin/generate_events 2 8 Test_10
ln -sf thr00.inc thr.inc
touch matrix*.f
../bin/generate_events 2 8 Test_00
