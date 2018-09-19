#!/bin/bash

cd rel=1e-3
./run_test_performance
cd ../rel=1e-4
./run_test_performance
cd ../rel=1e-5
./run_test_performance
cd ../rel=1e-6
./run_test_performance

