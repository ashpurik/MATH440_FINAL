#!/bin/bash

make
mpiexec -machinefile xxx_machinefile -np 4 ./parallel 3
