#!/bin/bash

# run fastStructure for full matrix, K2
python /data/elinck/pabu/faststructure/fastStructure/structure.py -K 2 --input=/data/elinck/pabu/structure/pabu_c60h60 --output=pabu_fs_K2 --format=str

# run fastStructure for full matrix, K3
python /data/elinck/pabu/faststructure/fastStructure/structure.py -K 3 --input=/data/elinck/pabu/structure/pabu_c60h60 --output=pabu_fs_K3 --format=str

# run fastStructure for western pop only matrix, K2
python /data/elinck/pabu/faststructure/fastStructure/structure.py -K 2 --input=/data/elinck/pabu/faststructure/pabu_c60h60_westonly --output=pabu_fs_westonly_K2 --format=str

# run fastStructure for western pop only matrix, K3
python /data/elinck/pabu/faststructure/fastStructure/structure.py -K 3 --input=/data/elinck/pabu/faststructure/pabu_c60h60_westonly --output=pabu_fs_westonly_K3 --format=str
