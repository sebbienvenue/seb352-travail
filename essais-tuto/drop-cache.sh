#!/bin/bash
#a besoin du sudo
clear
sudo sync && echo 3 | sudo tee /proc/sys/vm/drop_caches
