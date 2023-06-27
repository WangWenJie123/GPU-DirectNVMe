#!/bin/bash

# su passwd: #TCLOUDserver1411

if [ $1 == "unbind_nvme" ]
  then
  sudo sh -c 'sync && echo -n "0000:86:00.0" > \/sys\/bus\/pci\/devices\/0000\:86\:00.0\/driver\/unbind'
  command groups

  elif [ $1 == "bind_nvme" ]
  then
  sudo sh -c 'sync && echo -n "0000:86:00.0" > \/sys\/bus\/pci\/drivers\/nvme\/bind'
  command groups

  elif [ $1 == "bind_libnvm" ]
  then
  sudo sh -c 'sync && echo -n "0000:86:00.0" > \/sys\/bus\/pci\/drivers\/libnvm\ helper\/bind'
  command groups

  elif [ $1 == "unbind_libnvm" ]
  sudo sh -c 'sync && echo -n "0000:86:00.0" > \/sys\/bus\/pci\/drivers\/libnvm\ helper\/unbind'
  then
  command groups

else
  echo "Please input a operation!"
fi
