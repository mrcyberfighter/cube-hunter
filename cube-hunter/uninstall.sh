#!/bin/bash

sudo echo "Start desintallation of cube-hunter..."
if [[ $? != '0' ]] 
then
echo "This uninstall script must be run as root"
return 1
fi

sudo rm -R "/usr/share/cube-hunter/"
sudo rm "/usr/bin/cube_hunter.py"
sudo rm "/usr/share/applications/cube-hunter.desktop"

echo "cube-hunter successfull remove from your system !"