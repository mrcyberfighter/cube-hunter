#!/bin/bash

#Application:  cube-hunter. 

#Dependencies: python-pygame
#Dependencies: python-opengl


sudo echo 'Start installation from the programm cube-hunter...'
if [[ $? != '0' ]] 
then
echo "This installation script must be run as root"
return 1
fi

# Check if module pygame is installed on your system.
if [[ ! -d '/usr/lib/python2.7/dist-packages/pygame/' ]] 
then 
echo 'WARNING: Cannot find python module pygame !!!'
echo 'Python module pygame required !!! For cube-hunter installation.'
echo 'so check if you have the python-pygame module installed on your system.'
echo 'If you have python-pygame installed on your system and this WARNING.'
echo 'Follow the the instruction to continue installallation script'
echo 'Else break the install.sh script by following the instructions.'
echo 'And install the python-pygame package. then retry after.'
echo
echo -n 'You have the python-pygame package installed ? [y|n]'
read get_pygame ;

if [[ ${get_pygame} != 'y'  ]] 
then
  return 1
fi
fi



# Check if module pyopengl is installed on your system.
if [[ ! -d '/usr/lib/python2.7/dist-packages/OpenGL/' ]] 
then 
echo 'WARNING: Cannot find python module pyopengl !!!'
echo
echo 'Python module pyopengl required !!! For cube-hunter installation.'
echo
echo 'so check if you have the python-opengl package installed on your system.'
echo 'If you have python-opengl package installed on your system ignore this WARNING.'
echo 'And Follow the the instructions to continue installallation script'

echo 'Else break the install.sh script by following the instructions.'
echo 'And install the python-opengl package. then retry after.'
echo
echo -n 'You have the python-opengl package installed ? [y|n]'
read get_pyopengl ;

if [[ ${get_pyopengl} != 'y'  ]] 
then
  return 1
fi
fi

if [[ ! -d "/usr/share/cube-hunter/" ]] 
then
  
  sudo mkdir "/usr/share/cube-hunter/" 
  
  sudo cp -R  "."     "/usr/share/cube-hunter/" 
  
fi


if [[ ! -f /usr/bin/cube_hunter.py ]] 
then 

  sudo cp  "$PWD/Source/cube_hunter.py" "/usr/bin/cube_hunter.py"
  
  sudo chmod a+x "/usr/bin/cube_hunter.py" 
  
fi




if [[ ! -f /usr/share/applications/cube-hunter.desktop ]] 
then

  sudo echo 'Start shortcut creation'
  
  sudo touch /usr/share/applications/cube-hunter.desktop

  sudo chmod a+rwx /usr/share/applications/cube-hunter.desktop

  sudo echo "[Desktop Entry]" > /usr/share/applications/cube-hunter.desktop
  sudo echo "Type=Application" >> /usr/share/applications/cube-hunter.desktop
  sudo echo "Encoding=UTF-8" >> /usr/share/applications/cube-hunter.desktop
  sudo echo "Name=cube-hunter" >> /usr/share/applications/cube-hunter.desktop
  sudo echo "GenericName=cube-hunter" >> /usr/share/applications/cube-hunter.desktop
  sudo echo "Comment=an shoot-them-all game in a cubic meshed sliding playground where you must eliminate the number of cubes by shooting them in relationship to your choosen settings." >> /usr/share/applications/cube-hunter.desktop
  sudo echo "Icon=/usr/share/cube-hunter/Icon/cube_hunter_icon.png" >> /usr/share/applications/cube-hunter.desktop
  sudo echo "Exec=python2 /usr/bin/cube_hunter.py" >> /usr/share/applications/cube-hunter.desktop
  sudo echo "Terminal=false" >> /usr/share/applications/cube-hunter.desktop
  sudo echo "StartupNotify=true" >> /usr/share/applications/cube-hunter.desktop
  sudo echo "Categories=Application;Game;" >> /usr/share/applications/cube-hunter.desktop
  sudo echo "" >> /usr/share/applications/cube-hunter.desktop

fi

#You can remove the decrompressed cube-hunter directory now.
echo 'cube-hunter installation successfull'