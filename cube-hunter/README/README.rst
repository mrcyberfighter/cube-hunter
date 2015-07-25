
cube-hunter
===========

  :Writing by: Eddie Brüggemann                                                 
  
  :Language: python 2.7.6                                   
  
  :Contact: mrcyberfighter@gmail.com                                            
  
  :Credits: Thank's to my mother and to the doctors.   


  .. image:: ./cube_hunter_icon.png

Description
-----------

  :cube-hunter: an shoot-them-all game in a cubic meshed sliding playground     
               
                where you must eliminate the number of cubes by shooting them   
               
                in relationship to your choosen settings. 

                                                                              
  **cube-hunter** is an game in an 3 dimensional cubic sliding playground where
  
  
  you have to destroy all the cubes which number you can configurate so the    
  
  maximal number of cubes present simultaneous in the playground and the       
  
  number of life you have which decrease wenn a cube hit you 5 times.          
  
  To mark the life losing the pentahedron representing the player will         
  
  explose.                                                                     
                                                                              
  
  You can navigate throught the playground by orientating the pentahedron and  
  
  reconize the direction you will slide to by the pike of the pentahedron.     
  
  You can shoot in the orientation pike direction wenn you be an mesh          
  
  intersection from the sliding cubic playground. You cannot shoot wenn you    
  
  slide the player pentahedron or wenn a shoot ball has not reach the          
  
  playground end or hit an cube.                                               
                                                                              
  
  The cubes to shoot can have from 1 to 6 life points, wich decrease every you 
  
  touch an cube with an shoot. The cubes navigate randomly throught the mesh   
  
  of the playground. The best tactic to explose a cube is to hunt it because   
  
  the player moving speed is 4 times faster than the cube moving speed.        
                                                                              
  
  The keys to move the player pentahedron are:                                 
                                                                              
    * **Right arrow** to turn the pentahedron to the right.                         
                                                                              
    * **Left arrow** to turn the pentahedron to the left.                           
                                                                              
    * **Up arrow** to turn the pentahedron to the top.                              
                                                                              
    * **Down arrow** to turn the pentahedron to the bottom.                         
                                                                              
  You can slide in the current direction with the:                             
                                                                              
    * **'z'** key for AZERTY keyboards.                                             
   
    * **'w'** key for QWERTY keyboards.                                             
                                                                              
  You can shoot in the pike orienation direction with the:                     
                                                                              
    * **'s'** key.                                                                  
                                                                              
  Wenn the player explose and you still have lifes press the                   
                                                                              
    * **'d'** key to reappear in the cube mesh center wenn you wish.                
                                                                              
  :NOTE: This keys combinations are the gaming keys so set your hands           
	 
	 right on the keyboard to have a better control on the game.            
                                                                              
  By the configuration screen use the:                                         
                                                                              
    * **Up and Down arrows** to change the selection                                
                                                                              
    * **Right arrow** to increase the value.                                        
                                                                              
    * **Left arrow** to decrease the value.                                         
                                                                              
  Wenn the game is running you can switch to pause by hitting the space bar    
  
  then you can orientate the meshed cube view as you wish with:                
                                                                              
    * The arrows, the **'d'** key and the **'a'** or **'q'** key in relationship if your    
      
      keyboard is an AZERTY or QWERTY.                                          
      
      And you can reset to the default view by hitting the **'z'** or **'w'** key.      
                                                                              
    * you can increase or decrease the stroke width of all components from the  
     
      playground with the keypad plus and minus keys.                           
                                                                              
  So have fun by hunting and shooting all the cubes down.                      
                                                                              


Weel-known bugs
---------------

                                                                              
  I must admit that is my first game with OpenGL and in 3 dimension and so, it
  
  has some unwanted behavior bugs like:                                        
                                                                              
  1. By shooting an cube:                                                      
                                                                              
    1. The shoot ball pass through an cube instead of hitting it.              
                                                                              
    2. The cube life points can be decrease more than one per hit.             
                                                                              
  2. By cube hit the player:                                                  
                                                                              
    1. The life points subcounter can decrease more than one point by an hit.  
                                                                              
    2. the collision cannot be taking in charge if a cube hit you.             
                                                                              
  
  That's all: I decide to publish the game despite this unwanted bugs          
  
  only because the game is funny to play an the mentionned bugs doesn't        
  
  disturb really the interest of the game in consideration of the              
  
  frequence they appears in one game.                                          
                                                                              
  
  And i hope if you are an python you can if you want try to fix it:           
  
  The code is commented and under GPLv3 license: you are free to modify it and 
  
  to redistribute it.                                                          
                                                                              
  
  The code is well-commented but a little bit trashy at my own taste.          
                                                                              


Installation
------------

  To install cune hunter.
 
  Simply run the script install.sh as root                                    
 
  ..
 
    $ su root
 
    password:
 
    # . install.sh                                                               
 
  You can remove the extracted folder from the zip archive.                    
                                                                              
 
  To uninstall cube-hunter from your system,                                   
 
  run the uninstall script uninstall.sh as root.                               
 
  ..
 
    $ su root
  
    password:
 
    # . uninstall
  
  This will remove all the installed files from cube-hunter from your system   
                                                                              


Copyright
---------


                                                                              
 **cube-hunter** an shoot-them-all game in a cubic meshed playground where the   
 
 pentahedron player have to eliminate all the cubes.             
 
 
 Copyright (C) 2014 Bruggemann Eddie                                          
                                                                              
 
 This file is part of **cube-hunter**.                                            
 
 **cube-hunter** is free software: you can redistribute it and/or modify          
 
 it under the terms of the GNU General Public License as published by         
 
 the Free Software Foundation, either version 3 of the License, or            
 
 (at your option) any later version.                                          
 
 
 **cube-hunter** is distributed in the hope that it will be useful,               
 
 but WITHOUT ANY WARRANTY; without even the implied warranty of               
 
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                 
 
 GNU General Public License for more details.                                 
                                                                              
 
 You should have received a copy of the GNU General Public License            
 
 along with cube-hunter. If not, see <http://www.gnu.org/licenses/>
